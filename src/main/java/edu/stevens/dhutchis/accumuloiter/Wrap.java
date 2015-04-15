
package edu.stevens.dhutchis.accumuloiter;

import org.apache.accumulo.core.client.*;
import org.apache.accumulo.core.client.Scanner;
import org.apache.accumulo.core.data.Key;
import org.apache.accumulo.core.data.Mutation;
import org.apache.accumulo.core.data.Value;
import org.apache.accumulo.core.iterators.Combiner;
import org.apache.accumulo.core.iterators.IteratorUtil;
import org.apache.accumulo.core.iterators.LongCombiner;
import org.apache.accumulo.core.iterators.user.SummingCombiner;
import org.apache.accumulo.core.security.Authorizations;
import org.apache.hadoop.io.Text;
import org.junit.Assert;
import org.apache.accumulo.core.data.Range;

import java.util.*;
import java.util.Map.Entry;

/**
 *
 * @author dhutchis
 */
public class Wrap {

    public native boolean []    seqpass(String[] n , String hmm_path);

    private final String tableName = "TestTableIterator";
    private final String columnFamily="";
    private final String columnVisibility="";

    private void printList(Collection<?> list, String prefix) {
        System.out.print(prefix+": ");
        for (Object o : list) {
            System.out.print(o+", ");
        }
        System.out.println();
    }

    public static void main(String[] args)
    {
        System.out.println("HI im in main i updated again again");
	}

    void xin(Connector conn) throws AccumuloSecurityException, AccumuloException, TableNotFoundException
    {

        System.out.println("entasdfax to acc");
        // Setup BatchScanner to read rows that contain the accession numbers from TseqRaw, using 1 thread
        String TseqRaw = "TseqT";
        int numThreads = 1;
        Scanner scan = conn.createScanner(TseqRaw, Authorizations.EMPTY);
        scan.setRange(new Range("taxonomy" ,"taxonomy~"));

        Set<String> set = new HashSet<String>();

        // Do the scan

        for(Map.Entry<Key,Value> entry : scan)
        {

            String en = entry.getKey().toString();

            if(en.contains(";")) {
                set.add(en.substring(9, en.indexOf(';')));
            }

        }
        scan.close();

        for(String s : set)
        {
            System.out.println(s);
        }

    }

    //this is my new function

//the output should really be ranges instead of strings
    public List<Range> taxToAcc(Connector conn,List<String> taxaList) throws AccumuloSecurityException, AccumuloException, TableNotFoundException
    {
        System.out.println("entered tax to acc");
        // Setup BatchScanner to read rows that contain the accession numbers from TseqRaw, using 1 thread
        String TseqRaw = "TseqT";
        int numThreads = 1;
        Scanner scan = conn.createScanner(TseqRaw, Authorizations.EMPTY);
        scan.setRange(new Range("taxonomy|Bacteria; Cyanobacteria" ,"taxonomy|Bacteria; Cyanobacteria~"));

        List<Range> accList = new ArrayList<>();

        // Do the scan

        for(Map.Entry<Key,Value> entry : scan) {
                    String acc = entry.getKey().getColumnQualifier().toString();
                 // System.out.println("row is: " +entry.getKey().getRow().toString()+" acc is : " + acc );
                    accList.add(new Range(acc));
        }
        scan.close();

        return accList;
    }

    public String taxToRaw(Connector conn, String taxon) throws AccumuloSecurityException, AccumuloException, TableNotFoundException
    {
        long startTime = System.currentTimeMillis();

        System.load("/home/echerin/ppp/pipeline/src/main/java/edu/stevens/dhutchis/accumuloiter/Wrap.so");
        Wrap wrap = new Wrap();
        String hmm_path = "/home/echerin/48.hmm";


        int batchSize = 5000;
        System.out.println("entered tax to raw");
        // Setup BatchScanner to read rows that contain the accession numbers from TseqRaw, using 1 thread
        String TseqT = "TseqT";
        String TseqRaw = "TseqRaw";
        int numThreads = 1;

        Scanner scan = conn.createScanner(TseqT, Authorizations.EMPTY);
        scan.setRange(new Range(taxon ,taxon + "~"));

        BatchScanner batScan = conn.createBatchScanner(TseqRaw, Authorizations.EMPTY, numThreads);

        List<Range> accList = new ArrayList<>();
        Map<String,String> rawSeq = new HashMap<String,String>();

        // Do the scan
        int counter = 0;
        for(Map.Entry<Key,Value> entry : scan)
        {
            if(counter%batchSize == 0 && counter !=0)// when we have enough ranges
            {
                int co = 0;
                batScan.setRanges(accList);
                for(Map.Entry<Key,Value> batEntry : batScan)
                {
                    String seq = batEntry.getValue().toString();
                    String mykey = batEntry.getKey().toString();
                    rawSeq.put(mykey,seq);

                        System.out.println("counter is: : " + counter + " value is : " + seq + "\n");

                    co++;

                }

                wrap.seqpass( rawSeq.values().toArray(new String[rawSeq.size()]), hmm_path);
                accList  = new ArrayList<>();
            }

            String acc = entry.getKey().getColumnQualifier().toString();

            accList.add(new Range(acc));

            counter++;
        }


        //scan the last batch
        if(!accList.isEmpty())
        {
            batScan.setRanges(accList);
            for(Map.Entry<Key,Value> batEntry : batScan)
            {
                String seq = batEntry.getValue().toString();
                String mykey = batEntry.getKey().toString();
                rawSeq.put(mykey,seq);
            }
            wrap.seqpass( rawSeq.values().toArray(new String[rawSeq.size()]), hmm_path);
        }

        batScan.close();
        scan.close();
       String data =  Integer.toString(counter) + " " +  Long.toString(System.currentTimeMillis()-startTime);
        return data;
    }


    /** input is a list of Accession numbers ///used to be accession now its just con
     /   output is a map from accession numbers to sequences */
    public Map<String,String> accToRaw(Connector conn,List<Range> accessionList) throws AccumuloSecurityException, AccumuloException, TableNotFoundException
    {

        System.out.println("entered acc to raw");
        // Setup BatchScanner to read rows that contain the accession numbers from TseqRaw, using 1 thread
        String TseqRaw = "TseqRaw";
        int numThreads = 1;
        BatchScanner batScan = conn.createBatchScanner(TseqRaw, Authorizations.EMPTY, numThreads);

        /*
        List<Range> accessionRanges = new ArrayList<>(accessionList.size());
        for (String accession : accessionList)
        {
            accessionRanges.add(new Range(accession));
        }
        */


        Map<String,String> rawSeq = new HashMap<String,String>(accessionList.size());
        // Do the scan
        System.out.println("I am about to scan");

        batScan.setRanges(accessionList);
        for(Map.Entry<Key,Value> entry : batScan)
        {
            String seq = entry.getValue().toString();
            String mykey = entry.getKey().toString();
            //System.out.println("seq is: "  + seq + " key is : " + mykey);
            rawSeq.put(mykey,seq);
        }
        batScan.close();
        //rawSeq.values().toArray();
        return rawSeq;
    }


    public void testIter(Connector conn) throws AccumuloException, AccumuloSecurityException, TableNotFoundException, TableExistsException {

        System.out.println("I am in testIter");

        // check results
        Scanner scan = conn.createScanner("TseqDegT", Authorizations.EMPTY);
        scan.setRange(new Range("|Bacteria; Cyanobacteria"));
       // System.out.println("Scanner range: "+scan.getRange());
        int i = 0;
        for(Entry<Key,Value> entry : scan) {
           System.out.println(entry);
           // Assert.assertEquals("12", entry.getValue().toString());
            i++;
        }
        System.out.println("Number of results is: " + i );
        // try a custom iterator...

    }

    class ConcatIter extends Combiner {
        @Override
        public Value reduce(Key key, Iterator<Value> iter) {
            String ss = "";
            while (iter.hasNext()) {
                Value v = iter.next();
                ss += v.toString();
            }
            return new Value(ss.getBytes());
        }
    }
}

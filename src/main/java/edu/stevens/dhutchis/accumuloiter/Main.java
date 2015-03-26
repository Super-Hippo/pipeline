
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
public class Main {
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




    /** input is a list of Accession numbers
     /   output is a map from accession numbers to sequences */
    public Map<String,String> accToRaw(Connector conn,List<Range> accessionList) throws AccumuloSecurityException, AccumuloException, TableNotFoundException
    {

        System.out.println("entered acc to raw");
        // Setup BatchScanner to read rows that contain the accession numbers from TseqRaw, using 1 thread
        String TseqRaw = "TseqRaw";
        int numThreads = 1;
        BatchScanner scan = conn.createBatchScanner(TseqRaw, Authorizations.EMPTY, numThreads);


        scan.setRanges(accessionList);
System.out.println("size is: " + accessionList.size());

        Map<String,String> rawSeq = new HashMap<String,String>(accessionList.size());
        // Do the scan
        for(Map.Entry<Key,Value> entry : scan) {
            String seq = entry.getValue().toString();
            String mykey = entry.getKey().toString();
           //System.out.println("seq is: "  + seq + " key is : " + mykey);
            rawSeq.put(mykey,seq);
        }
        scan.close();
        //rawSeq.values().toArray();
        return rawSeq;
    }

    /*
    //old
        public List<String> accToRaw(Connector conn,List<String> accessionList) throws AccumuloSecurityException, AccumuloException, TableNotFoundException
    {

        System.out.println("entered acc to raw");
        // Setup BatchScanner to read rows that contain the accession numbers from TseqRaw, using 1 thread
        String TseqRaw = "TseqRaw";
        int numThreads = 1;
        BatchScanner scan = conn.createBatchScanner(TseqRaw, Authorizations.EMPTY, numThreads);

        List<Range> accessionRanges = new ArrayList<>(accessionList.size());
        for (String accession : accessionList)
        {
            accessionRanges.add(new Range(accession));
        }
        scan.setRanges(accessionRanges);

        List<String> rawSeq = new ArrayList<String>(accessionList.size());
        // Do the scan
        for(Map.Entry<Key,Value> entry : scan) {
            String seq = entry.getValue().toString();
           System.out.println(seq);
            rawSeq.add(seq);
        }
        scan.close();
        return rawSeq;
    }

     */

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

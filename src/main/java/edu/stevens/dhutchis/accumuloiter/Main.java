
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
	
    public void testIter(Connector conn) throws AccumuloException, AccumuloSecurityException, TableNotFoundException, TableExistsException {

        System.out.println("I am in testIter");

       // printList(conn.tableOperations().list(), "tables");


        /*Text rowID = new Text("row1");
        Text colFam = new Text("myColFam");
        Text colQual = new Text("myColQual");
        ColumnVisibility colVis = new ColumnVisibility();
        long timestamp = System.currentTimeMillis();
        Value value = new Value("myValue".getBytes());
        Mutation mutation = new Mutation(rowID);
        mutation.put(colFam, colQual, colVis, timestamp, value);*/
        /*
        if (conn.tableOperations().exists(tableName))
            conn.tableOperations().delete(tableName);
        conn.tableOperations().create(tableName);
        */

        String iterName = "summingIter";
		
        // Setup IteratorSetting
        IteratorSetting cfg = new IteratorSetting(1, iterName, SummingCombiner.class);
        LongCombiner.setEncodingType(cfg, LongCombiner.Type.STRING);
        // add columns to act on
        List<IteratorSetting.Column> combineColumns = new ArrayList<>();
        combineColumns.add(new IteratorSetting.Column(columnFamily, "leg"));
        Combiner.setColumns(cfg, combineColumns);

        // Add Iterator to table
        conn.tableOperations().attachIterator(tableName, cfg);
        // Verify successful add
        //Map<String,EnumSet<IteratorUtil.IteratorScope>> iterMap = conn.tableOperations().listIterators(tableName);
        //EnumSet<IteratorUtil.IteratorScope> iterScope = iterMap.get(iterName);
        //Assert.assertNotNull(iterScope);
        //Assert.assertTrue(iterScope.containsAll(EnumSet.allOf(IteratorUtil.IteratorScope.class)));
        
      //  Text row1 = new Text("row1");
       // Text cqleg = new Text("leg");
     //   Value[] vlegs = new Value[] {new Value("3".getBytes()), new Value("4".getBytes()), new Value("5".getBytes())  };

        // AH HA!! If you batch writes to the same key together, only the last one actually transmits to Accumulo.
        // Lesson: do each write separately.
//        Mutation m1 = new Mutation(row1);
//        for (Value vleg : vlegs)
//            m1.put(new Text(columnFamily), cqleg, vleg);
        
//        BatchWriterConfig config = new BatchWriterConfig();
  //      BatchWriter writer = conn.createBatchWriter(tableName, config);
//        writer.addMutation(m1);
//        writer.flush();

        /*
        for (Value vleg : vlegs) {
            Mutation m1 = new Mutation(row1);
            m1.put(new Text(columnFamily), cqleg, vleg);
            writer.addMutation(m1);
            writer.flush();
        }
        */


        // check results
        Scanner scan = conn.createScanner("Tseq", Authorizations.EMPTY);
        //scan.setRange(new Range("BAA22448.1"));
       // System.out.println("Scanner range: "+scan.getRange());
        for(Entry<Key,Value> entry : scan) {
           System.out.println(entry);
           // Assert.assertEquals("12", entry.getValue().toString());
        }

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

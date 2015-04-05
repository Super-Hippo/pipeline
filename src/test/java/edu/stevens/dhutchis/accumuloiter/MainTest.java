/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.stevens.dhutchis.accumuloiter;

import org.apache.accumulo.core.client.*;
import org.apache.accumulo.core.client.security.tokens.PasswordToken;
import org.apache.accumulo.core.data.Range;
import org.apache.accumulo.minicluster.MiniAccumuloCluster;
import org.junit.*;
import org.junit.rules.TemporaryFolder;

import java.io.File;
import java.io.*;
import java.util.*;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.Scanner;


public class MainTest {
    @Rule
    public TemporaryFolder tempFolder = new TemporaryFolder();

    private String username = "root";
    private String password = "secret";
    private static ClientConfiguration myconfig;
    static {
        String instance = "instance";
        String host = "b412srv.ece.stevens-tech.edu:2181";
        int timeout = 10000; // 10 seconds
        myconfig = ClientConfiguration.loadDefault().withInstance(instance).withZkHosts(host).withZkTimeout(timeout);
    }

    public Connector connectToAccumulo() throws AccumuloSecurityException, AccumuloException {
        Instance instance = new ZooKeeperInstance(myconfig);
        return instance.getConnector(username, new PasswordToken(password));
    }


    private void printList(Collection<?> list, String prefix) {
        System.out.print(prefix+": ");
        for (Object o : list) {
            System.out.print(o+", ");
        }
        System.out.println();
    }
    
    public MainTest() {
	
    }
    
    @BeforeClass
    public static void setUpClass() {
    }
    
    @AfterClass
    public static void tearDownClass() {
    }
    
    @Before
    public void setUp() {
    }
    
    @After
    public void tearDown() {
    }

    //@Test
    public void testMini() throws Exception {

        File tempDir = tempFolder.newFolder();
        MiniAccumuloCluster accumulo = new MiniAccumuloCluster(tempDir, "password");
        accumulo.start(); // doesn't work on Dylan's computer for some reason.  The OS closes the Zookeeper connection.

        Instance instance = new ZooKeeperInstance(accumulo.getInstanceName(), accumulo.getZooKeepers());
        Connector conn = instance.getConnector("root", new PasswordToken("password"));

      //  innerTest(instance, conn);

        accumulo.stop();
        tempDir.delete();
    }

    @Test
    public void testNormal() throws Exception {

        Connector conn = connectToAccumulo();
        System.out.println("hi");
        innerTest( conn);
    }

    private void innerTest( Connector conn) throws Exception {


        List<String> taxa = new ArrayList<>();
        taxa.add("Bacteria");
        taxa.add("Proteobacteria");
        String taxon = ""; //"taxonomy|Bacteria; Cyanobacteria";

        Main main = new Main();
        File f = new File("/home/echerin/ppp/pipeline/taxons.txt");
        java.util.Scanner s = new Scanner(f);
        List<String> data = new ArrayList<>();

        while(s.hasNextLine()) //assume taxonomy| ... has no spaces in family/genus names
        {

            taxon = s.nextLine();
            System.out.println("We are in outer while loop and using: " + taxon);

            java.util.Scanner wScan = new Scanner(taxon);
            String tInput = "taxonomy|";

            boolean first = true;
            while(wScan.hasNext())
            {
                if(first)
                {
                    tInput += wScan.next();
                    first = false;
                }
                else
                {
                    tInput += "; " + wScan.next();
                }
                System.out.println("scanning with : " + tInput);
                data.add(main.taxToRaw(conn, tInput));
            }

        }

        PrintWriter writer = new PrintWriter("graph" + new Date( ).toString() + ".txt", "UTF-8");

        for(String str : data)
        {
            writer.println(str);
        }
        writer.close();

/*
        String[]   s = result.values().toArray(new String[result.size()]);

        for(String str : s)
        {
            System.out.println(str);
        }
*/


        //Main main = new Main();
        //main.xin(conn);
    }
	

}

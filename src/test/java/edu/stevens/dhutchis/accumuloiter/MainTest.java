/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.stevens.dhutchis.accumuloiter;

import org.apache.accumulo.core.client.ClientConfiguration;
import org.apache.accumulo.core.client.Connector;
import org.apache.accumulo.core.client.Instance;
import org.apache.accumulo.core.client.ZooKeeperInstance;
import org.apache.accumulo.core.client.security.tokens.PasswordToken;
import org.apache.accumulo.minicluster.MiniAccumuloCluster;
import org.junit.*;
import org.junit.rules.TemporaryFolder;

import java.io.File;
import java.util.Collection;

/**
 *
 * @author dhutchis
 */
public class MainTest {
    @Rule
    public TemporaryFolder tempFolder = new TemporaryFolder();

    private String username = "root";
    private String password = "secret";
    private static ClientConfiguration myconfig;
    static {
        String instance = "instance";
        String host = "b412srv.ece.stevens-tech.edu:2181";
        int timeout = 100000;
        myconfig = new ClientConfiguration().withInstance(instance).withZkHosts(host).withZkTimeout(timeout);
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

    @Test
    public void testMini() throws Exception {

        File tempDir = tempFolder.newFolder();
        MiniAccumuloCluster accumulo = new MiniAccumuloCluster(tempDir, "password");
        accumulo.start(); // doesn't work on Dylan's computer for some reason.  The OS closes the Zookeeper connection.

        Instance instance = new ZooKeeperInstance(accumulo.getInstanceName(), accumulo.getZooKeepers());
        Connector conn = instance.getConnector("root", new PasswordToken("password"));

        innerTest(instance, conn);

        accumulo.stop();
        tempDir.delete();
    }

    @Test
    public void testNormal() throws Exception {
        Instance instance = new ZooKeeperInstance(myconfig.get(ClientConfiguration.ClientProperty.INSTANCE_NAME), myconfig.get(ClientConfiguration.ClientProperty.INSTANCE_ZK_HOST));
        //System.out.println("made instance : "+instance);
        Connector conn = instance.getConnector(username, new PasswordToken(password));
       // System.out.println("made connector: "+conn);

        innerTest(instance, conn);
    }

    private void innerTest(Instance instance, Connector conn) throws Exception {
        printList(conn.tableOperations().list(), "tables");
        printList(instance.getMasterLocations(), "master_locations");

        Main main = new Main();
        main.testIter(conn);
    }
	

}

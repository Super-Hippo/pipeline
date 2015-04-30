package edu.stevens.dhutchis.accumuloiter;

import org.apache.accumulo.core.data.ByteSequence;
import org.apache.accumulo.core.data.Key;
import org.apache.accumulo.core.data.Range;
import org.apache.accumulo.core.data.Value;
import org.apache.accumulo.core.iterators.IteratorEnvironment;
import org.apache.accumulo.core.iterators.SortedKeyValueIterator;
import org.apache.commons.lang.SerializationUtils;
import org.apache.log4j.Logger;

import java.io.File;
import java.io.IOException;
import java.util.*;

/**
 * Calls HMMER native library.
 * Takes in sequences from TseqRaw.
 * Batches into 5000 sequences each.
 * Calls HMMER, which returns an array of bool
 * Packages up results from all 5000 sequences.
 * Emits an entry with all those results.
 *
 * Todo: Timing with Watch
 */
public class HMMERIterator implements SortedKeyValueIterator<Key,Value> {
    private static final Logger log = Logger.getLogger(HMMERIterator.class);

    // Load native library
    static {
        // Check standard directories
//    List<File> directories = new ArrayList<File>(Arrays.asList(new File[]{new File("/usr/lib64"), new File("/usr/lib")}));
        // Check in ACCUMULO_HOME location, too
        String envAccumuloHome = System.getenv("ACCUMULO_HOME");
//    if (envAccumuloHome != null) {
//      directories.add(new File(envAccumuloHome + "/lib/native"));
//      directories.add(new File(envAccumuloHome + "/lib/native/map")); // old location, just in case somebody puts it here
//    }
        File libFile = new File(envAccumuloHome+"/lib/ext/Wrap.so");

        // Attempt to load from these directories, using standard names
        if (libFile.exists() && libFile.isFile()) {
            String errMsg = "Tried and failed to load native map library " + libFile;
            try {
                System.load(libFile.getAbsolutePath());
                log.info("Loaded native map shared library " + libFile);
            } catch (Exception | UnsatisfiedLinkError e) {
                log.error(errMsg, e);

                String ldLibraryPath = System.getProperty("java.library.path");
                errMsg = "Tried and failed to load native map library from " + ldLibraryPath;
                try {
                    System.loadLibrary("Wrap");
                    log.info("Loaded native map shared library from " + ldLibraryPath);
                } catch (Exception | UnsatisfiedLinkError e2) {
                    log.error(errMsg, e2);
                }

            }
        } else {
            log.debug("Native map library " + libFile + " not found or is not a file.");
        }
    }


    @SuppressWarnings("unchecked")
    protected Value hmmerAttachBool(String[] accIDs, String[] rawSeqs) {
        String hmm_path = "/home/echerin/48.hmm"; // TODO: put into init options
//    int batchSize = 5000; // TODO: make batch size option to init

        Wrap wrap = new Wrap();

        // TODO: at some point later, return the score/probability
        boolean[] booleans = wrap.seqpass(rawSeqs, hmm_path);

        HashMap<String,String> map = new HashMap<>(rawSeqs.length);
        for (int i = 0; i < rawSeqs.length; i++) {
            map.put(accIDs[i], (booleans[i] ? '1' : '0') + rawSeqs[i]);
        }

        byte[] bytes = SerializationUtils.serialize(map);

        return new Value(bytes);
    }


    @Override
    public void init(SortedKeyValueIterator<Key, Value> sortedKeyValueIterator, Map<String, String> map, IteratorEnvironment iteratorEnvironment) throws IOException {

    }

    @Override
    public boolean hasTop() {
        return false;
    }

    @Override
    public void next() throws IOException {

    }

    @Override
    public void seek(Range range, Collection<ByteSequence> collection, boolean b) throws IOException {

    }

    @Override
    public Key getTopKey() {
        return null;
    }

    @Override
    public Value getTopValue() {
        return null;
    }

    @Override
    public SortedKeyValueIterator<Key, Value> deepCopy(IteratorEnvironment iteratorEnvironment) {
        return null;
    }
}
package edu.stevens.dhutchis.accumuloiter;

import org.apache.accumulo.core.data.ByteSequence;
import org.apache.accumulo.core.data.Key;
import org.apache.accumulo.core.data.Range;
import org.apache.accumulo.core.data.Value;
import org.apache.accumulo.core.iterators.IteratorEnvironment;
import org.apache.accumulo.core.iterators.SortedKeyValueIterator;
import org.apache.commons.lang.SerializationUtils;
import org.apache.hadoop.io.Text;
import org.apache.log4j.Logger;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.concurrent.atomic.AtomicBoolean;

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

  private SortedKeyValueIterator<Key,Value> source;
  private String hmm_path = "/home/echerin/48.hmm";
  private Key topKey;
  private Value topValue;

//  private static final AtomicBoolean loadedNativeLibraries = new AtomicBoolean(false);

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




  //  @SuppressWarnings("unchecked")
  private Value hmmerAttachBool(String[] accIDs, String[] rawSeqs) {
    //    int batchSize = 5000; // TODO: make batch size option to init

    log.debug("hmmerAttachBool: rawSeqs.length= "+rawSeqs.length);
    // TODO: at some point later, return the score/probability
    boolean[] booleans = Wrap.seqpass(rawSeqs, hmm_path);

    HashMap<String,String> map = new HashMap<>(rawSeqs.length);
    for (int i = 0; i < rawSeqs.length; i++) {
      map.put(accIDs[i], (booleans[i] ? '1' : '0') + rawSeqs[i]);
    }

    byte[] bytes = SerializationUtils.serialize(map);

//    for (Map.Entry<Key, Value> entry : scanner) {
////        System.out.println("A Entry: "+entry.getKey() + " -> " + entry.getValue());
//      HashMap<String, String> map1 = (HashMap<String, String>) SerializationUtils.deserialize(entry.getValue().get());
//      for (Map.Entry<String, String> accToEncodedRawSeq : map1.entrySet()) {
//        String accID = accToEncodedRawSeq.getKey();
//        String tmp = accToEncodedRawSeq.getValue();
//        boolean b = tmp.charAt(0) != '0';
//        String rawSeq = tmp.substring(1);
//        // do something with accID, b, rawSeq
//      }
//    }

    return new Value(bytes);
  }

  private static final byte[] SEQ_COL = new Text("seq").copyBytes();

  private void prepareNextEntry() throws IOException {
    topKey = null;
    topValue = null;
    final int batchSize = Integer.MAX_VALUE;
    List<String> accIDs = new ArrayList<>(), rawSeqs = new ArrayList<>();
    Key k = new Key();
    for (int i = 0; i < batchSize && source.hasTop(); source.next(), i++) {
      // must have "seq" col
      if (!Arrays.equals(source.getTopKey().getColumnQualifierData().getBackingArray(), SEQ_COL))
        continue;
      k.set(source.getTopKey());
      accIDs.add(source.getTopKey().getRow().toString());
      rawSeqs.add(source.getTopValue().toString());
    }

    topKey = new Key(k);
    topValue = hmmerAttachBool(accIDs.toArray(new String[accIDs.size()]),
        rawSeqs.toArray(new String[rawSeqs.size()]));
  }


  @Override
  public void init(SortedKeyValueIterator<Key, Value> source, Map<String, String> map, IteratorEnvironment iteratorEnvironment) throws IOException {
    this.source = source;
    if (map != null && map.containsKey("hmm_path"))
      hmm_path = map.get("hmm_path");
    // TODO: add functionality from RemoteWriteIterator to control seeking
  }

  @Override
  public boolean hasTop() {
    return topKey != null;
  }

  @Override
  public void next() throws IOException {
    prepareNextEntry();
  }

  @Override
  public void seek(Range range, Collection<ByteSequence> columnFamilies, boolean inclusive) throws IOException {
//    this.seekRange = range;
    source.seek(range, columnFamilies, inclusive);
    prepareNextEntry();
  }

  @Override
  public Key getTopKey() {
    return topKey;
  }

  @Override
  public Value getTopValue() {
    return topValue;
  }

  @Override
  public SortedKeyValueIterator<Key, Value> deepCopy(IteratorEnvironment iteratorEnvironment) {
    HMMERIterator other = new HMMERIterator();
    other.hmm_path = hmm_path;
    other.source = source.deepCopy(iteratorEnvironment);
    return other;
  }
}

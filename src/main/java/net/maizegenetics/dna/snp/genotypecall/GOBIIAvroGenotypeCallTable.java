/*
 *  GOBIIAvroGenotypeCallTable
 * 
 *  Created on September 21, 2016
 */
package net.maizegenetics.dna.snp.genotypecall;

import com.google.common.cache.Cache;
import com.google.common.cache.CacheBuilder;
import java.nio.ByteBuffer;
import java.util.WeakHashMap;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ForkJoinPool;
import net.maizegenetics.dna.snp.NucleotideAlignmentConstants;
import net.maizegenetics.util.Tuple;
import org.apache.avro.generic.GenericData;
import org.apache.avro.generic.GenericRecord;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 *
 * @author Terry Casstevens
 *
 * This class implements the GenotypeCallTable interface for the specific
 * purpose to cache genotypes read from the GOBII Avro files (Resulting from the
 * genotypes stored in the MonetDB).
 *
 */
public class GOBIIAvroGenotypeCallTable extends AbstractGenotypeCallTable {

    private static final Logger myLogger = LogManager.getLogger(GOBIIAvroGenotypeCallTable.class);
    private static final int NUM_LOOK_AHEAD_BLOCKS = 5;
    public static final int GENOTYPE_BLOCK_SIZE = 256;

    private final GenericRecord myRecord;
    private final WeakHashMap<Thread, Tuple<Long, byte[][]>> myLastSite = new WeakHashMap<>();

    private final Cache<Long, byte[][]> myGenoCache;

    private final ConcurrentHashMap<Long, CompletableFuture<byte[][]>> myFutureQueue = new ConcurrentHashMap<>();

    private final ForkJoinPool myThreadPool;

    private GOBIIAvroGenotypeCallTable(int numTaxa, int numSites, boolean phased, GenericRecord record) {
        super(numTaxa, numSites, phased, NucleotideAlignmentConstants.NUCLEOTIDE_ALLELES);

        long oneThirdMemory = Runtime.getRuntime().maxMemory() / 65536;
        int maxCacheSize = (int) Math.min((long) (100 * Runtime.getRuntime().availableProcessors()), oneThirdMemory);

        myGenoCache = CacheBuilder.newBuilder()
                .initialCapacity(maxCacheSize)
                .maximumSize(maxCacheSize)
                .build();
        myThreadPool = new ForkJoinPool();

        myRecord = record;
    }

    public static GOBIIAvroGenotypeCallTable getInstance(int numTaxa, int numSites, boolean phased, GenericRecord record) {
        return new GOBIIAvroGenotypeCallTable(numTaxa, numSites, phased, record);
    }

    public static long getCacheKey(int taxon, int site) {
        return ((long) (taxon / GENOTYPE_BLOCK_SIZE) << 32) + (site / GENOTYPE_BLOCK_SIZE);
    }

    public static String getKey(int taxon, int site) {
        return "B" + String.valueOf(getCacheKey(taxon, site));
    }

    public static int myNumCacheMisses = 0;

    private byte[][] getFromCache(int taxon, int site) {

        //System.out.println("start getFromCache: taxon: " + taxon + "  site: " + site);
        long cacheKey = getCacheKey(taxon, site);

        byte[][] result = myGenoCache.getIfPresent(cacheKey);

        if (result == null) {
            myNumCacheMisses++;

            CompletableFuture<byte[][]> future = new CompletableFuture<>();
            CompletableFuture<byte[][]> temp = myFutureQueue.putIfAbsent(cacheKey, future);
            try {
                if (temp == null) {
                    myThreadPool.submit(new ReadBlocks(taxon, site));
                    return future.get();
                } else {
                    return temp.get();
                }

            } catch (Exception e) {
                myLogger.error(e.getMessage(), e);
            }

        }

        //System.out.println("end getFromCache: taxon: " + taxon + "  site: " + site);
        return result;

    }

    @Override
    public byte genotype(int taxon, int site) {
        try {
            Tuple<Long, byte[][]> temp = myLastSite.get(Thread.currentThread());
            long key = getCacheKey(taxon, site);
            if (temp != null) {
                if (temp.x == key) {
                    return temp.y[site % GENOTYPE_BLOCK_SIZE][taxon % GENOTYPE_BLOCK_SIZE];
                } else {
                    byte[][] result = getFromCache(taxon, site);
                    myLastSite.put(Thread.currentThread(), new Tuple<>(key, result));
                    return result[site % GENOTYPE_BLOCK_SIZE][taxon % GENOTYPE_BLOCK_SIZE];
                }
            } else {
                byte[][] result = getFromCache(taxon, site);
                myLastSite.put(Thread.currentThread(), new Tuple<>(key, result));
                return result[site % GENOTYPE_BLOCK_SIZE][taxon % GENOTYPE_BLOCK_SIZE];
            }
        } catch (Exception ex) {
            myLogger.error(ex.getMessage(), ex);
            throw new IllegalStateException("GOBIIGenotypeCallTable: genotype: Error getting genotype from cache: " + ex.getMessage());
        }
    }

    @Override
    public byte[] genotypeForAllTaxa(int site) {

        byte[] result = new byte[myTaxaCount];

        int count = 0;
        int siteOffset = site % GENOTYPE_BLOCK_SIZE;
        for (int t = 0; t < myTaxaCount; t += GENOTYPE_BLOCK_SIZE) {
            byte[][] temp = getFromCache(t, site);
            System.arraycopy(temp[siteOffset], 0, result, GENOTYPE_BLOCK_SIZE * count, temp[siteOffset].length);
            count++;
        }

        return result;

    }

    @Override
    public byte[] genotypeAllSites(int taxon) {

        byte[] result = new byte[mySiteCount];

        int count = 0;
        int taxaOffset = taxon % GENOTYPE_BLOCK_SIZE;
        for (int s = 0; s < mySiteCount; s += GENOTYPE_BLOCK_SIZE) {
            byte[][] temp = getFromCache(taxon, s);
            for (int s1 = 0, n = temp.length; s1 < n; s1++) {
                result[count++] = temp[s1][taxaOffset];
            }
        }

        return result;

    }

    @Override
    public String genotypeAsString(int taxon, int site) {
        return NucleotideAlignmentConstants.getNucleotideIUPAC(genotype(taxon, site));
    }

    @Override
    public String diploidAsString(int site, byte value) {
        return NucleotideAlignmentConstants.getNucleotideIUPAC(value);
    }

    @Override
    public void transposeData(boolean siteInnerLoop) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public boolean isSiteOptimized() {
        return true;
    }

    public static int myNumBlocksRead = 0;
    public static long myTimeReading = 0;
    public static int myNumProcesses = 0;

    private class ReadBlocks implements Runnable {

        private final int myStartTaxa;
        private final int myStartSite;
        private long myProcessBlock;

        public ReadBlocks(int taxon, int site) {
            myProcessBlock = getCacheKey(taxon, site);
            myStartTaxa = taxon / GENOTYPE_BLOCK_SIZE * GENOTYPE_BLOCK_SIZE;
            myStartSite = site / GENOTYPE_BLOCK_SIZE * GENOTYPE_BLOCK_SIZE;
        }

        @Override
        public void run() {
            //System.out.println("starting ReadBlocks: " + myProcessBlock);

            myNumProcesses++;

            if ((myStartTaxa >= myTaxaCount) || (myStartSite >= mySiteCount)) {
                return;
            }

            try {

                //GenericRecord record = null;
                byte[][] temp = myGenoCache.getIfPresent(myProcessBlock);
                if (temp != null) {
                    CompletableFuture<byte[][]> future = myFutureQueue.get(myProcessBlock);
                    if (future != null) {
                        future.complete(temp);
                        myFutureQueue.remove(myProcessBlock, future);
                    }
                } else {

                    myNumBlocksRead++;
                    //long previous = System.nanoTime();
                    GenericData.Array<ByteBuffer> byteBufferArray = (GenericData.Array<ByteBuffer>) myRecord.get(getKey(myStartTaxa, myStartSite));
                    if (byteBufferArray == null) {
                        throw new IllegalStateException("GOBIIAvroGenotypeCallTable: byte buffer array is null: " + getKey(myStartTaxa, myStartSite));
                    }
                    int numTaxa = byteBufferArray.size();
                    ByteBuffer current = byteBufferArray.get(0);
                    current.rewind();
                    int numSites = current.remaining();
                    byte[][] result = new byte[numTaxa][numSites];
                    for (int t = 0; t < numTaxa; t++) {
                        current = byteBufferArray.get(t);
                        current.rewind();
                        for (int s = 0; s < numSites; s++) {
                            result[t][s] = current.get();
                        }
                    }

                    myGenoCache.put(myProcessBlock, result);

                    CompletableFuture<byte[][]> future = myFutureQueue.get(myProcessBlock);
                    if (future != null) {
                        future.complete(result);
                        myFutureQueue.remove(myProcessBlock, future);
                    }

                    // This get to prevent early eviction from cache
                    myGenoCache.getIfPresent(myProcessBlock);
                    //myCurrentlyProcessingBlocks.remove(myProcessBlock);

                }

                int lookaheadTaxa = myStartTaxa;
                int lookaheadSite = myStartSite;

                for (int b = 0; b < NUM_LOOK_AHEAD_BLOCKS; b++) {

                    lookaheadTaxa += GENOTYPE_BLOCK_SIZE;
                    myProcessBlock = getCacheKey(lookaheadTaxa, myStartSite);
                    //System.out.println("look ahead ReadBlocks: " + myProcessBlock);
                    if ((lookaheadTaxa < myTaxaCount)
                            && (myGenoCache.getIfPresent(myProcessBlock) == null)
                            && (myFutureQueue.putIfAbsent(myProcessBlock, new CompletableFuture<>()) == null)) {

                        //System.out.println("look ahead inside check");
                        myNumBlocksRead++;
                        //previous = System.nanoTime();
                        //System.out.println("look ahead after record");
                        GenericData.Array<ByteBuffer> byteBufferArray = (GenericData.Array<ByteBuffer>) myRecord.get(getKey(lookaheadTaxa, myStartSite));
                        if (byteBufferArray == null) {
                            throw new IllegalStateException("GOBIIAvroGenotypeCallTable: byte buffer array is null: " + getKey(lookaheadTaxa, myStartSite));
                        }
                        int numTaxa = byteBufferArray.size();
                        ByteBuffer current = byteBufferArray.get(0);
                        current.rewind();
                        int numSites = current.remaining();
                        byte[][] result = new byte[numTaxa][numSites];
                        for (int t = 0; t < numTaxa; t++) {
                            current = byteBufferArray.get(t);
                            current.rewind();
                            for (int s = 0; s < numSites; s++) {
                                result[t][s] = current.get();
                            }
                        }
                        //myTimeReading += System.nanoTime() - previous;

                        myGenoCache.put(myProcessBlock, result);

                        CompletableFuture<byte[][]> future = myFutureQueue.get(myProcessBlock);
                        if (future != null) {
                            future.complete(result);
                            myFutureQueue.remove(myProcessBlock, future);
                        }

                        // This get to prevent early eviction from cache
                        myGenoCache.getIfPresent(myProcessBlock);
                        //myCurrentlyProcessingBlocks.remove(myProcessBlock);

                    }

                    lookaheadSite += GENOTYPE_BLOCK_SIZE;
                    myProcessBlock = getCacheKey(myStartTaxa, lookaheadSite);
                    //System.out.println("look ahead 2 ReadBlocks: " + myProcessBlock);
                    if ((lookaheadSite < mySiteCount)
                            && (myGenoCache.getIfPresent(myProcessBlock) == null)
                            && (myFutureQueue.putIfAbsent(myProcessBlock, new CompletableFuture<>()) == null)) {

                        //System.out.println("look ahead 2 inside check");
                        myNumBlocksRead++;
                        //previous = System.nanoTime();
                        GenericData.Array<ByteBuffer> byteBufferArray = (GenericData.Array<ByteBuffer>) myRecord.get(getKey(myStartTaxa, lookaheadSite));
                        if (byteBufferArray == null) {
                            throw new IllegalStateException("GOBIIAvroGenotypeCallTable: byte buffer array is null: " + getKey(myStartTaxa, lookaheadSite));
                        }
                        int numTaxa = byteBufferArray.size();
                        ByteBuffer current = byteBufferArray.get(0);
                        current.rewind();
                        int numSites = current.remaining();
                        byte[][] result = new byte[numTaxa][numSites];
                        for (int t = 0; t < numTaxa; t++) {
                            current = byteBufferArray.get(t);
                            current.rewind();
                            for (int s = 0; s < numSites; s++) {
                                result[t][s] = current.get();
                            }
                        }
                        //myTimeReading += System.nanoTime() - previous;

                        myGenoCache.put(myProcessBlock, result);

                        CompletableFuture<byte[][]> future = myFutureQueue.get(myProcessBlock);
                        if (future != null) {
                            future.complete(result);
                            myFutureQueue.remove(myProcessBlock, future);
                        }

                        // This get to prevent early eviction from cache
                        myGenoCache.getIfPresent(myProcessBlock);
                        //myCurrentlyProcessingBlocks.remove(myProcessBlock);

                    }

                }

            } catch (Exception e) {
                myLogger.error(e.getMessage(), e);
            }

            //System.out.println("ending ReadBlocks: " + myProcessBlock);
        }

    }

}

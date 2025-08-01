/*
 *  GOBIIGenotypeCallTable
 * 
 *  Created on July 11, 2016
 */
package net.maizegenetics.dna.snp.genotypecall;

import ch.systemsx.cisd.base.mdarray.MDArray;
import ch.systemsx.cisd.hdf5.HDF5Factory;
import ch.systemsx.cisd.hdf5.IHDF5Reader;
import com.google.common.cache.Cache;
import com.google.common.cache.CacheBuilder;
import java.util.WeakHashMap;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.CopyOnWriteArraySet;
import java.util.concurrent.ForkJoinPool;
import net.maizegenetics.dna.snp.NucleotideAlignmentConstants;
import net.maizegenetics.util.Tuple;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 *
 * @author Terry Casstevens
 *
 * This class implements the GenotypeCallTable interface for the specific
 * purpose to cache genotypes read from the GOBII HDF5 files (Resulting from the
 * genotypes stored in the MonetDB).
 *
 */
public class GOBIIGenotypeCallTable extends AbstractGenotypeCallTable {

    private static final Logger myLogger = LogManager.getLogger(GOBIIGenotypeCallTable.class);
    private static final int NUM_LOOK_AHEAD_BLOCKS = 103;

    private final String myFilename;
    private final int myNumLinesPerInterval = 100;
    private final int myMaxCacheSize;
    private final ConcurrentLinkedQueue<IHDF5Reader> myReaders = new ConcurrentLinkedQueue<>();
    private final CopyOnWriteArraySet<Integer> myCurrentlyProcessingBlocks = new CopyOnWriteArraySet<>();
    private final WeakHashMap<Thread, Tuple<Integer, byte[]>> myLastSite = new WeakHashMap<>();

    private final Cache<Integer, byte[][]> myGenoCache;

    private final ConcurrentHashMap<Integer, CompletableFuture<byte[]>> myFutureQueue = new ConcurrentHashMap<>();

    private final ForkJoinPool myThreadPool;

    private GOBIIGenotypeCallTable(int numTaxa, int numSites, boolean phased, String filename) {
        super(numTaxa, numSites, phased, NucleotideAlignmentConstants.NUCLEOTIDE_ALLELES);
        myFilename = filename;

        long oneThirdMemory = Runtime.getRuntime().maxMemory() / (numTaxa * 3);
        myMaxCacheSize = (int) Math.min((long) (110 * Runtime.getRuntime().availableProcessors()), oneThirdMemory);

        myGenoCache = CacheBuilder.newBuilder()
                .initialCapacity(myMaxCacheSize)
                .maximumSize(myMaxCacheSize)
                .build();
        myThreadPool = new ForkJoinPool();
    }

    public static GOBIIGenotypeCallTable getInstance(int numTaxa, int numSites, boolean phased, String filename) {
        return new GOBIIGenotypeCallTable(numTaxa, numSites, phased, filename);
    }

    private byte[] getFromCache(int site) {

        int blockNumber = site / myNumLinesPerInterval;

        byte[][] result = myGenoCache.getIfPresent(blockNumber);

        if (result == null) {

            CompletableFuture<byte[]> future = new CompletableFuture<>();
            CompletableFuture<byte[]> temp = myFutureQueue.putIfAbsent(site, future);
            if (temp != null) {
                future = temp;
            }
            if (myCurrentlyProcessingBlocks.add(blockNumber)) {
                myThreadPool.submit(new ProcessLines(site));
            }

            try {
                result = myGenoCache.getIfPresent(blockNumber);
                if (result != null) {
                    myFutureQueue.remove(site);
                    future.complete(result[site % myNumLinesPerInterval]);
                    return result[site % myNumLinesPerInterval];
                } else {
                    return future.get();
                }
            } catch (Exception e) {
                myLogger.error(e.getMessage(), e);
            }

        }

        return result[site % myNumLinesPerInterval];

    }

    private IHDF5Reader getReader() {
        IHDF5Reader reader = myReaders.poll();
        if (reader == null) {
            try {
                reader = HDF5Factory.openForReading(myFilename);
            } catch (Exception e) {
                myLogger.error(e.getMessage(), e);
            }
        }
        return reader;
    }

    @Override
    public byte genotype(int taxon, int site) {
        try {
            Tuple<Integer, byte[]> temp = myLastSite.get(Thread.currentThread());
            if (temp != null) {
                if (temp.x == site) {
                    return temp.y[taxon];
                } else {
                    byte[] result = getFromCache(site);
                    myLastSite.put(Thread.currentThread(), new Tuple<>(site, result));
                    return result[taxon];
                }
            } else {
                byte[] result = getFromCache(site);
                myLastSite.put(Thread.currentThread(), new Tuple<>(site, result));
                return result[taxon];
            }
        } catch (Exception ex) {
            myLogger.error(ex.getMessage(), ex);
            throw new IllegalStateException("GOBIIGenotypeCallTable: genotype: Error getting genotype from cache: " + ex.getMessage());
        }
    }

    @Override
    public byte[] genotypeForAllTaxa(int site) {
        byte[] result = new byte[myTaxaCount];
        System.arraycopy(getFromCache(site), 0, result, 0, myTaxaCount);
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

    /**
     * Parse line from GOBII HDF5 file to genotypes for a site.
     *
     * @param input input line
     * @param numTaxa number of taxa
     * @param site site
     * @param relativeSite relative site in input
     *
     * @return genotypes
     */
    private static byte[] parseLine(MDArray<String> input, int numTaxa, int site, int relativeSite) {

        if (input.dimensions()[1] != numTaxa) {
            throw new IllegalStateException("GOBIIGenotypeCallTable: Site: " + site + " has wrong number of taxa: " + input.size() + ".  Number of taxa: " + numTaxa);
        }

        byte[] data = new byte[numTaxa];
        for (int t = 0; t < numTaxa; t++) {
            data[t] = NucleotideAlignmentConstants.getNucleotideDiploidByte(input.get(relativeSite, t));
        }

        return data;
    }

    private class ProcessLines implements Runnable {

        private int myStartSite;
        private final int myProcessBlock;

        public ProcessLines(int site) {
            myProcessBlock = site / myNumLinesPerInterval;
            myStartSite = myProcessBlock * myNumLinesPerInterval;
        }

        @Override
        public void run() {

            if (myStartSite >= mySiteCount) {
                return;
            }

            IHDF5Reader reader = getReader();
            try {

                int numSites = Math.min(myNumLinesPerInterval, mySiteCount - myStartSite);
                byte[][] result = new byte[numSites][];
                MDArray<String> input = reader.string().readMDArrayBlockWithOffset("allelematrix", new int[]{numSites, myTaxaCount}, new long[]{myStartSite, 0});
                for (int i = 0; i < numSites; i++) {
                    result[i] = parseLine(input, myTaxaCount, myStartSite + i, i);
                    CompletableFuture<byte[]> future = myFutureQueue.remove(myStartSite + i);
                    if (future != null) {
                        future.complete(result[i]);
                    }
                }
                myGenoCache.put(myProcessBlock, result);
                // This get to prevent early eviction from cache
                myGenoCache.getIfPresent(myProcessBlock);
                myCurrentlyProcessingBlocks.remove(myProcessBlock);
                for (int i = 0; i < numSites; i++) {
                    CompletableFuture<byte[]> future = myFutureQueue.remove(myStartSite + i);
                    if (future != null) {
                        future.complete(result[i]);
                    }
                }
                myStartSite += myNumLinesPerInterval;
                if (myStartSite >= mySiteCount) {
                    return;
                }

                for (int b = 1; b < NUM_LOOK_AHEAD_BLOCKS; b++) {

                    if (myGenoCache.getIfPresent(myProcessBlock + b) != null) {
                        return;
                    }
                    if (!myCurrentlyProcessingBlocks.add(myProcessBlock + b)) {
                        return;
                    }

                    numSites = Math.min(myNumLinesPerInterval, mySiteCount - myStartSite);
                    result = new byte[numSites][];
                    input = reader.string().readMDArrayBlockWithOffset("allelematrix", new int[]{numSites, myTaxaCount}, new long[]{myStartSite, 0});
                    for (int i = 0; i < numSites; i++) {
                        result[i] = parseLine(input, myTaxaCount, myStartSite + i, i);
                    }
                    myGenoCache.put(myProcessBlock + b, result);
                    // This get to prevent early eviction from cache
                    myGenoCache.getIfPresent(myProcessBlock + b);
                    myCurrentlyProcessingBlocks.remove(myProcessBlock + b);
                    for (int i = 0; i < numSites; i++) {
                        CompletableFuture<byte[]> future = myFutureQueue.remove(myStartSite + i);
                        if (future != null) {
                            future.complete(result[i]);
                        }
                    }
                    myStartSite += myNumLinesPerInterval;
                    if (myStartSite >= mySiteCount) {
                        return;
                    }
                }

            } catch (Exception e) {
                myLogger.error(e.getMessage(), e);
            } finally {
                myReaders.add(reader);
            }

        }

    }

    public static void main(String[] args) {
        String filename = "/SSD/gobii/gobii_terry/DS_1.h5";
        IHDF5Reader reader = HDF5Factory.openForReading(filename);
        MDArray<String> input = reader.string().readMDArrayBlockWithOffset("allelematrix", new int[]{5, 282}, new long[]{1, 0});

        int[] dimensions = input.dimensions();
        for (int dim : dimensions) {
            System.out.println("dim: " + dim);
        }

        for (int i = 0; i < input.size(); i++) {
            System.out.println(i + ": " + input.get(i));
        }
    }

}

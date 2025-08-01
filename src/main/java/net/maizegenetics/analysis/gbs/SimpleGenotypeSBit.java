/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.analysis.gbs;
import ch.systemsx.cisd.hdf5.*;
import net.maizegenetics.dna.WHICH_ALLELE;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.ImportUtils;
import net.maizegenetics.util.OpenBitSet;

import java.io.*;
import java.util.Arrays;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

/**
 * Store sBit for tag genetic mapping and LD detection, optimized for really fast calculation
 * @author Fei Lu
 */
public class SimpleGenotypeSBit {
    public int taxaNum;
    public int chrNum;
    public int siteNum;
    public int wordNum;
    
    public String[] taxaNames;
    public int[] chromosomeNumber;
    /**First SNP index of gt chromosome in the whole SNP list */
    public int[] chrStartIndex;
    /**Last SNP index of gt chromosome in the whole SNP list, exclusive */
    public int[] chrEndIndex;
    
    public int[] position;
    public double[] maf;
    public OpenBitSet[] obsMajor;
    public OpenBitSet[] obsMinor;
    
    /**HDF5 setting*/
    private static final int BITS_TO_SHIFT_FOR_CHUNK = 8;
    private static final int CHUNK_SIZE = 1 << BITS_TO_SHIFT_FOR_CHUNK;
    int maxTaxaNameLength = 100;
    /**Attributes*/
    String ROOT = "/";
    String TAXANUM = "taxaNum";
    String CHRNUM = "chrNum";
    String SITENUM = "siteNum";
    String WORDNUM = "wordNum";
    /**Path*/
    String TAXANAMES = "taxaNames";
    String CHROMOSOME = "chromosome";
    String CHRSTARTINDEX = "chrStartIndex";
    String CHRENDINDEX = "chrEndIndex";
    String POSITION = "position";
    String MAF = "maf";
    String OBSMAJOR = "obsMajor";
    String OBSMINOR = "obsMinor";
    /**compression level*/
    HDF5IntStorageFeatures intFeature = HDF5IntStorageFeatures.createDeflation(1);
    HDF5GenericStorageFeatures genericFeature = HDF5GenericStorageFeatures.createDeflation(1);
    HDF5FloatStorageFeatures floatFeature = HDF5FloatStorageFeatures.createDeflation(1);
    
    /**
     * Convert HDF5 Alignment/Genotype file to SimpleGenotypeSBit
     * @param genotypeH5FileS
     * @param sBitFileS 
     */
    public SimpleGenotypeSBit (String genotypeH5FileS, String sBitFileS) {
        long lastTimePoint  = System.nanoTime();
        GenotypeTable gt = ImportUtils.readGuessFormat(genotypeH5FileS);
        taxaNames = new String[gt.numberOfTaxa()];
        for (int i = 0; i < taxaNames.length; i++) {
            taxaNames[i] = gt.taxaName(i);
        }
        taxaNum = taxaNames.length;
        int[] chrOffSet = gt.chromosomesOffsets();
        chromosomeNumber = new int[gt.numChromosomes()];
        chrStartIndex = new int[chromosomeNumber.length];
        chrEndIndex = new int[chromosomeNumber.length];
        for (int i = 0; i < chromosomeNumber.length; i++) {
            chromosomeNumber[i] = gt.chromosomes()[i].getChromosomeNumber();
            chrStartIndex[i] = chrOffSet[i];
            chrEndIndex[i] = chrOffSet[i] + gt.chromosomeSiteCount(gt.chromosomes()[i]);
        }
        chrNum = this.chromosomeNumber.length;
        position = gt.physicalPositions();
        siteNum = this.position.length;
        System.out.println("This genotype has " + taxaNum + " taxa, " + chromosomeNumber.length + " chromosomes, " + siteNum + " sites");
        System.out.println("Will be transformed to sBit with " + this.getChunkNum() + " chunks, each chunk has " + this.getChunkSize()+" sites");
        maf = new double[gt.numberOfSites()];
        obsMajor = new OpenBitSet[gt.numberOfSites()];
        obsMinor = new OpenBitSet[gt.numberOfSites()];    
        OpenBitSet bsMa, bsMi, het;
        long majorCount, minorCount;
        
        for (int i = 0; i < position.length; i++) {
            het = new OpenBitSet(gt.allelePresenceForAllTaxa(i, WHICH_ALLELE.Major).getBits().clone());
            bsMa = new OpenBitSet(gt.allelePresenceForAllTaxa(i, WHICH_ALLELE.Major).getBits());
            bsMi = new OpenBitSet(gt.allelePresenceForAllTaxa(i, WHICH_ALLELE.Minor).getBits());
            het.and(bsMi);
            bsMa.xor(het);
            bsMi.xor(het);
            obsMajor[i] = bsMa;
            obsMinor[i] = bsMi;
            majorCount = bsMa.cardinality();
            minorCount = bsMi.cardinality();
            maf[i] = (double)minorCount/(majorCount+minorCount);
            if (i%10000 == 0) System.out.println(String.valueOf(i+1)+" sites are converted to bits. " +String.valueOf((double)(i+1)/position.length) + " completed");
        }
        wordNum = obsMajor[0].getNumWords();
        System.out.println("Transform Genotype to sBit took " + this.getTimeSpanSecond(lastTimePoint) +  " seconds");
        this.writeH5File(sBitFileS);
    }
    
    /**
     * Read in HDF5 SimpleGenotypeSBit
     * @param inputfileS 
     */
    public SimpleGenotypeSBit (String inputfileS) {
        this.readH5File(inputfileS);
    }
    
    public void writeH5File (String outputFileS) {
        long lastTimePoint  = this.getCurrentTimeNano();
        IHDF5WriterConfigurator config = HDF5Factory.configure(new File(outputFileS));
        config.overwrite();
        config.useUTF8CharacterEncoding();
        IHDF5Writer h5 = config.writer();
        h5.int32().setAttr(ROOT, TAXANUM, taxaNum);
        h5.int32().setAttr(ROOT, CHRNUM, chrNum);
        h5.int32().setAttr(ROOT, SITENUM, siteNum);
        h5.int32().setAttr(ROOT, WORDNUM, wordNum);
        h5.string().createArray(TAXANAMES, maxTaxaNameLength, taxaNames.length, genericFeature);
        h5.string().writeArray(TAXANAMES, taxaNames, genericFeature);
        h5.int32().createArray(CHROMOSOME, chrNum, intFeature);
        h5.int32().writeArray(CHROMOSOME, chromosomeNumber, intFeature);
        h5.int32().createArray(CHRSTARTINDEX, chrNum, intFeature);
        h5.int32().writeArray(CHRSTARTINDEX, chrStartIndex, intFeature);
        h5.int32().createArray(CHRENDINDEX, chrNum, intFeature);
        h5.int32().writeArray(CHRENDINDEX, chrEndIndex, intFeature);
        h5.int32().createArray(POSITION, siteNum, intFeature);
        h5.int32().writeArray(POSITION, position, intFeature);
        h5.float64().createArray(MAF, siteNum, floatFeature);
        h5.float64().writeArray(MAF, maf, floatFeature);
        long[][] dis = new long[siteNum][];
        for (int i = 0; i < siteNum; i++) {
            dis[i] = obsMajor[i].getBits();
        }
        h5.int64().createMatrix(OBSMAJOR, siteNum, wordNum, this.getChunkSize(), wordNum, intFeature);
        for (int i = 0; i < this.getChunkNum(); i++) {
            long[][] chunk = this.getSubLongMatrix(dis, this.getChunkSize(), i);
            h5.int64().writeMatrixBlock(OBSMAJOR, chunk, i, 0);
        }
        for (int i = 0; i < siteNum; i++) {
            dis[i] = obsMinor[i].getBits();
        }
        h5.int64().createMatrix(OBSMINOR, siteNum, wordNum, this.getChunkSize(), wordNum, intFeature);
        for (int i = 0; i < this.getChunkNum(); i++) {
            long[][] chunk = this.getSubLongMatrix(dis, this.getChunkSize(), i);
            h5.int64().writeMatrixBlock(OBSMINOR, chunk, i, 0);
        }
        h5.file().flush();
        h5.close();
        System.out.println("Write to SimpleGenotypeSBit took " + this.getTimeSpanSecond(lastTimePoint) +  " seconds");
    }
    
    private long[][] getSubLongMatrix (long[][] dis, int actualChunkSize, int blockNumberX) {
        long[][] result = new long[actualChunkSize][];
        int chunkStartSiteIndex = actualChunkSize*blockNumberX;
        if (chunkStartSiteIndex+actualChunkSize > this.siteNum) {
            actualChunkSize = -(chunkStartSiteIndex-this.siteNum);
            for (int i = actualChunkSize; i < this.getChunkSize(); i++) {
                result[i] = new long[dis[0].length];
                for (int j = 0; j < dis[0].length; j++) {
                    result[i][j] = Long.MIN_VALUE;
                }
            }
        }
        for (int i = 0; i < actualChunkSize; i++) {
            result[i] = dis[chunkStartSiteIndex+i];
        }
        return result;
    }
    
    
    private long[][] readInSBitChunk (IHDF5Writer h5, String path) {
        long[][] dis = new long[this.siteNum][this.wordNum];
        long[][] sub = new long[this.getChunkSize()][this.wordNum];
        int actualSize = this.getChunkSize();
        int chunkStartSiteIndex;
        for (int i = 0; i < this.getChunkNum(); i++) {
            chunkStartSiteIndex = this.getChunkSize()*i;
            actualSize = this.getChunkSize();
            sub = h5.int64().readMatrixBlock(path, this.getChunkSize(), this.wordNum, i, 0);
            if (chunkStartSiteIndex + this.getChunkSize() > this.siteNum) {
                actualSize = siteNum-chunkStartSiteIndex;
            }
            for (int j = 0; j < actualSize; j++) {
                dis[chunkStartSiteIndex+j] = sub[j];
            }
        }
        return dis;
    }
    
    private class MTReadInSBitChunk implements Runnable {
        IHDF5Writer h5;
        String path;
        long[][] dis;
        int chunkIndex;
        public MTReadInSBitChunk (IHDF5Writer h5, String path, long[][] dis, int chunkIndex) {
            this.h5 = h5;
            this.path = path;
            this.dis = dis;
            this.chunkIndex = chunkIndex;
        }
        
        @Override
        public void run() {
            int chunkStartSiteIndex = getChunkSize()*chunkIndex;
            int actualSize = getChunkSize();
            long[][] sub = h5.int64().readMatrixBlock(path, getChunkSize(), wordNum, chunkIndex, 0);
            if (chunkStartSiteIndex + getChunkSize() > siteNum) {
                actualSize = siteNum-chunkStartSiteIndex;
            }
            for (int i = 0; i < actualSize; i++) {
                dis[chunkStartSiteIndex+i] = sub[i];
            }
            h5.close();
        }
        
    }
    
    public void readH5File (String inputFileS) {
        long lastTimePoint = this.getCurrentTimeNano();
        IHDF5Writer h5 = HDF5Factory.open(inputFileS);
        taxaNum = h5.int32().getAttr(ROOT, TAXANUM);
        chrNum = h5.int32().getAttr(ROOT, CHRNUM);
        siteNum = h5.int32().getAttr(ROOT, SITENUM);
        wordNum = h5.int32().getAttr(ROOT, WORDNUM);
        taxaNames = h5.readStringArray(TAXANAMES);
        chromosomeNumber = h5.readIntArray(CHROMOSOME);
        chrStartIndex = h5.readIntArray(CHRSTARTINDEX);
        chrEndIndex = h5.readIntArray(CHRENDINDEX);
        position = h5.readIntArray(POSITION);
        maf = h5.readDoubleArray(MAF);
        obsMajor = new OpenBitSet[siteNum];
        obsMinor = new OpenBitSet[siteNum];

        
//MT readin##############################################/        
        h5.close();
        ExecutorService pool = Executors.newFixedThreadPool(Runtime.getRuntime().availableProcessors());
        long[][] dis = new long[siteNum][];
        for (int i = 0; i < this.getChunkNum(); i++) {
            pool.execute(new SimpleGenotypeSBit.MTReadInSBitChunk(HDF5Factory.open(inputFileS), OBSMAJOR, dis, i));
            //new MTReadInSBitChunk(h5, OBSMAJOR, dis, i).run();
        }
        pool.shutdown();
        try {
            pool.awaitTermination(Integer.MAX_VALUE, TimeUnit.SECONDS);
        } 
        catch (InterruptedException e) {
           System.out.println(e.toString());
        }
        for (int i = 0; i < siteNum; i++) {
            obsMajor[i] = new OpenBitSet(dis[i]);
        }
        
        pool = Executors.newFixedThreadPool(Runtime.getRuntime().availableProcessors());
        dis = new long[siteNum][];
        for (int i = 0; i < this.getChunkNum(); i++) {
            pool.execute(new SimpleGenotypeSBit.MTReadInSBitChunk(HDF5Factory.open(inputFileS), OBSMINOR, dis, i));
            //new MTReadInSBitChunk(h5, OBSMINOR, dis, i).run();
        }
        pool.shutdown();
        try {
            pool.awaitTermination(Integer.MAX_VALUE, TimeUnit.SECONDS);
        } 
        catch (InterruptedException e) {
           System.out.println(e.toString());
        }
        for (int i = 0; i < siteNum; i++) {
            obsMinor[i] = new OpenBitSet(dis[i]);
        }
//#################################################################/
/*        
//Readin###########################################################/        
        long[][] dis = this.readInSBitChunk(h5, OBSMAJOR);
        for (int i = 0; i < siteNum; i++) {
            obsMajor[i] = new OpenBitSet(dis[i]);
        }
        dis = this.readInSBitChunk(h5, OBSMINOR);
        for (int i = 0; i < siteNum; i++) {
            obsMinor[i] = new OpenBitSet(dis[i]);
        }
        h5.close();
//#################################################################/
*/ 
        System.out.println("Read SimpleGenotypeSBit took " + this.getTimeSpanSecond(lastTimePoint) +  " seconds");
    }
    
    public int getSiteNum () {
        return this.siteNum;
    }
    
    public int getTaxaNum () {
        return this.taxaNum;
    }
    
    public int getTaxonIndex (String taxonName) {
        for (int i = 0; i < taxaNum; i++) {
            if (taxaNames[i].equals(taxonName)) return i;
        }
        return -1;
    }
    
    private double getTimeSpanSecond (long lastTimePoint) {
        return (double)this.getTimeSpanNano(lastTimePoint)/1000000000;
    }
    
    private long getTimeSpanNano (long lastTimePoint) {
        return this.getCurrentTimeNano()- lastTimePoint;
    }
    
    private long getCurrentTimeNano () {
        return System.nanoTime();
    }
    
    /**
     * Return index of nearest site of given position
     * @param chr
     * @param pos
     * @return 
     */
    public int getSiteIndex (int chr, int pos) {
        int chrIndex = Arrays.binarySearch(this.chromosomeNumber, chr);
        if (chrIndex < 0) return Integer.MIN_VALUE;
        int siteIndex = Arrays.binarySearch(position, chrStartIndex[chrIndex], chrEndIndex[chrIndex], pos);
        if (siteIndex < 0) {
            siteIndex = - siteIndex - 2;
        }
        if (siteIndex < chrStartIndex[chrIndex]) siteIndex = chrStartIndex[chrIndex];
        if (siteIndex > chrEndIndex[chrIndex]) siteIndex = chrEndIndex[chrIndex];
        return siteIndex;
    }
    
    public int[] getAdjacentSiteIndexRange (int siteIndex, int siteNum) {
        int currentChr = this.getChr(siteIndex);
        int half = siteNum/2;
        int[] indexRange = new int[2];
        indexRange[0] = siteIndex-half;
        indexRange[1] = siteIndex+half+1;
        if (indexRange[0] < 0 || this.getChr(indexRange[0]) != currentChr) {
            indexRange[0] = this.chrStartIndex[this.getChrIndex(currentChr)];
        }
        if (indexRange[1] >= this.getSiteNum() || this.getChr(indexRange[1]) != currentChr) {
            indexRange[1] = this.chrEndIndex[this.getChrIndex(currentChr)];
        }
        return indexRange;
    }
    
    public int getChrIndex (int chr) {
        return Arrays.binarySearch(this.chromosomeNumber, chr);
    }
    
    public int getChr (int index) {
        int chrIndex = Arrays.binarySearch(this.chrStartIndex, index);
        if (chrIndex < 0) chrIndex = -chrIndex - 2;
        return this.chromosomeNumber[chrIndex];
    }
    
   /**
     * Return the total number of chunks
     * @return 
     */
    public int getChunkNum () {
        int num = siteNum/this.getChunkSize();
        if (siteNum% this.getChunkSize() == 0) return num;
        else return num+1;
    }
    
    /**
     * Return the chunk size (Number of tags in gt chunk)
     * @return 
     */
    public int getChunkSize () {
        return this.CHUNK_SIZE;
    }
    
    public int getPosition (int siteIndex) {
        return position[siteIndex];
    }
    
    public int getChromosomeNumber (int siteIndex) {
        int hit = Arrays.binarySearch(this.chrStartIndex, siteIndex);
        if (hit < 0) hit = -hit-2;
        return this.chromosomeNumber[hit];
    }
    
    public void writeBinaryFile (String outputFileS) {
        try {
            DataOutputStream dos = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outputFileS), 65536));
            dos.writeInt(taxaNames.length);
            dos.writeInt(chrNum);
            dos.writeInt(siteNum);
            dos.write(wordNum);
            for (int i = 0; i < taxaNames.length; i++) {
                dos.writeUTF(taxaNames[i]);
            }
            for (int i = 0; i < this.chromosomeNumber.length; i++) {
                dos.writeInt(this.chromosomeNumber[i]);
                dos.writeInt(this.chrStartIndex[i]);
                dos.writeInt(this.chrEndIndex[i]);
            }
            long[] bits;
            for (int i = 0; i < maf.length; i++) {
                dos.writeInt(position[i]);
                dos.writeDouble(maf[i]);
                bits = obsMajor[i].getBits();
                for (int j = 0; j < bits.length; j++) {
                    dos.writeLong(bits[j]);
                }
                bits = obsMinor[i].getBits();
                for (int j = 0; j < bits.length; j++) {
                    dos.writeLong(bits[j]);
                }
            }
            dos.flush();
            dos.close();
        }
        catch (Exception e) {
            System.out.println(e.toString());
            System.exit(1);
        }
    }
    
    public void readBinaryFile (String inputFileS) {
        try {
            DataInputStream dis = new DataInputStream(new BufferedInputStream(new FileInputStream(inputFileS), 65536));
            taxaNum = dis.readInt();
            chrNum = dis.readInt();
            siteNum = dis.readInt();
            wordNum = dis.readInt();
            this.initialize();
            for (int i = 0; i < taxaNum; i++) {
                taxaNames[i] = dis.readUTF();
            }
            for (int i = 0; i < chrNum; i++) {
                this.chromosomeNumber[i] = dis.readInt();
                this.chrStartIndex[i] = dis.readInt();
                this.chrEndIndex[i] = dis.readInt();
            }
            long[] bits = new long[wordNum];
            for (int i = 0; i < siteNum; i++) {
                position[i] = dis.readInt();
                maf[i] = dis.readDouble();
                for (int j = 0; j < wordNum; j++) {
                    bits[j] = dis.readLong();
                }
                obsMajor[i] = new OpenBitSet(bits);
                for (int j = 0; j < wordNum; j++) {
                    bits[j] = dis.readLong();
                }
                obsMinor[i] = new OpenBitSet(bits);
            }
        }
        catch (Exception e) {
            System.out.println(e.toString());
            System.exit(1);
        }
    }
    
    private void initialize () {
        taxaNames = new String[taxaNum];
        chromosomeNumber = new int[chrNum];
        chrStartIndex = new int[chrNum];
        chrEndIndex = new int[chrNum];
        position = new int[siteNum];
        maf = new double[siteNum];
        obsMajor = new OpenBitSet[siteNum];
        obsMinor = new OpenBitSet[siteNum];
    }
}

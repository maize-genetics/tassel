package net.maizegenetics.dna.tag;

import ch.systemsx.cisd.hdf5.*;
import net.maizegenetics.constants.GBSConstants;
import net.maizegenetics.constants.GeneralConstants;
import org.junit.Test;
import org.junit.Ignore;

import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Random;

/**
 * Created by edbuckler on 5/19/14.
 */
public class TBTCompressionTest {
    private String testInTBT= GBSConstants.GBS_TEMP_SEQ_TO_TBT_HDF5_PLUGIN_FILE;
    //private String testInTBT= "/Users/edbuckler/SolexaAnal/GBS27/tbt/node11_ZM_c10_TBT_20120623.h5";
    private static final String outName= GeneralConstants.TEMP_DIR+"RANDcompTest.h5";
    private static String outFile;

    private static int step=8;

    static {
        Random r=new Random();
        //outFile=outName.replace("RAND",""+r.nextInt(1000));
        outFile=outName.replace("RAND",""+876);
    }


    @Test
    @Ignore
    public void testStandardSize() throws Exception {
        TagsByTaxaByteHDF5TaxaGroups expectedTBT = new TagsByTaxaByteHDF5TaxaGroups(testInTBT);
        System.out.printf("Input File:  Taxa:%d Tags:%d%n", expectedTBT.getTaxaCount(), expectedTBT.getTagCount());
        Path path= Paths.get(testInTBT);
        long origFileSize=Files.size(path);
        System.out.printf("File size:%d%n", origFileSize);
        IHDF5WriterConfigurator config = HDF5Factory.configure(outFile);
        config.dontUseExtendableDataTypes();
        IHDF5Writer writer = config.writer();
        int blockSize=(1<<16 < expectedTBT.getTagCount())?1<<16:expectedTBT.getTagCount();


        HDF5IntStorageFeatures intDeflation = HDF5IntStorageFeatures.createDeflation(1);
        writer.int64().createMatrix("tags", expectedTBT.getTagSizeInLong(), expectedTBT.getTagCount(),
                expectedTBT.getTagSizeInLong(), expectedTBT.getTagCount());
        writer.writeLongMatrix("tags", expectedTBT.tags);
        writer.int8().createMatrix("base", expectedTBT.getTaxaCount(), expectedTBT.getTagCount(),
                step, blockSize, intDeflation);
        for (int tx = 0; tx < expectedTBT.getTaxaCount(); tx+=step) {
            byte[][] reads=new byte[step][];
            for (int i = 0; i < step; i++) {
                reads[i]=expectedTBT.getReadCountDistributionForTaxon(tx+i);
            }

            //System.out.println(tx+":"+reads[0][0]);
            //System.out.println(expectedTBT.getReadCountForTagTaxon(0,tx));
            writer.int8().writeMatrixBlock("base", reads, tx/step, 0);
        }
        writer.close();
        Path pathNew= Paths.get(outFile);
        long newFileSize=Files.size(pathNew);
        System.out.printf("File size:%d%n", newFileSize);
    }

    @Ignore
    @Test
    public void testReadRate() throws Exception {
        TagsByTaxaByteHDF5TaxaGroups expectedTBT = new TagsByTaxaByteHDF5TaxaGroups(testInTBT);
        System.out.printf("Input File:  Taxa:%d Tags:%d%n", expectedTBT.getTaxaCount(), expectedTBT.getTagCount());
        Path path= Paths.get(testInTBT);
        long origFileSize=Files.size(path);
        System.out.printf("File size:%d%n", origFileSize);
        long times=System.nanoTime();
        long sum=0;
        for (int tx = 0; tx < expectedTBT.getTaxaCount(); tx++) {
            byte[] reads=expectedTBT.getReadCountDistributionForTaxon(tx);
            for (byte read : reads) {
                sum+=read;
            }
        }
        System.out.printf("Sum of all reads:%d Time:%d%n",sum, (System.nanoTime()-times));
        times=System.nanoTime();
        System.out.println("Reading 85compTest.h5");
        IHDF5Reader reader=HDF5Factory.openForReading(outFile);
        sum=0;
        for (int tx = 0; tx < expectedTBT.getTaxaCount(); tx+=step) {
            //System.out.println("tx:"+tx);
            byte[][] reads = reader.int8().readMatrixBlock("base",step,expectedTBT.getTagCount(),tx/step,0l);
            for (byte[] readt : reads) {
                for (byte read : readt) {
                    sum+=read;
                }
            }

        }
        System.out.printf("Sum of all reads:%d Time:%d%n", sum, (System.nanoTime() - times));
//        Random r=new Random();
//        String outFile=outName.replace("RAND",""+r.nextInt(1000));
//        IHDF5WriterConfigurator config = HDF5Factory.configure(outFile);
//        config.dontUseExtendableDataTypes();
//        IHDF5Writer writer = config.writer();
//        int blockSize=(1<<16 < expectedTBT.getTagCount())?1<<16:expectedTBT.getTagCount();
//        writer.createByteMatrix("base", expectedTBT.getTaxaCount(), expectedTBT.getTagCount(),
//                1, blockSize, Tassel5HDF5Constants.intDeflation);
//        for (int tx = 0; tx < expectedTBT.getTaxaCount(); tx++) {
//            byte[][] reads=new byte[1][];
//            reads[0]=expectedTBT.getReadCountDistributionForTaxon(tx);
//            System.out.println(tx+":"+reads[0][0]);
//            //System.out.println(expectedTBT.getReadCountForTagTaxon(0,tx));
//            writer.writeByteMatrixBlockWithOffset("base", reads, tx, 0);
//        }
//        writer.close();
//        Path pathNew= Paths.get(outFile);
//        long newFileSize=Files.size(pathNew);
//        System.out.printf("File size:%d%n", newFileSize);
    }
}

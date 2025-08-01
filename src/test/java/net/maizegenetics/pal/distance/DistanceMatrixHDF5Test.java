package net.maizegenetics.pal.distance;

import ch.systemsx.cisd.hdf5.*;
import org.apache.commons.lang3.ArrayUtils;
import org.junit.Test;
import ch.systemsx.cisd.base.mdarray.MDIntArray;
import java.io.File;
import java.util.Random;
import static org.junit.Assert.assertArrayEquals;

/**
 * User: dek29
 * Date: 8/22/13
 * Time: 1:07 PM
 */
public class DistanceMatrixHDF5Test {


    /**
     * Contains simplifications over previous tests (square only) and looping
     */
    @Test
    public void fifthWriteTest(){
        System.out.println("Running DistanceMatrixHDF5Test: Writing FIFTH matrix block test");
        int testNumber = 5;

        int[] deflationLevel = { 0,1,2,3,5,7,9 };
        int scaling = 3;
//        String scaling = "none";
        // number of objects in each dimension of a block of data
        // arrays are aligned so as you give 10240 x 10240 matrix
        int[] dataDim =  { 32, 64, 128, 256, 512, 1024, 2048, 10240 };
        // number of blocks in each dimension
        int[] blockCount = { 320, 160, 80, 40, 20, 10, 5, 1 };

        for(int h = 0; h < deflationLevel.length; h++) {
            for (int i = 0; i < dataDim.length; i++) {
                // number of object in each dimension across full file
                int fileDim = dataDim[i] * blockCount[i];

                int totalObjects = fileDim * fileDim;

                int totalWrites = (int) Math.ceil((double) totalObjects / (double) (dataDim[i] * dataDim[i] * 2)); //2 for triangular matrix


                float[][] mydata = new float[dataDim[i]][dataDim[i]];

                String dir = "/local/hdf5/test/";

                String filename = dataDim[i] + "_" + dataDim[i] + "matrix" + "writeBlocks" + blockCount[i] + "_" + blockCount[i] + "def_Scl" +deflationLevel[h] + "_" + scaling +"test" + testNumber + ".h5";
                File hdf5File = new File(dir + filename);
                IHDF5WriterConfigurator config = HDF5Factory.configure(hdf5File);
                config.overwrite();
                config.dontUseExtendableDataTypes();
                IHDF5Writer h5w = config.writer();
                HDF5FloatStorageFeatures features = HDF5FloatStorageFeatures.createDeflateAndFloatScaling(deflationLevel[h], scaling);
//            HDF5FloatStorageFeatures features=HDF5FloatStorageFeatures.createDeflation(deflationLevel[h]);

                // Define the block size
                h5w.float32().createMatrix("mydata", fileDim, fileDim, dataDim[i], dataDim[i], features);


                long time = System.nanoTime();
                for (int bx = 0; bx < blockCount[i]; bx++) {
                    for (int by = 0; by < blockCount[i]; by++) {
                        if (bx > by) continue;
                        fillMatrix( mydata);
                        int xOffset = bx * dataDim[i];
                        int yOffset = by * dataDim[i];
                        //                System.out.println("Offsets x: " + xOffset + " y: " + yOffset);
                        h5w.float32().writeMatrixBlockWithOffset("mydata", mydata, dataDim[i], dataDim[i], xOffset, yOffset);
                    }
                }


                long totalTime = System.nanoTime() - time;
                double timePerObj = (double) totalTime / (double) (totalObjects);

                time = System.nanoTime();
                for (int bx = 0; bx < blockCount[i]; bx++) {
                    for (int by = 0; by < blockCount[i]; by++) {
                        if (bx > by) continue;
                        int xOffset = bx * dataDim[i];
                        int yOffset = by * dataDim[i];
                        //                System.out.println("Offsets x: " + xOffset + " y: " + yOffset);
                        float[][] in = h5w.float32().readMatrixBlockWithOffset("mydata", dataDim[i], dataDim[i], xOffset, yOffset);
                        //            assertArrayEquals(mydata[i],in[0],0.001f);
                    }
                }
                h5w.close();
                totalTime = System.nanoTime() - time;
                System.out.print("TotalObjects: " + fileDim + "x" + fileDim + "\tBlockSize: " + dataDim[i] + "x" + dataDim[i]);
                System.out.print("\tTotalWrites: " + totalWrites);
                long fileSizeBytes = hdf5File.length();
                long fileSizeMB = fileSizeBytes / (1024 * 1024);
                System.out.print("\tFileName: " + filename + "\tFileSize(MB):\t" + fileSizeMB + "\tDeflationLevel:\t" + deflationLevel[h] +
                        "\tScaling:\t" + scaling +"\tHDF5WritingTime(s):\t" + totalTime / 1e9 + "\tAvgPerObj(ns):\t" + timePerObj);
                timePerObj = (double) totalTime / (double) (totalObjects);
                System.out.println("\tHDF5ReadingTime(s):\t" + totalTime / 1e9 + "\tAvgPerObj(ns):\t" + timePerObj);

            }
        }
    }



    public void secondWriteTest(){
        System.out.println("Running DistanceMatrixHDF5Test: Writing SECOND matrix block test");
        int testNumber = 2;
        int xFileBlocks = 16;
        int yFileBlocks = 16;
        int dimension = 8;
        int x_block = 128;
        int y_block = 128;

        int xDataBlockSize = dimension*x_block;
        int yDataBlockSize = dimension*y_block;
        int xFileSize = xDataBlockSize*xFileBlocks;
        int yFileSize = yDataBlockSize*yFileBlocks;

        Random rng = new Random();
//        float[][] mydata = new float[dimension*x_block][dimension*y_block];
        float[][] mydata = new float[dimension*x_block][dimension*y_block];
        fillMatrix( mydata);

        String filename = xDataBlockSize+ "_" + yDataBlockSize + "matrix" + dimension + "writeBlocks" + xFileBlocks + "_" + yFileBlocks + "test" + testNumber + ".h5";
        File hdf5File = new File(filename);
        IHDF5WriterConfigurator config = HDF5Factory.configure(hdf5File);
        config.overwrite();
        config.dontUseExtendableDataTypes();
        IHDF5Writer h5w = config.writer();
        //  HDF5FloatStorageFeatures features=HDF5FloatStorageFeatures.createDeflateAndFloatScaling(9,3);
        HDF5FloatStorageFeatures features=HDF5FloatStorageFeatures.createDeflation(0);
        double totalWrites = (double)(xFileSize * yFileSize)/ (double)(xDataBlockSize * yDataBlockSize);
        // Define the block size
        h5w.float32().createMatrix("mydata", xFileSize, yFileSize, xDataBlockSize, yDataBlockSize, features);
        System.out.println("Created Float Matrix - overall size: " + xFileSize + "x" + yFileSize + " & block sizes: " + xDataBlockSize + "x" + yDataBlockSize);
        System.out.println("Total writes: " + totalWrites);

        long time=System.nanoTime();
        for (int bx = 0; bx < xFileBlocks; bx++){
            for (int by = 0; by < yFileBlocks; by++){
                if(bx > by) continue;
                int xOffset = bx* xDataBlockSize;
                int yOffset = by* yDataBlockSize;
//                System.out.println("Offsets x: " + xOffset + " y: " + yOffset);
               h5w.float32().writeMatrixBlockWithOffset("mydata",mydata,xDataBlockSize, yDataBlockSize,xOffset, yOffset);
            }
        }

        long totalTime=System.nanoTime()-time;
        double timePerObj=(double)totalTime/(double)(dimension*dimension);
        System.out.println(filename + " HDF5WritingTime:"+totalTime/1e9+"s  AvgPerObj:"+timePerObj+"ns");
        time=System.nanoTime();
        for (int bx = 0; bx < xFileBlocks; bx++){
            for (int by = 0; by < yFileBlocks; by++){
                if(bx > by) continue;
                int xOffset = bx* xDataBlockSize;
                int yOffset = by* yDataBlockSize;
//                System.out.println("Offsets x: " + xOffset + " y: " + yOffset);
                float[][] in = h5w.float32().readMatrixBlockWithOffset("mydata",xDataBlockSize, yDataBlockSize,xOffset, yOffset);
  //            assertArrayEquals(mydata[i],in[0],0.001f);
            }
        }
//        for (int i = 0; i < dimension; i++) {
//            float[][] out=new float[1][];
//            out[0]=mydata[i];
//            float[][] in=h5w.readFloatMatrixBlockWithOffset("mydata", 1, dimension, i, 0);
////            assertArrayEquals(mydata[i],in[0],0.001f);
//            //System.out.println(i);
//        }
        totalTime=System.nanoTime()-time;
        timePerObj=(double)totalTime/(double)(dimension*dimension);
        System.out.println(filename + " HDF5ReadingTime:"+totalTime/1e9+"s  AvgPerObj:"+timePerObj+"ns");
        h5w.close();
    }

    public void fillMatrix( float[][] mydata)
    {
        Random rng = new Random();

        for (int i = 0; i < mydata.length; ++i)
        {
            for (int j = 0; j < mydata[i].length; ++j)
            {
                /*if(i>j) {mydata[i][j]=0;}
                else { */mydata[i][j] = (float)rng.nextInt(10000)/10000f;//}
            }
        }
    }


    //    @Ignore
//    @Test
//    public void firstWriteTest(){
//        System.out.println("Running DistanceMatrixHDF5Test: Writing FIRST matrix block test");
//        int testNumber = 1;
//        Random rng = new Random();
//        float[][] mydata = new float[dimension][dimension];
//        fillMatrix(rng, mydata);
//
//        // Write the float matrix.
//        String filename = dimension + "_" + dimension + "matrix" + x_block  + "-"+ y_block + "test" + testNumber + ".h5";
//        IHDF5Writer writer = HDF5Factory.open(filename);
//        // Define the block size
//        writer.createFloatMatrix("mydata", dimension, dimension);
//        // Write 5 x 7 blocks.
//        long time=System.nanoTime();
//        for (int bx = 0; bx < x_block; ++bx)
//        {
//            for (int by = 0; by < y_block; ++by)
//            {
//                writer.writeFloatMatrixBlock("mydata", mydata, bx, by);
//            }
//        }
//        writer.close();
//        long totalTime=System.nanoTime()-time;
//        double timePerObj=(double)totalTime/(double)(dimension*dimension);
//        System.out.println(filename + " HDF5WritingTime:"+totalTime/1e9+"s  AvgPerObj:"+timePerObj+"ns");
//    }

    /**
     *    From http://svncisd.ethz.ch/repos/cisd/jhdf5/tags/release/12.02.x/12.02.0/jhdf5/sourceExamples/java/ch/systemsx/cisd/hdf5/examples/BlockwiseMatrixExample.java
     */
    @Test
    public void startingCodeTest()
    {
        Random rng = new Random();
        int[][] mydata = new int[10][10];

        // Write the integer matrix.
        IHDF5Writer writer = HDF5Factory.open("largeimatrix.h5");
        // Define the block size as 10 x 10.
        writer.int32().createMatrix("mydata", 10, 10);
        // Write 5 x 7 blocks.
        for (int bx = 0; bx < 5; ++bx)
        {
            for (int by = 0; by < 7; ++by)
            {
                fillMatrix(rng, mydata);
                writer.int32().writeMatrixBlock("mydata", mydata, bx, by);
            }
        }
        writer.close();

        // Read the matrix in again, using the "natural" 10 x 10 blocks.
        IHDF5Reader reader = HDF5Factory.openForReading("largeimatrix.h5");
        for (HDF5MDDataBlock<MDIntArray> block : reader.int32().getMDArrayNaturalBlocks("mydata"))
        {
            System.out.println(ArrayUtils.toString(block.getIndex()) + " -> "
                    + block.getData().toString());
        }
        reader.close();
    }

    static void fillMatrix(Random rng, int[][] mydata)
    {
        for (int i = 0; i < mydata.length; ++i)
        {
            for (int j = 0; j < mydata[i].length; ++j)
            {
                mydata[i][j] = rng.nextInt();
            }
        }
    }
}
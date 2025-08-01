/*
 *  AlleleDepthUtilTest
 * 
 *  Created on Feb 24, 2015
 */
package net.maizegenetics.dna.snp.score;

import org.junit.Assert;
import org.junit.Test;

import java.util.Random;

/**
 * @author Ed Buckler
 * @author Terry Casstevens
 */
public class AlleleDepthUtilTest {

    private static final int[] TRANSLATED_DEPTH_VALUES = new int[]{0, 1, 2, 3,
        4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22,
        23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40,
        41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58,
        59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76,
        77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94,
        95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109,
        110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123,
        124, 125, 126, 127, -1, 9763, 9094, 8471, 7892, 7353, 6851, 6384, 5950,
        5545, 5169, 4819, 4493, 4190, 3908, 3645, 3401, 3174, 2962, 2765, 2582,
        2411, 2253, 2105, 1968, 1840, 1721, 1610, 1507, 1411, 1322, 1239, 1161,
        1089, 1023, 960, 902, 848, 798, 752, 708, 668, 630, 595, 562, 532, 504,
        478, 453, 430, 409, 390, 371, 354, 338, 323, 310, 297, 285, 274, 264,
        254, 245, 237, 229, 222, 215, 209, 203, 198, 193, 188, 184, 182, 181,
        180, 179, 178, 177, 176, 175, 174, 173, 172, 171, 170, 169, 168, 167,
        166, 165, 164, 163, 162, 161, 160, 159, 158, 157, 156, 155, 154, 153,
        152, 151, 150, 149, 148, 147, 146, 145, 144, 143, 142, 141, 140, 139,
        138, 137, 136, 135, 134, 133, 132, 131, 130, 129, 128};

    @Test
    /**
     * Test roundtrip of byte-int conversions. Notice the error double by added
     * together.
     */
    public void testByteIntDepthConversion() throws Exception {
        Random r = new Random(1);
        for (int i = 0; i < 9760; i++) {
            byte i2b = AlleleDepthUtil.depthIntToByte(i);
            int i2b2i = AlleleDepthUtil.depthByteToInt(i2b);
            //System.out.printf("%d %d %d%n", i, i2b, i2b2i);
            Assert.assertEquals("Not close for depth:" + i, i, i2b2i, 0.04 * i);

            int rDepth = r.nextInt(9760 - i);
            byte ri2b = AlleleDepthUtil.depthIntToByte(rDepth);
            byte addIandR = AlleleDepthUtil.addByteDepths(i2b, ri2b);
            int addIandRtoI = AlleleDepthUtil.depthByteToInt(addIandR);
            //System.out.printf("%d %d %d%n", i, i+rDepth, addIandRtoI);
            Assert.assertEquals("Not close for added depth:" + i, i + rDepth, addIandRtoI, 0.08 * (i + rDepth));
        }

    }

    @Test
    public void testByteToIntDepth() {
        for (int i = 0; i < 256; i++) {
            Assert.assertEquals("Depth Byte To Int Fails: ", TRANSLATED_DEPTH_VALUES[i], AlleleDepthUtil.depthByteToInt((byte) i));
        }
    }

    @Test
    public void testIntToByteDepth() {
        for (int i = 0; i < 256; i++) {
            Assert.assertEquals("Depth Int To Byte Fails: ", (byte) i, AlleleDepthUtil.depthIntToByte(TRANSLATED_DEPTH_VALUES[i]));
        }
    }

    @Test
    public void testDepthOrder() {
        int current = 0;
        for (int i = 0; i < 30000; i++) {
            int next = AlleleDepthUtil.depthByteToInt(AlleleDepthUtil.depthIntToByte(i));
            Assert.assertTrue("Depths out of order.", next >= current);
            current = next;
        }
    }

}

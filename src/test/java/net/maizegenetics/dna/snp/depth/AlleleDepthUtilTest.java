package net.maizegenetics.dna.snp.depth;

import java.util.Random;
import net.maizegenetics.dna.snp.score.AlleleDepthUtil;
import org.junit.Assert;
import org.junit.Test;

/**
 * Test roundtrip of byte-int conversions. Notice the error double by added
 * together.
 *
 * @author Ed Buckler
 */
public class AlleleDepthUtilTest {

    @Test
    public void testByteIntDepthConversion() throws Exception {
        Random r = new Random(1);
        for (int i = 0; i < 9763; i++) {
            byte i2b = AlleleDepthUtil.depthIntToByte(i);
            int i2b2i = AlleleDepthUtil.depthByteToInt(i2b);
            //System.out.printf("%d %d %d%n", i, i2b, i2b2i);
            Assert.assertEquals("Not close for depth:" + i, i, i2b2i, 0.04 * i);

            int rDepth = r.nextInt(9763 - i);
            byte ri2b = AlleleDepthUtil.depthIntToByte(rDepth);
            byte addIandR = AlleleDepthUtil.addByteDepths(i2b, ri2b);
            int addIandRtoI = AlleleDepthUtil.depthByteToInt(addIandR);
            //System.out.printf("%d %d %d%n", i, i+rDepth, addIandRtoI);
            Assert.assertEquals("Not close for added depth:" + i, i + rDepth, addIandRtoI, 0.08 * (i + rDepth));
        }

    }

}

/**
 *
 */
package net.maizegenetics.dna.map;

import org.junit.Test;

import static org.junit.Assert.assertTrue;

/**
 * @author lcj34
 */
public class ChromosomeTest {
    @Test
    public void ChromosomeCompareToTest() throws Exception {

        System.out.println("Running ChromosomeCompareToTest\n");

        Chromosome one = new Chromosome("1");
        Chromosome nine = new Chromosome("9");
        Chromosome nineA = new Chromosome("9A");
        Chromosome oneA = new Chromosome("1A");
        Chromosome oneD = new Chromosome("1D");
        Chromosome ten = new Chromosome("10");
        Chromosome tenA = new Chromosome("10A");
        Chromosome tenD = new Chromosome("10D");
        Chromosome sevenD = new Chromosome("7D");
        Chromosome a = new Chromosome("A");
        Chromosome chrA = new Chromosome("chrA");
        Chromosome chr1 = new Chromosome("chr1");
        Chromosome chromosome9 = new Chromosome("chromosome9");

        assertTrue(one.compareTo(nine) < 0);
        assertTrue(nine.compareTo(one) > 0);
        assertTrue(oneA.compareTo(nineA) < 0);
        assertTrue(oneA.compareTo(oneD) < 0);
        assertTrue(nine.compareTo(tenA) < 0);
        assertTrue(nine.compareTo(ten) < 0);
        assertTrue(nineA.compareTo(tenA) < 0);
        assertTrue(one.compareTo(a) < 0);
        assertTrue(tenA.compareTo(tenD) < 0);
        assertTrue(sevenD.compareTo(tenD) < 0);

        assertTrue(one.compareTo(chr1) == 0);
        assertTrue(chromosome9.compareTo(a) < 0);
        assertTrue(chrA.compareTo(a) == 0);
        assertTrue(sevenD.compareTo(nine) < 0);
        assertTrue(sevenD.compareTo(nineA) < 0);
        assertTrue(nineA.compareTo(sevenD) > 0);

    }
}

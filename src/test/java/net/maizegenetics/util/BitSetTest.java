/*
 * BitSetTest
 */
package net.maizegenetics.util;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Random;
import java.util.TreeSet;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author Terry Casstevens
 */
public class BitSetTest {

    private final int myNumBits = 5000;
    private int myNumBitsSet = 100;
    private final BitSet myBitSet;
    private final TreeSet<Integer> mySetBits = new TreeSet<Integer>();

    public BitSetTest() {
        myBitSet = new OpenBitSet(myNumBits);
        for (int i = 0; i < myNumBitsSet; i++) {
            int index = (int) Math.round(Math.random() * myNumBits);
            myBitSet.fastSet(index);
            mySetBits.add(index);
        }
        myNumBitsSet = mySetBits.size();
    }

    /**
     * Test of capacity method, of class BitSet.
     */
    public void testCapacity() {
        System.out.println("Testing BitSet Capacity...");
        assertEquals("Capacity Mismatch: ", myNumBits, myBitSet.capacity());
    }

    /**
     * Test of size method, of class BitSet.
     */
    public void testSize() {
        System.out.println("Testing BitSet Size...");
        assertEquals("Size Mismatch: ", myNumBits, myBitSet.size());
    }

    /**
     * Test of isEmpty method, of class BitSet.
     */
    @Test
    public void testIsEmpty() {
        System.out.println("Testing BitSet IsEmpty...");
        assertEquals("IsEmpty Mismatch: ", false, myBitSet.isEmpty());
    }

    /**
     * Test of getNumWords method, of class BitSet.
     */
    @Test
    public void testGetNumWords() {
        System.out.println("Testing BitSet GetNumWords...");
        assertEquals("GetNumWords Mismatch: ", BitUtil.bits2words(myNumBits), myBitSet.getNumWords());
    }

    /**
     * Test of flip method, of class BitSet.
     */
    @Test
    public void testFlip_long_long() {
        System.out.println("Testing BitSet Flip_long_long");
        int numBits = 489;
        BitSet instance = new OpenBitSet(numBits);
        Random random = new Random();
        List<Integer> listSetBits = new ArrayList<>();
        for (int i = 0; i < 100; i++) {
            int current = random.nextInt(100);
            listSetBits.add(current);
            instance.fastSet(current);
        }
        instance.flip(0, numBits);
        for (int i = 0; i < numBits; i++) {
            assertEquals("Flip Test: bit not flipped", instance.fastGet(i), !listSetBits.contains(i));
        }
    }

    /**
     * Test of nextSetBit method, of class BitSet.
     */
    @Test
    public void testNextSetBit_int() {
        System.out.println("Testing BitSet NextSetBit_int...");
        Iterator<Integer> itr = mySetBits.iterator();
        int next = itr.next();
        for (int i = 0; i <= next; i++) {
            assertEquals(next, myBitSet.nextSetBit(i));
        }
        int previous = next;
        while (itr.hasNext()) {
            next = itr.next();
            for (int i = previous + 1; i <= next; i++) {
                assertEquals(next, myBitSet.nextSetBit(i));
            }
            previous = next;
        }
        for (int i = previous + 1; i < myNumBits; i++) {
            assertEquals(-1, myBitSet.nextSetBit(i));
        }
    }

    /**
     * Test of nextSetBit method, of class BitSet.
     */
    @Test
    public void testNextSetBit_long() {
        System.out.println("Testing BitSet NextSetBit_long...");
        Iterator<Integer> itr = mySetBits.iterator();
        long next = itr.next();
        for (long i = 0; i <= next; i++) {
            assertEquals(next, myBitSet.nextSetBit(i));
        }
        long previous = next;
        while (itr.hasNext()) {
            next = itr.next();
            for (long i = previous + 1; i <= next; i++) {
                assertEquals(next, myBitSet.nextSetBit(i));
            }
            previous = next;
        }
        for (long i = previous + 1; i < myNumBits; i++) {
            assertEquals(-1L, myBitSet.nextSetBit(i));
        }
    }

    /**
     * Test of previousSetBit method, of class BitSet.
     */
    @Test
    public void testPreviousSetBit_int() {
        System.out.println("Testing BitSet PreviousSetBit_int...");
        Iterator<Integer> itr = mySetBits.iterator();
        int previous = itr.next();
        for (int i = 0; i < previous; i++) {
            assertEquals(-1, myBitSet.previousSetBit(i));
        }
        while (itr.hasNext()) {
            int next = itr.next();
            for (int i = previous; i < next; i++) {
                assertEquals(previous, myBitSet.previousSetBit(i));
            }
            previous = next;
        }
        for (int i = previous; i < myNumBits; i++) {
            assertEquals(previous, myBitSet.previousSetBit(i));
        }
    }

    /**
     * Test of previousSetBit method, of class BitSet.
     */
    @Test
    public void testPreviousSetBit_long() {
        System.out.println("Testing BitSet PreviousSetBit_long...");
        Iterator<Integer> itr = mySetBits.iterator();
        long previous = itr.next();
        for (long i = 0; i < previous; i++) {
            assertEquals(-1L, myBitSet.previousSetBit(i));
        }
        while (itr.hasNext()) {
            long next = itr.next();
            for (long i = previous; i < next; i++) {
                assertEquals(previous, myBitSet.previousSetBit(i));
            }
            previous = next;
        }
        for (long i = previous; i < myNumBits; i++) {
            assertEquals(previous, myBitSet.previousSetBit(i));
        }
    }

    /**
     * Test of indexOfNthSetBit method, of class BitSet.
     */
    @Test
    public void testIndexOfNthSetBit() {
        System.out.println("Testing BitSet IndexOfNthSetBit...");
        assertEquals("IndexOfNthSetBit: ", -1, myBitSet.indexOfNthSetBit(0));
        Iterator<Integer> itr = mySetBits.iterator();
        int count = 1;
        while (itr.hasNext()) {
            int index = itr.next();
            assertEquals("IndexOfNthSetBit: ", index, myBitSet.indexOfNthSetBit(count));
            count++;
        }
        assertEquals("IndexOfNthSetBit: ", -1, myBitSet.indexOfNthSetBit(count));
    }

    /**
     * Test of getIndicesOfSetBits method, of class BitSet.
     */
    @Test
    public void testGetIndicesOfSetBits() {
        System.out.println("Testing BitSet GetIndicesOfSetBits...");
        int[] indices = myBitSet.getIndicesOfSetBits();
        assertEquals("Number bits set differ: ", myNumBitsSet, indices.length);
        for (int i = 0; i < indices.length; i++) {
            assertEquals("Indices differ: ", true, mySetBits.contains(indices[i]));
        }
    }

}

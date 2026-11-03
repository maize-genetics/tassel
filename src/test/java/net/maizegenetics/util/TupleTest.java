package net.maizegenetics.util;

import org.junit.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotEquals;
import static org.junit.Assert.assertThrows;
import static org.junit.Assert.assertTrue;

public class TupleTest {

    @Test
    public void testCompareToOrdersByXThenY() {
        Tuple<Double, Byte> first = new Tuple<>(0.1, (byte) -3);
        Tuple<Double, Byte> second = new Tuple<>(0.1, (byte) 5);
        Tuple<Double, Byte> third = new Tuple<>(0.2, (byte) -7);

        assertTrue(first.compareTo(second) < 0);
        assertTrue(second.compareTo(first) > 0);
        assertTrue(second.compareTo(third) < 0);
        assertEquals(0, first.compareTo(first));
        assertEquals(0, first.compareTo(new Tuple<>(0.1, (byte) -3)));

        List<Tuple<Double, Byte>> tuples = new ArrayList<>(Arrays.asList(third, second, first));
        Collections.sort(tuples);
        assertEquals(Arrays.asList(first, second, third), tuples);
    }

    @Test
    public void testCompareToUsesYWhenXIsNotComparable() {
        Object x = new Object();
        Tuple<Object, Integer> smaller = new Tuple<>(x, 1);
        Tuple<Object, Integer> larger = new Tuple<>(x, 2);

        assertTrue(smaller.compareTo(larger) < 0);
        assertTrue(larger.compareTo(smaller) > 0);
    }

    @Test
    public void testCompareToRejectsTuplesWithNoComparableValue() {
        Tuple<Object, Object> tuple = new Tuple<>(new Object(), new Object());
        assertThrows(IllegalStateException.class, () -> tuple.compareTo(new Tuple<>(new Object(), new Object())));
    }

    @Test
    public void testEqualsHashCodeAndAccessors() {
        Tuple<String, Integer> tuple = new Tuple<>("chr1", 42);

        assertEquals("chr1", tuple.getX());
        assertEquals(Integer.valueOf(42), tuple.getY());
        assertEquals("(chr1,42)", tuple.toString());
        assertEquals(new Tuple<>("chr1", 42), tuple);
        assertEquals(new Tuple<>("chr1", 42).hashCode(), tuple.hashCode());
        assertNotEquals(new Tuple<>("chr1", 43), tuple);
        assertNotEquals(new Tuple<>("chr2", 42), tuple);
        assertFalse(tuple.equals(null));
        assertFalse(tuple.equals("chr1"));
    }
}

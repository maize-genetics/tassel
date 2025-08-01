/*
 *  ManhattanNumberFormat
 * 
 *  Created on Dec 8, 2016
 */
package net.maizegenetics.analysis.chart;

import java.text.FieldPosition;
import java.text.NumberFormat;
import java.text.ParsePosition;
import java.util.List;

/**
 *
 * @author Terry Casstevens
 */
public class ManhattanNumberFormat extends NumberFormat {

    private final List<Long> myActualPosition;
    private final NumberFormat myBase;

    public ManhattanNumberFormat(NumberFormat base, List<Long> actualPosition) {
        myBase = base;
        myActualPosition = actualPosition;
    }

    @Override
    public StringBuffer format(double number, StringBuffer toAppendTo, FieldPosition pos) {
        return myBase.format(translate((long) number), toAppendTo, pos);
    }

    @Override
    public StringBuffer format(long number, StringBuffer toAppendTo, FieldPosition pos) {
        return myBase.format(translate(number), toAppendTo, pos);
    }

    @Override
    public Number parse(String source, ParsePosition parsePosition) {
        return myBase.parse(source, parsePosition);
    }

    private long translate(long number) {
        long result = 0l;
        for (long current : myActualPosition) {
            if (number >= current) {
                result = current;
            } else {
                break;
            }
        }
        return number - result;
    }

}

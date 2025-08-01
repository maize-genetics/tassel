/*
 *  ReadBedfile
 *
 *  Created on Feb 15, 2017
 */
package net.maizegenetics.dna.snp.io;

import com.google.common.collect.Range;
import com.google.common.collect.RangeMap;
import com.google.common.collect.RangeSet;
import com.google.common.collect.TreeRangeMap;
import com.google.common.collect.TreeRangeSet;
import net.maizegenetics.dna.map.Position;
import net.maizegenetics.dna.map.PositionList;
import net.maizegenetics.dna.map.PositionListBuilder;
import net.maizegenetics.util.Utils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.io.BufferedReader;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

import static java.util.stream.Collectors.collectingAndThen;

/**
 * @author Terry Casstevens
 */
public class ReadBedfile {

    private static final Logger myLogger = LogManager.getLogger(ReadBedfile.class);

    private ReadBedfile() {
        // utility
    }

    /**
     * Function to parse the bedFile and create a List of BedFileRanges.
     *
     * The positions stored in BedFileRange are 1-based inclusive exclusive.
     * This is done by adding 1 to both the start and end position from the BED file.
     * Bedfiles are 0-based.
     *
     * @param bedFile filename of bed file
     * @return
     */
    public static List<BedFileRange> getRanges(String bedFile) {

        List<BedFileRange> result = new ArrayList<>();

        String line = null;
        try (BufferedReader reader = Utils.getBufferedReader(bedFile)) {
            int lineNum = 1;
            line = reader.readLine();
            while (line != null) {

                if (line.startsWith("#")) {
                    line = reader.readLine();
                    lineNum++;
                    continue;
                }

                String[] tokens = line.trim().split("\t");
                if (tokens.length < 3) {
                    throw new IllegalStateException("getRanges: Expecting at least 3 columns on line: " + lineNum);
                }

                // tokens[0] is chromosome
                // tokens[1] is start position from bed file.
                // plus one because bed files are 0-base
                int startPos = Integer.parseInt(tokens[1]) + 1;

                // tokens[2] is start position from bed file.
                // plus one because bed files are 0-base
                int endPos = Integer.parseInt(tokens[2]) + 1;

                // tokens[3] is name from bed file
                String name = null;
                if (tokens.length > 3) {
                    if (tokens[3] == null || tokens[3].isEmpty()) {
                        name = null;
                    } else {
                        name = tokens[3];
                    }
                }

                result.add(new BedFileRange(tokens[0], startPos, endPos, name));

                line = reader.readLine();
                lineNum++;
            }
        } catch (Exception e) {
            myLogger.debug(e.getMessage(), e);
            throw new IllegalStateException("getRanges: problem reading: " + bedFile + " line: " + line);
        }

        return result;

    }

    /**
     * Gets position list from specified bed file.
     *
     * @return position list
     */
    public static PositionList getPositionList(String bedfile) {
        PositionListBuilder builder = new PositionListBuilder();
        getRanges(bedfile).stream().forEach(range -> {
            for (int pos = range.start(); pos < range.end(); pos++) {
                builder.add(Position.of(range.chr(), pos));
            }
        });
        return builder.build();
    }

    /**
     * Function that returns the 1-based Position ranges from a BED file as a RangeSet of Positions.
     * NOTE: getRanges(bedFile) will be called which will shift the start and end positions in the BED file up by 1.
     *       Because of this the ranges returned will be 1-based Closed-Open(Inclusive-Exclusive).
     *       This is NOT returning ranges in BED specification(0-based Inclusive-Exclusive).
     * @param bedfile
     * @return
     */
    public static RangeSet<Position> getRangesAsPositions(String bedfile) {
        return getRanges(bedfile).stream()
                //This needs to be closedOpen as we retain inclusive-exclusive.
                //Also using bedFileRange.myChr instead of bedFileRange.myChrInt as converting a String -> Int -> String is lossy when the string is non-numeric.
                .map(bedFileRange -> Range.closedOpen(Position.of(bedFileRange.myChr, bedFileRange.myStartPos),
                        Position.of(bedFileRange.myChr, bedFileRange.myEndPos)))
                .collect(collectingAndThen(Collectors.toSet(), TreeRangeSet::create));
    }

    /**
     * Function that returns the 1-based closed Position ranges from a BED file as a RangeSet of Positions.
     * NOTE: getRanges(bedFile) will be called which will shift the start in the BED file up by 1.
     *       Because of this the ranges returned will be 1-based Closed(Inclusive-Inclusive).
     *       This is NOT returning ranges in BED specification(0-based Inclusive-Exclusive).
     * @param bedfile
     * @return
     */
    public static RangeSet<Position> getClosedRangesAsPositions(String bedfile) {
        return getRanges(bedfile).stream()
                //This needs to be closedOpen as we retain inclusive-exclusive.
                //Also using bedFileRange.myChr instead of bedFileRange.myChrInt as converting a String -> Int -> String is lossy when the string is non-numeric.
                .map(bedFileRange -> Range.closed(Position.of(bedFileRange.myChr, bedFileRange.myStartPos),
                        Position.of(bedFileRange.myChr, bedFileRange.myEndPos-1)))
                .collect(collectingAndThen(Collectors.toSet(), TreeRangeSet::create));
    }

    /**
     * Function that returns the 1-based Position ranges from a BED file as a RangeMap of Positions to the annotated name of the region.
     * NOTE: getRanges(bedFile) will be called which will shift the start and end positions in the BED file up by 1.
     *       Because of this the ranges returned will be 1-based Closed-Open(Inclusive-Exclusive).
     *       This is NOT returning ranges in BED specification(0-based Inclusive-Exclusive).
     * @param bedfile
     * @return
     */
    public static RangeMap<Position, String> getRangesAsPositionMap(String bedfile) {
        TreeRangeMap<Position, String> positionNameRangeMap = TreeRangeMap.create();
        for (BedFileRange bedFileRange : getRanges(bedfile)) {
            //This needs to be closedOpen as we retain inclusive-exclusive.
            //Also using bedFileRange.myChr instead of bedFileRange.myChrInt as converting a String -> Int -> String is lossy when the string is non-numeric.
            positionNameRangeMap.put(Range.closedOpen(Position.of(bedFileRange.myChr, bedFileRange.myStartPos),
                    Position.of(bedFileRange.myChr, bedFileRange.myEndPos)),
                    bedFileRange.myName);
        }
        return positionNameRangeMap;
    }

    public static class BedFileRange implements Comparable<BedFileRange> {

        private final String myChr;
        private final int myChrInt;
        private final int myStartPos;
        private final int myEndPos;
        private final String myName;

        public BedFileRange(String chr, int startPos, int endPos, String name) {
            myChr = chr;
            int temp;
            try {
                temp = Integer.parseInt(chr);
            } catch (Exception e) {
                temp = -1;
            }
            myChrInt = temp;
            myStartPos = startPos;
            myEndPos = endPos;
            myName = name;
        }

        /**
         * Return chromosome
         *
         * @return chromosome
         */
        public String chr() {
            return myChr;
        }

        /**
         * Returns start position (inclusive)
         *
         * @return start position
         */
        public int start() {
            return myStartPos;
        }

        /**
         * Returns end position (exclusive)
         *
         * @return end position
         */
        public int end() {
            return myEndPos;
        }

        public String name() {
            return myName;
        }

        @Override
        public int compareTo(BedFileRange o) {

            if (myChrInt != -1) {
                if (myChrInt < o.myChrInt) {
                    return -1;
                } else if (myChrInt > o.myChrInt) {
                    return 1;
                }
            } else if (!myChr.equals(o.myChr)) {
                return myChr.compareTo(o.myChr);
            }

            if (myStartPos < o.myStartPos) {
                return -1;
            } else if (myStartPos > o.myStartPos) {
                return 1;
            }

            if (myEndPos < o.myEndPos) {
                return -1;
            } else if (myEndPos > o.myEndPos) {
                return 1;
            } else {
                return 0;
            }

        }

    }

}

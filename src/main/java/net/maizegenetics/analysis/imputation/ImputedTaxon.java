package net.maizegenetics.analysis.imputation;

import com.google.common.collect.ImmutableSortedSet;
import com.google.common.collect.Ordering;
import net.maizegenetics.dna.map.DonorHaplotypes;
import net.maizegenetics.util.BitSet;

import java.util.Arrays;
import java.util.NavigableSet;
import java.util.TreeSet;

/**
 * Created with IntelliJ IDEA.
 * User: edbuckler
 * Date: 2/3/14
 * Time: 2:06 PM
 * To change this template use File | Settings | File Templates.
 */
class ImputedTaxon {
    private final int taxon;
    private final boolean isProjection;
    private boolean segmentSolved=false; //change to true when done
    private final byte[] origGeno;
    byte[] impGeno;  //the best imputed estimate, orginal sequence is not part of this
    byte[] resolveGeno;  //what to set the alignment with combination of original and imp.
    byte[] chgHis; //Viterbi is negative, blockNN is positive
    BitSet[] modBitsOfTarget;
    //byte[][][] allDist;
  //  TreeMap<Integer,int[]> breakPoints;
    private final NavigableSet<DonorHaplotypes> breakPoints;
    private int blocksSolved=0;

    public ImputedTaxon(int taxon, byte[] origGeno, boolean isProjection) {
        this.taxon = taxon;
        this.origGeno = origGeno;
        this.isProjection=isProjection;
        impGeno= Arrays.copyOf(origGeno, origGeno.length); //imputed sequence
        resolveGeno=Arrays.copyOf(origGeno, origGeno.length); //imputed sequence
        chgHis=new byte[origGeno.length];
        breakPoints=new TreeSet<>();
        //breakPoints.put(0, new int[]{-1,-1});
    }

    public void addBreakPoint(DonorHaplotypes dh) {
        breakPoints.add(dh);
    }

    public int taxon() {
        return taxon;
    }

    public NavigableSet<DonorHaplotypes> getBreakPoints() {
        ImmutableSortedSet.Builder<DonorHaplotypes> result=new ImmutableSortedSet.Builder<>(Ordering.natural());
        if(breakPoints.size()==0) return result.build();
        DonorHaplotypes currDH=breakPoints.first();
        for (DonorHaplotypes dhEn : breakPoints) {
            if((currDH.getParent1index()==dhEn.getParent1index())&&(currDH.getParent2index()==dhEn.getParent2index())&&
                (currDH.getChromosome()==dhEn.getChromosome())) {
                currDH=DonorHaplotypes.getMergedInstance(currDH,dhEn);
            } else {
                result.add(currDH);
                currDH=dhEn;
            }
        }
        result.add(currDH);
        return  result.build();
    }

    public byte[] getOrigGeno() {
        return origGeno;
    }

    public byte getOrigGeno(int site) {
        return origGeno[site];
    }

    public byte[] getImpGeno() {
        return impGeno;
    }

    public byte getImpGeno(int site) {
        return impGeno[site];
    }

    public void setImpGeno(int site, byte diploidGenotype) {
        impGeno[site] = diploidGenotype;
    }

    public boolean isProjection() {
        return isProjection;
    }

    public boolean isSegmentSolved() {
        return segmentSolved;
    }

    public void setSegmentSolved(boolean segmentSolved) {
        this.segmentSolved = segmentSolved;
    }

    public int getBlocksSolved() {
        return blocksSolved;
    }

    public void incBlocksSolved() {
        blocksSolved++;
    }



}

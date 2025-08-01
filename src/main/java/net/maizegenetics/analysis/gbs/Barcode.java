
package net.maizegenetics.analysis.gbs;

import java.util.Arrays;

import net.maizegenetics.dna.BaseEncoder;
import net.maizegenetics.taxa.Taxon;

/**
 * Container class for storing information on GBS barcodes.
 * 
 * @author Ed Buckler
 */
public class Barcode implements Comparable<Barcode> {
    /**Barcode sequence */
    String barcodeS;
    /**Overhang sequence from the restriction enzyme */
    String[] overhangS;
    /**Taxon (sample) name */
    String taxaName;
    /**Tissue name */
    String tissueName;
    /**Flowcell name */
    String flowcell;
     /**Flowcell lane name */
    String lane;
    /**Barcode with overhang sequence*/
    String[] barWOverhangS;
    /**Barcode  encoded in 2-bit long*/
    long[] barOverLong;
    /**Length of barcode plus overhang*/
    int barOverLength;
     /**Length of barcode*/
    int barLength;
    /**Global index of taxa based on the key file (first time taxa encountered.*/
    int taxaIndex;
    /** Global index of tissue based on the key file */
    int tissueIndex;

    private Taxon taxon;


    /**
     * Constructor creating a barcode
     * @param barcodeS barcode sequence
     * @param overhangSunsort overhang sequence array unsorted (array size is 1 for
     * non-degerate restriction enzymes, degenerate enzymes >1)
     * @param taxa name of taxon (sample) 
     * @param flowcell name of the flowcell
     * @param lane name of the lane
     */
    public Barcode(String barcodeS, String[] overhangSunsort, String taxa, int globalTaxaIndex, String flowcell, String lane) {
        this.barcodeS = barcodeS;
        Arrays.sort(overhangSunsort);
        this.overhangS = overhangSunsort;
        this.flowcell = flowcell;
        this.lane = lane;
        this.taxaName = taxa;
        this.tissueName = ""; // null or empty string ??
        taxon=new Taxon(taxaName);
        this.taxaIndex=globalTaxaIndex;
        this.tissueIndex = -1;
        barOverLong = new long[overhangS.length];
        barWOverhangS = new String[overhangS.length];
        for (int i = 0; i < overhangS.length; i++) {
            barWOverhangS[i]=barcodeS + overhangS[i];
            barOverLong[i] = BaseEncoder.getLongFromSeq(barWOverhangS[i]);
        }
        barOverLength = barcodeS.length() + overhangS[0].length();
        barLength = barcodeS.length();
    }

    /**
     * Constructor creating a barcode
     * @param barcodeS barcode sequence
     * @param overhangSunsort overhang sequence array unsorted (array size is 1 for
     * non-degerate restriction enzymes, degenerate enzymes >1)
     * @param taxa name of taxon (sample) 
     * @param tissue name identifying the tissue
     * @param flowcell name of the flowcell
     * @param lane name of the lane
     */
    public Barcode(String barcodeS, String[] overhangSunsort, String taxa,  int globalTaxaIndex, 
            String tissue, int globalTissueIndex, String flowcell, String lane) {
        this.barcodeS = barcodeS;
        Arrays.sort(overhangSunsort);
        this.overhangS = overhangSunsort;
        this.flowcell = flowcell;
        this.lane = lane;
        this.taxaName = taxa;
        this.tissueName = null; // null or empty string ??
        taxon=new Taxon(taxaName);
        this.taxaIndex=globalTaxaIndex;
        this.tissueIndex=globalTissueIndex;
        barOverLong = new long[overhangS.length];
        barWOverhangS = new String[overhangS.length];
        for (int i = 0; i < overhangS.length; i++) {
            barWOverhangS[i]=barcodeS + overhangS[i];
            barOverLong[i] = BaseEncoder.getLongFromSeq(barWOverhangS[i]);
        }
        barOverLength = barcodeS.length() + overhangS[0].length();
        barLength = barcodeS.length();
    }

    /**
     * Return the minimum sequence divergence between a query sequence and 
     * barcode with its overhang 
     * @param queryLong query sequence encoded in 2-bit long
     * @param maxDivCheck maximum divergence to search upto
     * @return minimum divergence between barcode and query
     */
    public int compareSequence(long queryLong, int maxDivCheck) {
        int div = barOverLength;
        for (long targetLong : barOverLong) {
            int c = BaseEncoder.seqDifferencesForSubset(targetLong, queryLong, barOverLength, maxDivCheck);
            if (c < div) {
                div = c;
            }
        }
        return div;
    }

    @Override
    public int compareTo(Barcode anotherBarcode) {
        if (this.barOverLong[0] < anotherBarcode.barOverLong[0]) {
            return -1;
        }
        if (this.barOverLong[0] > anotherBarcode.barOverLong[0]) {
            return 1;
        }
        return 0;
    }

    public String getTaxaName() {
        return taxaName;
    }

    public String getBarcodeString() {
        return barcodeS;
    }

    public long[] getBarWOverHangLong() {
        return barOverLong;
    }

    public String[] getBarWOverHang() {
        return barWOverhangS;
    }

    public int getBarWithOverHangLength() {
        return barOverLength;
    }

    public int getBarLength() {
        return barLength;
    }

    public int getTaxaIndex() {
        return taxaIndex;
    }
    
    public int getTissueIndex() {
        return tissueIndex;
    }

    public Taxon getTaxon() {
        return taxon;
    }

    @Override
    public String toString() {
        return "Barcode{" +
                "barcodeS='" + barcodeS + '\'' +
                ", overhangS=" + Arrays.toString(overhangS) +
                ", taxaName='" + taxaName + '\'' +
                ", tissueName='" + tissueName + '\'' +
                ", flowcell='" + flowcell + '\'' +
                ", lane='" + lane + '\'' +
                ", taxaIndex=" + taxaIndex +
                '}';
    }
}

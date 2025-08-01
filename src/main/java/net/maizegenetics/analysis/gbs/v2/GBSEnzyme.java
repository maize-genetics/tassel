package net.maizegenetics.analysis.gbs.v2;

/**
 * Determines which cut sites to look for, and sets them, based on the
 * enzyme used to generate the GBS library. For two-enzyme GBS both enzymes
 * MUST be specified and separated by a dash "-". e.g. PstI-MspI, SbfI-MspI
 * The enzyme pair "PstI-EcoT22I" uses the Elshire common adapter while
 * PstI-MspI, PstI-TaqI, and SbfI-MspI use a Y adapter (Poland et al. 2012)
 *
 * @deprecated use EnzymeList.Enzyme.
 */
public class GBSEnzyme {
    private final String theEnzyme;
    private final String[] initialCutSiteRemnant;
    private final String[] likelyReadEnd;
    private final int readEndCutSiteRemnantLength;
    
    /**
     * @param enzyme The name of the enzyme (case insensitive)
     * @deprecated use EnzymeList.getEnzyme(enzyme)
     */
    public GBSEnzyme(String enzyme) {
        EnzymeList list = new EnzymeList();
        EnzymeList.Enzyme e = list.getEnzyme(enzyme);
        if (e != null) {
            theEnzyme = e.name;
            initialCutSiteRemnant = e.initialCutSiteRemnant;
            likelyReadEnd = e.likelyReadEnd; // full cut site (from partial digest or chimera) or common adapter start
            readEndCutSiteRemnantLength = e.readEndCutSiteRemnantLength;
            return;
        }
        //leaving the original code as a failsafe, but it is redundant with the defaults in the EnzymeList class

        // Check for case-insensitive (?i) match to a known enzyme
        // The common adapter is: [readEndCutSiteRemnant]AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG
        if (enzyme.matches("(?i)apek[i1]")) {
            theEnzyme = "ApeKI";
            initialCutSiteRemnant = new String[]{"CAGC", "CTGC"};
            likelyReadEnd = new String[]{"GCAGC", "GCTGC", "GCAGAGAT", "GCTGAGAT"}; // full cut site (from partial digest or chimera) or common adapter start
            readEndCutSiteRemnantLength = 4;
        } else if (enzyme.matches("(?i)pst[i1]")) {
            theEnzyme = "PstI";
            initialCutSiteRemnant = new String[]{"TGCAG"};
            likelyReadEnd = new String[]{"CTGCAG", "CTGCAAGAT"}; // full cut site (from partial digest or chimera) or common adapter start
            readEndCutSiteRemnantLength = 5;
        } else if (enzyme.matches("(?i)ecot22[i1]")) {
            theEnzyme = "EcoT22I";
            initialCutSiteRemnant = new String[]{"TGCAT"};
            likelyReadEnd = new String[]{"ATGCAT", "ATGCAAGAT"}; // full cut site (from partial digest or chimera) or common adapter start
            readEndCutSiteRemnantLength = 5;
        } else if (enzyme.matches("(?i)pas[i1]")) {
            theEnzyme = "PasI";
            initialCutSiteRemnant = new String[]{"CAGGG", "CTGGG"};
            likelyReadEnd = new String[]{"CCCAGGG", "CCCTGGG", "CCCTGAGAT", "CCCAGAGAT"}; // full cut site (from partial digest or chimera) or common adapter start
            readEndCutSiteRemnantLength = 5;
        } else if (enzyme.matches("(?i)hpaii|(?i)hpa2")) {
            theEnzyme = "HpaII";
            initialCutSiteRemnant = new String[]{"CGG"};
            likelyReadEnd = new String[]{"CCGG", "CCGAGATCGG"}; // full cut site (from partial digest or chimera) or common adapter start
            readEndCutSiteRemnantLength = 3;
        } else if (enzyme.matches("(?i)msp[i1]")) {
            theEnzyme = "MspI";  // MspI and HpaII are isoschizomers (same recognition seq and overhang)
            initialCutSiteRemnant = new String[]{"CGG"};
            likelyReadEnd = new String[]{"CCGG", "CCGAGATCGG"}; // full cut site (from partial digest or chimera) or common adapter start
            readEndCutSiteRemnantLength = 3;
        } else if (enzyme.matches("(?i)pst[i1]-apek[i1]")) {
            theEnzyme = "PstI-ApeKI";
            initialCutSiteRemnant = new String[]{"TGCAG"};
            likelyReadEnd = new String[]{"GCAGC", "GCTGC", "CTGCAG", "GCAGAGAT", "GCTGAGAT"}; // look for ApeKI site, PstI site, or common adapter for ApeKI
            readEndCutSiteRemnantLength = 4;
        } else if (enzyme.matches("(?i)pst[i1]-bfa[i1]")) {
            theEnzyme = "PstI-BfaI";  //  PstI: CTGCA^G  BfaI: C^TAG
            initialCutSiteRemnant = new String[]{"TGCAG"};
            likelyReadEnd = new String[]{"CTGCAG", "CTAG", "CTAAGATC"}; // look for PstI site, BfaI site, or common adapter start for BfaI
            readEndCutSiteRemnantLength = 3;
        } else if (enzyme.matches("(?i)pst[i1]-ecot22[i1]")) {
            theEnzyme = "PstI-EcoT22I";
            initialCutSiteRemnant = new String[]{"TGCAG", "TGCAT"};
            likelyReadEnd = new String[]{"ATGCAT", "CTGCAG", "CTGCAAGAT", "ATGCAAGAT"}; // look for EcoT22I site, PstI site, or common adapter for PstI/EcoT22I
            readEndCutSiteRemnantLength = 5;
        } else if (enzyme.matches("(?i)pst[i1]-msp[i1]")) {
            theEnzyme = "PstI-MspI";
            initialCutSiteRemnant = new String[]{"TGCAG"};
            likelyReadEnd = new String[]{"CCGG", "CTGCAG", "CCGAGATC"}; // look for MspI site, PstI site, or common adapter for MspI
            readEndCutSiteRemnantLength = 3;
        } else if (enzyme.matches("(?i)pst[i1]-msp[i1]-GDFcustom")) {
            theEnzyme = "PstI-MspI-GDFcustom";
            initialCutSiteRemnant = new String[]{"TGCAG"};
            // changed from  CCGAGAT to CCGCTCAGG, as IGD/GDF used a custom Y adapter for MspI
            likelyReadEnd = new String[]{"CCGG", "CTGCAG", "CCGCTCAGG"}; // look for MspI site, PstI site, or GDF custom common adapter for MspI
            readEndCutSiteRemnantLength = 3;
        } else if (enzyme.matches("(?i)pst[i1]-taq[i1]")) {
            theEnzyme = "PstI-TaqI";
            initialCutSiteRemnant = new String[]{"TGCAG"};
            likelyReadEnd = new String[]{"TCGA", "CTGCAG", "TCGAGATC"}; // look for TaqI site, PstI site, or common adapter for TaqI
            readEndCutSiteRemnantLength = 3;
        } else if (enzyme.matches("(?i)nsi[i1]-msp[i1]")) {
            theEnzyme = "NsiI-MspI";  //  ATGCA^T   C^CGG
            initialCutSiteRemnant = new String[]{"TGCAT"};
            likelyReadEnd = new String[]{"CCGG", "ATGCAT", "CCGAGATC"}; // look for MspI site, NsiI site, or common adapter for MspI
            readEndCutSiteRemnantLength = 3;
        } else if(enzyme.matches("(?i)PaeR7[i1]-Hha[i1]")) {
            theEnzyme = "PaeR7I-HhaI";
            // Requested by Ye, Songqing, use same Y adapter as Polland paper  -QS
            initialCutSiteRemnant=new String[]{"TCGAG"};
            likelyReadEnd = new String[]{"GCGC", "CTCGAG", "GCGAGATC"}; // look for HhaI site, PaeR7I site, or common adapter for HhaI
            readEndCutSiteRemnantLength = 3;
        } else if (enzyme.matches("(?i)sbf[i1]-msp[i1]")) {
            theEnzyme = "SbfI-MspI";  // CCTGCA^GG  C^CGG
            initialCutSiteRemnant = new String[]{"TGCAGG"};
            likelyReadEnd = new String[]{"CCGG", "CCTGCAGG", "CCGAGATC"}; // look for MspI site, SbfI site, or common adapter for MspI
            readEndCutSiteRemnantLength = 3;
        } else if (enzyme.matches("(?i)sbf[i1]-hpaii|(?i)sbf[i1]-hpa2")) {
            theEnzyme = "SbfI-HpaII";  // CCTGCA^GG  C^CGG  Nb: HpaII is an isoschizomer of MspI
            initialCutSiteRemnant = new String[]{"TGCAGG"};
            likelyReadEnd = new String[]{"CCGG", "CCTGCAGG", "CCGAGATC"}; // look for HpaII site, SbfI site, or common adapter for HpaII
            readEndCutSiteRemnantLength = 3;
        } else if (enzyme.matches("(?i)sbf[i1]-bfa[i1]")) {
            theEnzyme = "SbfI-BfaI";  // CCTGCA^GG  C^TAG
            initialCutSiteRemnant = new String[]{"TGCAGG"};
            likelyReadEnd = new String[]{"CTAG", "CCTGCAGG", "CTAAGATC"}; // look for BfaI site, SbfI site, or common adapter for BfaI
            readEndCutSiteRemnantLength = 3;
        } else if (enzyme.matches("(?i)sph[i1]-ecor[i1]")) {
            theEnzyme = "SphI-EcoRI";  // GCATG^C  G^AATTC
            initialCutSiteRemnant = new String[]{"CATGC"};  // SphI overhang from the genomic DNA fragment (top strand)
            likelyReadEnd = new String[]{"GCATGC", "GAATTC", "GAATTAGATC"}; // look for SphI site, EcoRI site, or common adapter for EcoRI
            readEndCutSiteRemnantLength = 5;
        } else if (enzyme.matches("(?i)asis[i1]-msp[i1]")) {
            theEnzyme = "AsiSI-MspI";
            initialCutSiteRemnant = new String[]{"ATCGC"};
            likelyReadEnd = new String[]{"CCGG", "GCGATCGC", "CCGAGATC"}; // look for MspI site, AsiSI site, or common adapter for MspI
            readEndCutSiteRemnantLength = 3;
        } else if (enzyme.matches("(?i)bsshii-msp[i1]|(?i)bssh2-msp[i1]")) {
            theEnzyme = "BssHII-MspI";
            initialCutSiteRemnant = new String[]{"CGCGC"};
            likelyReadEnd = new String[]{"CCGG", "GCGCGC", "CCGAGATC"}; // look for MspI site, BssHII site, or common adapter for MspI
            readEndCutSiteRemnantLength = 3;
        } else if (enzyme.matches("(?i)fse[i1]-msp[i1]")) {
            theEnzyme = "FseI-MspI";
            initialCutSiteRemnant = new String[]{"CCGGCC"};
            likelyReadEnd = new String[]{"CCGG", "GGCCGGCC", "CCGAGATC"}; // look for MspI site, FseI site, or common adapter for MspI
            readEndCutSiteRemnantLength = 3;
        } else if (enzyme.matches("(?i)sal[i1]-msp[i1]")) {
            theEnzyme = "SalI-MspI";
            initialCutSiteRemnant = new String[]{"TCGAC"};
            likelyReadEnd = new String[]{"CCGG", "GTCGAC", "CCGAGATC"}; // look for MspI site, SalI site, or common adapter for MspI
            readEndCutSiteRemnantLength = 3;
        } else if (enzyme.matches("(?i)ecor[i1]-msp[i1]")) {
            theEnzyme = "EcoRI-MspI";   //  G^AATTC  C^CGG
            initialCutSiteRemnant = new String[]{"AATTC"};
            likelyReadEnd = new String[]{"CCGG", "GAATTC", "CCGAGATC"}; // look for MspI site, EcoRI site, or Poland et al. 2012 Y-adapter for MspI
            readEndCutSiteRemnantLength = 3;
        } else if (enzyme.matches("(?i)hindiii-msp[i1]|(?i)hind3-msp[i1]")) {
            theEnzyme = "HindIII-MspI"; // A^AGCTT   C^CGG
            initialCutSiteRemnant = new String[]{"AGCTT"};
            likelyReadEnd = new String[]{"CCGG", "AAGCTT", "CCGAGATC"}; // look for MspI site, HindIII site, or Poland et al. 2012 Y-adapter for MspI
            readEndCutSiteRemnantLength = 3;
        } else if (enzyme.matches("(?i)hindiii-nlaiii|(?i)hind3-nla3")) {
            theEnzyme = "HindIII-NlaIII"; // A^AGCTT   ^CATG (not blunt)
            initialCutSiteRemnant = new String[]{"AGCTT"};
            likelyReadEnd = new String[]{"CATG", "AAGCTT", "CATGAGATC"}; // look for NlaIII site, HindIII site, or Poland et al. 2012 Y-adapter for Nla3
            readEndCutSiteRemnantLength = 4;
        } else if (enzyme.matches("(?i)sexa[i1]-sau3a[i1]")) {
            theEnzyme = "SexAI-Sau3AI";  // A^CCWGGT   ^GATC (not blunt)
            initialCutSiteRemnant = new String[]{"CCAGGT", "CCTGGT"};
            likelyReadEnd = new String[]{"GATC", "ACCAGGT", "ACCTGGT", "GATCAGATC"}; // look for SexAI site, Sau3AI site, or Poland et al. 2012 Y-adapter for Sau3AI
            readEndCutSiteRemnantLength = 4;
        } else if (enzyme.matches("(?i)bamh[i1l]-mluc[i1]")) {
            theEnzyme = "BamHI-MluCI";  // G^GATCC   ^AATT (not blunt)
            initialCutSiteRemnant = new String[]{"GATCC"};
            likelyReadEnd = new String[]{"AATT", "GGATCC", "AATTAGATC"}; // look for MluCI site, BamHI site, or Poland et al. 2012 Y-adapter for MluCI
            readEndCutSiteRemnantLength = 4;
        } else if (enzyme.matches("(?i)pst[i1]-mluc[i1]")) {
            theEnzyme = "PstI-MluCI";  // CTGCA^G   ^AATT (not blunt)
            initialCutSiteRemnant = new String[]{"TGCAG"};
            likelyReadEnd = new String[]{"AATT", "CTGCAG", "AATTAGATC"}; // look for MluCI site, PstI site, or Poland et al. 2012 Y-adapter for MluCI
            readEndCutSiteRemnantLength = 4;
        } else if (enzyme.matches("(?i)psti-msei|(?i)pst1-mse1")) {
            theEnzyme = "PstI-MseI"; // CTGCA^G   T^TAA
            initialCutSiteRemnant = new String[]{"TGCAG"};
            likelyReadEnd = new String[]{"TTAA", "CTGCAG", "TTAAGATC"}; // look for MseI site, PstI site, or Poland et al. 2012 Y-adapter for MseI
            readEndCutSiteRemnantLength = 3;
        } else if (enzyme.matches("(?i)avaii-msei|(?i)ava2-mse1")) {
            theEnzyme = "AvaII-MseI"; // G^GWCC   T^TAA  W=AorT
            initialCutSiteRemnant = new String[]{"GACC", "GTCC"};
            likelyReadEnd = new String[]{"TTAA", "GGACC", "GGTCC", "TTAAGATC"}; // look for MseI site, AvaII site, or Poland et al. 2012 Y-adapter for MseI
            readEndCutSiteRemnantLength = 3;
        } else if (enzyme.matches("(?i)ecori-msei|(?i)ecor1-mse1")) {
            theEnzyme = "EcoRI-MseI"; // G^AATTC   T^TAA
            initialCutSiteRemnant = new String[]{"AATTC"};
            likelyReadEnd = new String[]{"TTAA", "GAATTC", "TTAAGATC"}; // look for MseI site, EcoRI site, or Poland et al. 2012 Y-adapter for MseI
            readEndCutSiteRemnantLength = 3;
        } else if (enzyme.matches("(?i)ecori-avaii|(?i)ecor1-ava2")) {
            theEnzyme = "EcoRI-AvaII"; // G^AATTC   G^GWCC
            initialCutSiteRemnant = new String[]{"AATTC"};
            likelyReadEnd = new String[]{"GGACC", "GGTCC", "GAATTC", "GGACAGATC", "GGTCAGATC"}; // look for AvaII site, EcoRI site, or Poland et al. 2012 Y-adapter for AvaII
            readEndCutSiteRemnantLength = 4;
        } else if (enzyme.matches("(?i)ecori-hinfi|(?i)ecor1-hinf1")) {
            theEnzyme = "EcoRI-HinfI"; // G^AATTC   G^ANTC
            initialCutSiteRemnant = new String[]{"AATTC"};
            likelyReadEnd = new String[]{"GAATC", "GACTC", "GAGTC", "GATTC", "GAATTC", "GAATAGATC", "GACTAGATC", "GAGTAGATC", "GATTAGATC"}; // look for HinfI site, EcoRI site, or Poland et al. 2012 Y-adapter for HinfI
            readEndCutSiteRemnantLength = 4;
        } else if (enzyme.matches("(?i)bbvci-mspi|(?i)bbvc1-msp1")) {
            theEnzyme = "BbvCI-MspI"; // CCTCAGC (-5/-2)   C^CGG
            initialCutSiteRemnant = new String[]{"TCAGC"};
            likelyReadEnd = new String[]{"CCGG", "CCTCAGC", "CCGAGATC"}; // look for MspI site, BbvCI site, or Poland et al. 2012 Y-adapter for MspI
            readEndCutSiteRemnantLength = 3;
        } else if (enzyme.matches("(?i)msp[i1]-apek[i1]")) {
            theEnzyme = "MspI-ApeKI";  // C^CGG  G^CWGC
            initialCutSiteRemnant = new String[]{"CGG", "CAGC", "CTGC"};
            likelyReadEnd = new String[]{"CCGG", "GCAGC", "GCTGC", "CCGAGATCGG", "GCAGAGAT", "GCTGAGAT"}; 
            readEndCutSiteRemnantLength = 3; 
        } else if(enzyme.matches("(?i)apo[i1]")){
            theEnzyme = "ApoI";
            initialCutSiteRemnant=new String[]{"AATTC","AATTT"};
            likelyReadEnd = new String[]{"AAATTC","AAATTT","GAATTC","GAATTT","AAATTAGAT","GAATTAGAT"}; // full cut site (from partial digest or chimera) or common adapter start
            readEndCutSiteRemnantLength = 5;
        } else if (enzyme.matches("(?i)BamH[i1l]")) {
            theEnzyme = "BamHI";
            initialCutSiteRemnant = new String[]{"GATCC"};
            likelyReadEnd = new String[]{"GGATCC", "GGATCAGAT"}; // full cut site (from partial digest or chimera) or common adapter start
            // likelyReadEnd = new String[]{"GGATCC", "AGATCGGAA", "AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG"}; // <-- corrected from this by Jeff Glaubitz on 2012/09/12
            readEndCutSiteRemnantLength = 5;
        } else if (enzyme.matches("(?i)mse[i1]")) {
            theEnzyme = "MseI";
            initialCutSiteRemnant = new String[]{"TAA"};
            likelyReadEnd = new String[]{"TTAA", "TTAAGAT"}; // full cut site (from partial digest or chimera) or common adapter start
            readEndCutSiteRemnantLength = 3;
        } else if(enzyme.matches("(?i)Sau3A[i1]")){
            theEnzyme = "Sau3AI";
            initialCutSiteRemnant=new String[]{"GATC"};
            likelyReadEnd = new String[]{"GATC","GATCAGAT"}; // full cut site (from partial digest or chimera) or common adapter start
            readEndCutSiteRemnantLength = 4;
        } else if(enzyme.matches("(?i)nde[i1]")){
            theEnzyme = "NdeI";
            initialCutSiteRemnant=new String[]{"TATG"};
            likelyReadEnd = new String[]{"CATATG","CATAAGAT"}; // full cut site (from partial digest or chimera) or common adapter start
            readEndCutSiteRemnantLength = 4;
        } else if(enzyme.matches("(?i)hinp1[i1]")){
            theEnzyme = "HinP1I";
            initialCutSiteRemnant=new String[]{"CGC"};
            likelyReadEnd = new String[]{"GCGC","GCGAGAT"}; // full cut site (from partial digest or chimera) or common adapter start
            readEndCutSiteRemnantLength = 3;
        } else if(enzyme.matches("(?i)sbf[i1]")){
            theEnzyme = "SbfI";
            initialCutSiteRemnant=new String[]{"TGCAGG"};
            likelyReadEnd = new String[]{"CCTGCAGG","CCTGCAAGAT"}; // full cut site (from partial digest or chimera) or common adapter start
            readEndCutSiteRemnantLength = 6;
        } else if (enzyme.matches("(?i)hindiii|(?i)hind3")) {
            theEnzyme = "HindIII"; // A^AGCTT
            initialCutSiteRemnant=new String[]{"AGCTT"};
            likelyReadEnd = new String[]{"AAGCTT","AAGCTAGAT"}; // full cut site (from partial digest or chimera) or common adapter start
            readEndCutSiteRemnantLength = 5;
        } else if(enzyme.matches("(?i)ecor[i1]")) {
            theEnzyme = "EcoRI";  // G^AATTC
            initialCutSiteRemnant= new String[]{"AATTC"};
            likelyReadEnd = new String[]{"GAATTC","GAATTAGAT"}; // full cut site (from partial digest or chimera) or common adapter start
            readEndCutSiteRemnantLength = 5;
        } else if(enzyme.matches("(?i)cviq[i1]")){  
            theEnzyme = "CviQI";  // CviQI and Csp6I are isoschizomers (same recognition seq and overhang)
            initialCutSiteRemnant=new String[]{"TAC"};
            likelyReadEnd = new String[]{"GTAC","GTAAGATCGG"}; // full cut site (from partial digest or chimera) or common adapter start
            readEndCutSiteRemnantLength = 3;
        } else if(enzyme.matches("(?i)csp6[i1]")){  
            theEnzyme = "Csp6I";  // Csp6I and CviQI are isoschizomers (same recognition seq and overhang)
            initialCutSiteRemnant=new String[]{"TAC"};
            likelyReadEnd = new String[]{"GTAC","GTAAGATCGG"}; // full cut site (from partial digest or chimera) or common adapter start
            readEndCutSiteRemnantLength = 3;
        } else if (enzyme.matches("(?i)nlaiii|(?i)nla3")) {
            theEnzyme = "NlaIII"; // CATG^
            initialCutSiteRemnant=new String[]{"CATG"};
            likelyReadEnd = new String[]{"CATG","CATGAGAT"}; // full cut site (from partial digest or chimera) or common adapter start
            readEndCutSiteRemnantLength = 4;
        } else if(enzyme.matches("(?i)sph[i1]")){
            theEnzyme = "SphI";  // GCATG^C
            initialCutSiteRemnant=new String[]{"CATGC"};
            likelyReadEnd = new String[]{"GCATGC","GCATGAGAT"}; // full cut site (from partial digest or chimera) or common adapter start
            readEndCutSiteRemnantLength = 5;
        } else if(enzyme.matches("(?i)nsp[i1]")){
            theEnzyme = "NspI";  // RCATG^Y
            initialCutSiteRemnant=new String[]{"CATGC","CATGT"};
            likelyReadEnd = new String[]{"ACATGT","GCATGC","ACATGAGAT","GCATGAGAT"}; // full cut site (from partial digest or chimera) or common adapter start
            readEndCutSiteRemnantLength = 5;
        } else if(enzyme.matches("(?i)kpn[i1]")){
            theEnzyme = "KpnI";  // GGTAC^C
            initialCutSiteRemnant=new String[]{"GTACC"};
            likelyReadEnd = new String[]{"GGTACC","GGTACAGAT"}; // full cut site (from partial digest or chimera) or common adapter start
            readEndCutSiteRemnantLength = 5;
        } else if(enzyme.matches("(?i)sty[i1]")){
            theEnzyme = "StyI";  // C^CWWGG
            initialCutSiteRemnant=new String[]{"CAAGG","CATGG","CTAGG","CTTGG"};
            likelyReadEnd = new String[]{"CCAAGG","CCATGG","CCTAGG","CCTTGG","CCAAGAGAT","CCATGAGAT","CCTAGAGAT","CCTTGAGAT"}; // full cut site (from partial digest or chimera) or common adapter start
            readEndCutSiteRemnantLength = 5;
        } else if(enzyme.matches("(?i)styi-msei|(?i)sty1-mse1")){
            theEnzyme = "StyI-MseI";  // C^CWWGG & T^TAA
            initialCutSiteRemnant=new String[]{"CAAGG","CATGG","CTAGG","CTTGG"};
            likelyReadEnd = new String[]{"TTAA","CCAAGG","CCATGG","CCTAGG","CCTTGG","TTAAGAT"}; // full cut site (from partial digest or chimera) or common adapter start
            readEndCutSiteRemnantLength = 3;
        } else if(enzyme.matches("(?i)fse[i1]")){
            theEnzyme = "FseI";  // GGCCGG^CC
            initialCutSiteRemnant=new String[]{"CCGGCC"};
            likelyReadEnd = new String[]{"GGCCGGCC","AGATCGGAAG"}; // full cut site (from partial digest or chimera) or Morishige et al (BMC Genomics, 2013) T adapter start
            readEndCutSiteRemnantLength = 0;  // assumes that common T adapter is far more likely than a second full cut site
        } else if(enzyme.matches("(?i)NgoMIV|(?i)NgoM4")){
            theEnzyme = "NgoMIV";  // G^CCGGC
            initialCutSiteRemnant=new String[]{"CCGGC"};
            likelyReadEnd = new String[]{"GCCGGC","AGATCGGAAG"}; // full cut site (from partial digest or chimera) or Morishige et al (BMC Genomics, 2013) T adapter start
            readEndCutSiteRemnantLength = 0;  // assumes that common T adapter is far more likely than a second full cut site
        } else if(enzyme.matches("(?i)msl[i1]")){
            theEnzyme = "MslI";  // CAYNN^NNRTG  -- has 32 different cut sites (assuming constrained to palindromic YNN -- 32^2 otherwise)
            initialCutSiteRemnant=new String[]{""};
            likelyReadEnd = new String[]{"AGATCGGA"}; // common adapter start only (too many possible cut sites!)
            readEndCutSiteRemnantLength = 0;  
        } else if(enzyme.matches("(?i)ase[i1]")){
            theEnzyme = "AseI";  // AT^TAAT
            initialCutSiteRemnant=new String[]{"TAAT"};
            likelyReadEnd = new String[]{"ATTAAT","ATTAAGAT"};  // full cut site (from partial digest or chimera) or common adapter start
            readEndCutSiteRemnantLength = 4;  
        } else if(enzyme.matches("(?i)avaii|(?i)ava2")){
            theEnzyme = "AvaII";  // G^GWCC
            initialCutSiteRemnant=new String[]{"GACC","GTCC"};
            likelyReadEnd = new String[]{"GGACC","GGTCC","GGACAGAT","GGTCAGAT"};  // full cut site (from partial digest or chimera) or common adapter start
            readEndCutSiteRemnantLength = 4;  
        } else if (enzyme.matches("(?i)kpn[i1]-msp[i1]")) {
            theEnzyme = "KpnI-MspI";  //  KpnI: GGTAC^C  MspI: C^CGG
            initialCutSiteRemnant = new String[]{"GTACC"};
            likelyReadEnd = new String[]{"CCGG", "GGTACC", "CCGAGATC"}; // look for MspI site, KpnI site, or common adapter start for MspI
            readEndCutSiteRemnantLength = 3;
        } else if (enzyme.matches("(?i)RBSTA")) {
            theEnzyme = "RBSTA";
            initialCutSiteRemnant = new String[]{"TA"};
            likelyReadEnd = new String[]{"TTAA", "GTAC", "CTAG", "TTAAGAT", "GTAAGAT", "CTAAGAT"}; // full cut site (from partial digest or chimera) of MseI, CVIQi, XspI or common adapter start
            readEndCutSiteRemnantLength = 3;
        } else if (enzyme.matches("(?i)RBSCG")) {
            theEnzyme = "RBSCG";
            initialCutSiteRemnant = new String[]{"CG"};
            likelyReadEnd = new String[]{"CCGC", "TCGA", "GCGC", "CCGG", "ACGT", "CCGAGAT", "TCGAGAT", "GCGAGAT", "ACGAGAT"}; // full cut site (from partial digest or chimera) of AciI, TaqaI, HinpI, HpaII, HpyCH4IV or common adapter start
            readEndCutSiteRemnantLength = 3;
        } else if(enzyme.matches("(?i)ignore")){
            theEnzyme = "unspecified";  // can be used for new enzymes -- only looks for barcodes and common adapter starts
            initialCutSiteRemnant=new String[]{""};
            likelyReadEnd = new String[]{"AGATCGGA"}; // common adapter start only
            readEndCutSiteRemnantLength = 0;  
        } else {
            enyzmeErrorMessage();
            theEnzyme = null;
            initialCutSiteRemnant = null;
            likelyReadEnd = null; // full cut site (from partial digest or chimera) of AciI, TaqaI, HinpI, HpaII, HpyCH4IV or common adapter start
            readEndCutSiteRemnantLength = -1;
        }
        System.out.println("Enzyme: " + enzyme);
    }

    private void enyzmeErrorMessage() {
            System.out.println("The software didn't recognize your restriction enzyme (-e option).\n"
                    +"Currently, only the following enzymes are recognized for single enzyme digests:\n"
                    +"  ApeKI"    +"\n"
                    +"  ApoI"     +"\n"
                    +"  AseI"     +"\n"
                    +"  AvaII"    +"\n"
                    +"  BamHI"    +"\n"
                    +"  Csp6I"    +"\n"
                    +"  CviQI"    +"\n"
                    +"  EcoRI"    +"\n"
                    +"  EcoT22I"  +"\n"
                    +"  FseI"     +"\n"
                    +"  HindIII"  +"\n"
                    +"  HinP1I"   +"\n"
                    +"  HpaII"    +"\n"
                    +"  KpnI"     +"\n"
                    +"  MseI"     +"\n"
                    +"  MslI"     +"\n"
                    +"  MspI"     +"\n"
                    +"  NdeI"     +"\n"
                    +"  NgoMIV"   +"\n"
                    +"  NlaIII"   +"\n"
                    +"  NspI"     +"\n"
                    +"  PasI"     +"\n"
                    +"  PstI"     +"\n"
                    +"  Sau3AI"   +"\n"
                    +"  SbfI"     +"\n"
                    +"  SphI"     +"\n"
                    +"  StyI"     +"\n"
                    +"  RBSTA"    +"\n"
                    +"  RBSCG"    +"\n"
                    +"  ignore"    +"\n"
                    +"Or the following for two-enzyme digests:\n"
                    +"  AsiSI-MspI"   +"\n"
                    +"  AvaII-MseI"   +"\n"
                    +"  BamHI-MluCI"  +"\n"
                    +"  BbvCI-MspI"   +"\n"
                    +"  BssHII-MspI"  +"\n"
                    +"  EcoRI-AvaII"  +"\n"
                    +"  EcoRI-HinfI"  +"\n"
                    +"  EcoRI-MseI"   +"\n"
                    +"  EcoRI-MspI"   +"\n"
                    +"  FseI-MspI"    +"\n"
                    +"  HindIII-MspI" +"\n"
                    +"  HindIII-NlaIII" +"\n"
                    +"  KpnI-MspI"    +"\n"
                    +"  MspI-ApeKI"   +"\n"
                    +"  NsiI-MspI"    +"\n"
                    +"  PaeR7I-HhaI"  +"\n"
                    +"  PstI-ApeKI"   +"\n"
                    +"  PstI-BfaI"    +"\n"
                    +"  PstI-EcoT22I" +"\n"
                    +"  PstI-MluCI"   +"\n"
                    +"  PstI-MseI"    +"\n"
                    +"  PstI-MspI"    +"\n"
                    +"  PstI-MspI-GDFcustom"+"\n"
                    +"  PstI-TaqI"    +"\n"
                    +"  SbfI-BfaI"    +"\n"
                    +"  SbfI-HpaII"   +"\n"
                    +"  SalI-MspI"    +"\n"
                    +"  SbfI-MspI"    +"\n"
                    +"  SexAI-Sau3AI" +"\n"
                    +"  SphI-EcoRI"   +"\n"
                    +"  StyI-MseI"    +"\n"
                    +"  ignore"    +"\n"
            );
            System.out.println("For two-enzyme digest, enzyme names should be separated by a dash, e.g. -e PstI-MspI");
            System.out.println("\nIf your enzyme is not on the above list you can use \"-e ignore\". In this case\n"
                    +"   barcodes and common adapter start sequences will be recognized, but chimeric DNA\n"
                    +"   fragments (or partial digests) will not be trimmed.");
    }

    public String enzyme() {
        return theEnzyme;
    }

    public String[] initialCutSiteRemnant() {
        return initialCutSiteRemnant;
    }

    public String[] likelyReadEnd() {
        return likelyReadEnd;
    }

    public int readEndCutSiteRemnantLength() {
        return readEndCutSiteRemnantLength;
    }
}

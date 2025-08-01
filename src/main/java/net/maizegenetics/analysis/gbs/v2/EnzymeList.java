package net.maizegenetics.analysis.gbs.v2;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.net.URISyntaxException;
import java.net.URL;
import java.util.LinkedHashMap;
import java.util.Set;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.ini4j.Ini;
import org.ini4j.Profile;

/**
 * This class contains a list of all the enzymes supported for GBS.
 * @author dpavelec
 */
public class EnzymeList {
    
    private static final Logger myLogger = LogManager.getLogger(EnzymeList.class);
    /**
     * This is the application level default that will initiate using the EnzymeList 
     * default constructor.
     */
    public static final EnzymeList defaultCache = new EnzymeList();
    
    /**
     * Looks for "enzymes.ini" config file in the lib path. 
     * Loads default Enzymes if not found 
     */
    public EnzymeList() {
        String parentPath = getJarPath(Ini.class).getPath();
        File iniFile = new File(parentPath + "/enzymes.ini");
        if (!loadFromFile(iniFile)) {
            loadDefaults();
        }
    }

    /**
     * Custom Config file
     * @param iniFile 
     */
    public EnzymeList(File iniFile) {
        if (!loadFromFile(iniFile)) {
            loadDefaults();
        }
    }

    private synchronized boolean loadFromFile(File iniFile) {
        if(map==null){
            map = new LinkedHashMap();
        }
        try {
            Ini ini = new Ini(iniFile);
            if (!iniFile.exists()) {
                return false;
            }
            ini.load();
            
            //List<Profile.Section> all = ini.getAll(ini);
            Set<String> keySet = ini.keySet();
            String empty = null;
            for(String key : keySet){
                Profile.Section section = (Profile.Section) ini.get(key);
                String name = section.get("name", empty);
                String initialCutSiteRemnant = section.get("initialCutSiteRemnant", empty);
                String likelyReadEnd = section.get("likelyReadEnd", empty);
                String readLenString = section.get("readEndCutSiteRemnantLength", "-1");
                int readEndCutSiteRemnantLength = Integer.parseInt(readLenString);
                if(name==null|initialCutSiteRemnant==null|likelyReadEnd==null|readEndCutSiteRemnantLength<0){
                    myLogger.warn("WARNING! Cannot load Enzyme section" + key);
                    continue;
                }else{
                    map.put(key.trim().toUpperCase(), new Enzyme(name, initialCutSiteRemnant.split(","), likelyReadEnd.split(","), readEndCutSiteRemnantLength));
                }
            }
        } catch (Exception ex) {
            myLogger.warn("ERROR! Cannot load Enzyme List -- " + iniFile.getPath());
            return false;
        }
        return true;
    }
    
    private void loadDefaults() {
        if(map==null){
            map = new LinkedHashMap();
        }
        map.put("APEKI", new Enzyme("ApeKI"
            ,new String[]{"CAGC", "CTGC"}
            ,new String[]{"GCAGC", "GCTGC", "GCAGAGAT", "GCTGAGAT"} // full cut site (from partial digest or chimera) or common adapter start
            ,4));
        map.put("PSTI", new Enzyme("PstI"
            ,new String[]{"TGCAG"}
            ,new String[]{"CTGCAG", "CTGCAAGAT"} // full cut site (from partial digest or chimera) or common adapter start
            ,5));
        map.put("ECOT22I", new Enzyme("EcoT22I"
            ,new String[]{"TGCAT"}
            ,new String[]{"ATGCAT", "ATGCAAGAT"} // full cut site (from partial digest or chimera) or common adapter start
            ,5));
        map.put("PASI", new Enzyme("PasI"
            ,new String[]{"CAGGG", "CTGGG"}
            ,new String[]{"CCCAGGG", "CCCTGGG", "CCCTGAGAT", "CCCAGAGAT"} // full cut site (from partial digest or chimera) or common adapter start
            ,5));
        map.put("HPAII", new Enzyme("HpaII"
            ,new String[]{"CGG"}
            ,new String[]{"CCGG", "CCGAGATCGG"} // full cut site (from partial digest or chimera) or common adapter start
            ,3));
        map.put("MSPI", new Enzyme("MspI"  // MspI and HpaII are isoschizomers (same recognition seq and overhang)
            ,new String[]{"CGG"}
            ,new String[]{"CCGG", "CCGAGATCGG"} // full cut site (from partial digest or chimera) or common adapter start
            ,3));
        map.put("PSTI-APEKI", new Enzyme("PstI-ApeKI"
            ,new String[]{"TGCAG"}
            ,new String[]{"GCAGC", "GCTGC", "CTGCAG", "GCAGAGAT", "GCTGAGAT"} // look for ApeKI site, PstI site, or common adapter for ApeKI
            ,4));
        map.put("PSTI-BFAI", new Enzyme("PstI-BfaI"  //  PstI: CTGCA^G  BfaI: C^TAG
            ,new String[]{"TGCAG"}
            ,new String[]{"CTGCAG", "CTAG", "CTAAGATC"} // look for PstI site, BfaI site, or common adapter start for BfaI
            ,3));
        map.put("PSTI-ECOT22I", new Enzyme("PstI-EcoT22I"
            ,new String[]{"TGCAG", "TGCAT"}
            ,new String[]{"ATGCAT", "CTGCAG", "CTGCAAGAT", "ATGCAAGAT"} // look for EcoT22I site, PstI site, or common adapter for PstI/EcoT22I
            ,5));
        map.put("PSTI-MSPI", new Enzyme("PstI-MspI"
            ,new String[]{"TGCAG"}
            ,new String[]{"CCGG", "CTGCAG", "CCGAGATC"} // look for MspI site, PstI site, or common adapter for MspI
            ,3));
        map.put("PSTI-MSPI-GDF-CUSTOM", new Enzyme("PstI-MspI-GDFcustom"
            ,new String[]{"TGCAG"}
            // changed from  CCGAGAT to CCGCTCAGG, as IGD/GDF used a custom Y adapter for MspI
            ,new String[]{"CCGG", "CTGCAG", "CCGCTCAGG"} // look for MspI site, PstI site, or GDF custom common adapter for MspI
            ,3));
        map.put("PSTI-TAQI", new Enzyme("PstI-TaqI"
            ,new String[]{"TGCAG"}
            ,new String[]{"TCGA", "CTGCAG", "TCGAGATC"} // look for TaqI site, PstI site, or common adapter for TaqI
            ,3));
        map.put("NSII-MSPI", new Enzyme("NsiI-MspI"  //  ATGCA^T   C^CGG
            ,new String[]{"TGCAT"}
            ,new String[]{"CCGG", "ATGCAT", "CCGAGATC"} // look for MspI site, NsiI site, or common adapter for MspI
            ,3));
        map.put("PAER7I-HHAI", new Enzyme("PaeR7I-HhaI"
            // Requested by Ye, Songqing, use same Y adapter as Polland paper  -QS
            ,new String[]{"TCGAG"}
            ,new String[]{"GCGC", "CTCGAG", "GCGAGATC"} // look for HhaI site, PaeR7I site, or common adapter for HhaI
            ,3));
        map.put("SBFI-MSPI", new Enzyme("SbfI-MspI"  // CCTGCA^GG  C^CGG
            ,new String[]{"TGCAGG"}
            ,new String[]{"CCGG", "CCTGCAGG", "CCGAGATC"} // look for MspI site, SbfI site, or common adapter for MspI
            ,3));
        map.put("SBFI-HPAII", new Enzyme("SbfI-HpaII"  // CCTGCA^GG  C^CGG  Nb: HpaII is an isoschizomer of MspI
            ,new String[]{"TGCAGG"}
            ,new String[]{"CCGG", "CCTGCAGG", "CCGAGATC"} // look for HpaII site, SbfI site, or common adapter for HpaII
            ,3));
        map.put("SBFI-BFAI", new Enzyme("SbfI-BfaI"  // CCTGCA^GG  C^TAG
            ,new String[]{"TGCAGG"}
            ,new String[]{"CTAG", "CCTGCAGG", "CTAAGATC"} // look for BfaI site, SbfI site, or common adapter for BfaI
            ,3));
        map.put("SPHI-ECORI", new Enzyme("SphI-EcoRI"  // GCATG^C  G^AATTC
            ,new String[]{"CATGC"}  // SphI overhang from the genomic DNA fragment (top strand)
            ,new String[]{"GCATGC", "GAATTC", "GAATTAGATC"} // look for SphI site, EcoRI site, or common adapter for EcoRI
            ,5));
        map.put("ASISI-MSPI", new Enzyme("AsiSI-MspI"
            ,new String[]{"ATCGC"}
            ,new String[]{"CCGG", "GCGATCGC", "CCGAGATC"} // look for MspI site, AsiSI site, or common adapter for MspI
            ,3));
        map.put("BSSHII-MSPI", new Enzyme("BssHII-MspI"
            ,new String[]{"CGCGC"}
            ,new String[]{"CCGG", "GCGCGC", "CCGAGATC"} // look for MspI site, BssHII site, or common adapter for MspI
            ,3));
        map.put("FSEI-MSPI", new Enzyme("FseI-MspI"
            ,new String[]{"CCGGCC"}
            ,new String[]{"CCGG", "GGCCGGCC", "CCGAGATC"} // look for MspI site, FseI site, or common adapter for MspI
            ,3));
        map.put("SALI-MSPI", new Enzyme("SalI-MspI"
            ,new String[]{"TCGAC"}
            ,new String[]{"CCGG", "GTCGAC", "CCGAGATC"} // look for MspI site, SalI site, or common adapter for MspI
            ,3));
        map.put("ECORI-MSPI", new Enzyme("EcoRI-MspI"   //  G^AATTC  C^CGG
            ,new String[]{"AATTC"}
            ,new String[]{"CCGG", "GAATTC", "CCGAGATC"} // look for MspI site, EcoRI site, or Poland et al. 2012 Y-adapter for MspI
            ,3));
        map.put("HINDIII-MSPI", new Enzyme("HindIII-MspI" // A^AGCTT   C^CGG
            ,new String[]{"AGCTT"}
            ,new String[]{"CCGG", "AAGCTT", "CCGAGATC"} // look for MspI site, HindIII site, or Poland et al. 2012 Y-adapter for MspI
            ,3));
        map.put("HINDIII-NLAIII", new Enzyme("HindIII-NlaIII" // A^AGCTT   ^CATG (not blunt)
            ,new String[]{"AGCTT"}
            ,new String[]{"CATG", "AAGCTT", "CATGAGATC"} // look for NlaIII site, HindIII site, or Poland et al. 2012 Y-adapter for Nla3
            ,4));
        map.put("SEXAI-SAU3AI", new Enzyme("SexAI-Sau3AI"  // A^CCWGGT   ^GATC (not blunt)
            ,new String[]{"CCAGGT", "CCTGGT"}
            ,new String[]{"GATC", "ACCAGGT", "ACCTGGT", "GATCAGATC"} // look for SexAI site, Sau3AI site, or Poland et al. 2012 Y-adapter for Sau3AI
            ,4));
        map.put("BAMHI-MLUCI", new Enzyme("BamHI-MluCI"  // G^GATCC   ^AATT (not blunt)
            ,new String[]{"GATCC"}
            ,new String[]{"AATT", "GGATCC", "AATTAGATC"} // look for MluCI site, BamHI site, or Poland et al. 2012 Y-adapter for MluCI
            ,4));
        map.put("PSTI-MLUCI", new Enzyme("PstI-MluCI"  // CTGCA^G   ^AATT (not blunt)
            ,new String[]{"TGCAG"}
            ,new String[]{"AATT", "CTGCAG", "AATTAGATC"} // look for MluCI site, PstI site, or Poland et al. 2012 Y-adapter for MluCI
            ,4));
        map.put("PSTI-MSEI", new Enzyme("PstI-MseI" // CTGCA^G   T^TAA
            ,new String[]{"TGCAG"}
            ,new String[]{"TTAA", "CTGCAG", "TTAAGATC"} // look for MseI site, PstI site, or Poland et al. 2012 Y-adapter for MseI
            ,3));
        map.put("AVAII-MSEI", new Enzyme("AvaII-MseI" // G^GWCC   T^TAA  W=AorT
            ,new String[]{"GACC", "GTCC"}
            ,new String[]{"TTAA", "GGACC", "GGTCC", "TTAAGATC"} // look for MseI site, AvaII site, or Poland et al. 2012 Y-adapter for MseI
            ,3));
        map.put("ECORI-MSEI", new Enzyme("EcoRI-MseI" // G^AATTC   T^TAA
            ,new String[]{"AATTC"}
            ,new String[]{"TTAA", "GAATTC", "TTAAGATC"} // look for MseI site, EcoRI site, or Poland et al. 2012 Y-adapter for MseI
            ,3));
        map.put("ECORI-AVAII", new Enzyme("EcoRI-AvaII" // G^AATTC   G^GWCC
            ,new String[]{"AATTC"}
            ,new String[]{"GGACC", "GGTCC", "GAATTC", "GGACAGATC", "GGTCAGATC"} // look for AvaII site, EcoRI site, or Poland et al. 2012 Y-adapter for AvaII
            ,4));
        map.put("ECORI-HINFI", new Enzyme("EcoRI-HinfI" // G^AATTC   G^ANTC
            ,new String[]{"AATTC"}
            ,new String[]{"GAATC", "GACTC", "GAGTC", "GATTC", "GAATTC", "GAATAGATC", "GACTAGATC", "GAGTAGATC", "GATTAGATC"} // look for HinfI site, EcoRI site, or Poland et al. 2012 Y-adapter for HinfI
            ,4));
        map.put("BBVCI-MSPI", new Enzyme("BbvCI-MspI" // CCTCAGC (-5/-2)   C^CGG
            ,new String[]{"TCAGC"}
            ,new String[]{"CCGG", "CCTCAGC", "CCGAGATC"} // look for MspI site, BbvCI site, or Poland et al. 2012 Y-adapter for MspI
            ,3));
        map.put("MSPI-APEKI", new Enzyme("MspI-ApeKI"  // C^CGG  G^CWGC
            ,new String[]{"CGG", "CAGC", "CTGC"}
            ,new String[]{"CCGG", "GCAGC", "GCTGC", "CCGAGATCGG", "GCAGAGAT", "GCTGAGAT"} 
            ,3 ));
        map.put("APOI", new Enzyme("ApoI"
            ,new String[]{"AATTC","AATTT"}
            ,new String[]{"AAATTC","AAATTT","GAATTC","GAATTT","AAATTAGAT","GAATTAGAT"} // full cut site (from partial digest or chimera) or common adapter start
            ,5));
        map.put("BAMHI", new Enzyme("BamHI"
            ,new String[]{"GATCC"}
            ,new String[]{"GGATCC", "GGATCAGAT"} // full cut site (from partial digest or chimera) or common adapter start
            // ,new String[]{"GGATCC", "AGATCGGAA", "AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG"} // <-- corrected from this by Jeff Glaubitz on 2012/09/12
            ,5));
        map.put("MSEI", new Enzyme("MseI"
            ,new String[]{"TAA"}
            ,new String[]{"TTAA", "TTAAGAT"} // full cut site (from partial digest or chimera) or common adapter start
            ,3));
        map.put("SAU3AI", new Enzyme("Sau3AI"
            ,new String[]{"GATC"}
            ,new String[]{"GATC","GATCAGAT"} // full cut site (from partial digest or chimera) or common adapter start
            ,4));
        map.put("NDEI", new Enzyme("NdeI"
            ,new String[]{"TATG"}
            ,new String[]{"CATATG","CATAAGAT"} // full cut site (from partial digest or chimera) or common adapter start
            ,4));
        map.put("HINP1I", new Enzyme("HinP1I"
            ,new String[]{"CGC"}
            ,new String[]{"GCGC","GCGAGAT"} // full cut site (from partial digest or chimera) or common adapter start
            ,3));
        map.put("SBFI", new Enzyme("SbfI"
            ,new String[]{"TGCAGG"}
            ,new String[]{"CCTGCAGG","CCTGCAAGAT"} // full cut site (from partial digest or chimera) or common adapter start
            ,6));
        map.put("HINDIII", new Enzyme("HindIII" // A^AGCTT
            ,new String[]{"AGCTT"}
            ,new String[]{"AAGCTT","AAGCTAGAT"} // full cut site (from partial digest or chimera) or common adapter start
            ,5));
        map.put("ECORI", new Enzyme("EcoRI"  // G^AATTC
            ,new String[]{"AATTC"}
            ,new String[]{"GAATTC","GAATTAGAT"} // full cut site (from partial digest or chimera) or common adapter start
            ,5));
        map.put("CVI1I", new Enzyme("CviQI"  // CviQI and Csp6I are isoschizomers (same recognition seq and overhang)
            ,new String[]{"TAC"}
            ,new String[]{"GTAC","GTAAGATCGG"} // full cut site (from partial digest or chimera) or common adapter start
            ,3));
        map.put("CSP6I", new Enzyme("Csp6I"  // Csp6I and CviQI are isoschizomers (same recognition seq and overhang)
            ,new String[]{"TAC"}
            ,new String[]{"GTAC","GTAAGATCGG"} // full cut site (from partial digest or chimera) or common adapter start
            ,3));
        map.put("NLAIII", new Enzyme("NlaIII" // CATG^
            ,new String[]{"CATG"}
            ,new String[]{"CATG","CATGAGAT"} // full cut site (from partial digest or chimera) or common adapter start
            ,4));
        map.put("SPHI", new Enzyme("SphI"  // GCATG^C
            ,new String[]{"CATGC"}
            ,new String[]{"GCATGC","GCATGAGAT"} // full cut site (from partial digest or chimera) or common adapter start
            ,5));
        map.put("NSPI", new Enzyme("NspI"  // RCATG^Y
            ,new String[]{"CATGC","CATGT"}
            ,new String[]{"ACATGT","GCATGC","ACATGAGAT","GCATGAGAT"} // full cut site (from partial digest or chimera) or common adapter start
            ,5));
        map.put("KPNI", new Enzyme("KpnI"  // GGTAC^C
            ,new String[]{"GTACC"}
            ,new String[]{"GGTACC","GGTACAGAT"} // full cut site (from partial digest or chimera) or common adapter start
            ,5));
        map.put("STYI", new Enzyme("StyI"  // C^CWWGG
            ,new String[]{"CAAGG","CATGG","CTAGG","CTTGG"}
            ,new String[]{"CCAAGG","CCATGG","CCTAGG","CCTTGG","CCAAGAGAT","CCATGAGAT","CCTAGAGAT","CCTTGAGAT"} // full cut site (from partial digest or chimera) or common adapter start
            ,5));
        map.put("STYI-MSEI", new Enzyme("StyI-MseI"  // C^CWWGG & T^TAA
            ,new String[]{"CAAGG","CATGG","CTAGG","CTTGG"}
            ,new String[]{"TTAA","CCAAGG","CCATGG","CCTAGG","CCTTGG","TTAAGAT"} // full cut site (from partial digest or chimera) or common adapter start
            ,3));
        map.put("FSEI", new Enzyme("FseI"  // GGCCGG^CC
            ,new String[]{"CCGGCC"}
            ,new String[]{"GGCCGGCC","AGATCGGAAG"} // full cut site (from partial digest or chimera) or Morishige et al (BMC Genomics, 2013) T adapter start
            ,0));  // assumes that common T adapter is far more likely than a second full cut site));
        map.put("NgoMIV", new Enzyme("NgoMIV"  // G^CCGGC
            ,new String[]{"CCGGC"}
            ,new String[]{"GCCGGC","AGATCGGAAG"} // full cut site (from partial digest or chimera) or Morishige et al (BMC Genomics, 2013) T adapter start
            ,0));  // assumes that common T adapter is far more likely than a second full cut site));
        map.put("MSLI", new Enzyme("MslI"  // CAYNN^NNRTG  -- has 32 different cut sites (assuming constrained to palindromic YNN -- 32^2 otherwise)
            ,new String[]{""}
            ,new String[]{"AGATCGGA"} // common adapter start only (too many possible cut sites!)
            ,0  ));
        map.put("ASEI", new Enzyme("AseI"  // AT^TAAT
            ,new String[]{"TAAT"}
            ,new String[]{"ATTAAT","ATTAAGAT"}  // full cut site (from partial digest or chimera) or common adapter start
            ,4  ));
        map.put("AVAII", new Enzyme("AvaII"  // G^GWCC
            ,new String[]{"GACC","GTCC"}
            ,new String[]{"GGACC","GGTCC","GGACAGAT","GGTCAGAT"}  // full cut site (from partial digest or chimera) or common adapter start
            ,4  ));
        map.put("KPNI-MSPI", new Enzyme("KpnI-MspI"  //  KpnI: GGTAC^C  MspI: C^CGG
            ,new String[]{"GTACC"}
            ,new String[]{"CCGG", "GGTACC", "CCGAGATC"} // look for MspI site, KpnI site, or common adapter start for MspI
            ,3));
        map.put("RBSTA", new Enzyme("RBSTA"
            ,new String[]{"TA"}
            ,new String[]{"TTAA", "GTAC", "CTAG", "TTAAGAT", "GTAAGAT", "CTAAGAT"} // full cut site (from partial digest or chimera) of MseI, CVIQi, XspI or common adapter start
            ,3));
        map.put("RBSCG", new Enzyme("RBSCG"
            ,new String[]{"CG"}
            ,new String[]{"CCGC", "TCGA", "GCGC", "CCGG", "ACGT", "CCGAGAT", "TCGAGAT", "GCGAGAT", "ACGAGAT"} // full cut site (from partial digest or chimera) of AciI, TaqaI, HinpI, HpaII, HpyCH4IV or common adapter start
            ,3));
        map.put("IGNORE", new Enzyme("unspecified"  // can be used for new enzymes -- only looks for barcodes and common adapter starts
            ,new String[]{""}
            ,new String[]{"AGATCGGA"} // common adapter start only
            ,0));
    }

    /**
     * Gets the file location of the jar file containing the specified class.
     * <br>NOTE: may not work with remote files.
     *
     * @param clazz class to find
     * @return the location of the jar
     */
    static File getJarPath(Class clazz) {
        try {
            URL location = clazz.getProtectionDomain().getCodeSource().getLocation();
            File f = new File(location.toURI().getPath());
            if (!f.exists()) {
                File nf = new File(location.getPath());
                if (nf.exists()) {
                    f = nf;
                }
            }
            if (f.getName().toLowerCase().endsWith(".jar")) {
                return f.getParentFile();
            }
            return f;
        } catch (URISyntaxException ex) {
            return new File(clazz.getProtectionDomain().getCodeSource().getLocation().getPath());
        }
    }

    /**
     * Get the Enzyme by the section name <br>
     * Replaces the previous new GBSEnzyme(String enzyme)
     * @param key The section name specified the enzyme configuration file
     * @return Enzyme
     */
    public synchronized Enzyme getEnzyme(String key) {
        if(map==null){
            loadDefaults();
        }
        return map.get(key.trim().toUpperCase());
    }
    
    /**
     * Add Enzyme to a configuration file
     * @param ini the config file 
     * @param enzyme The enzyme to add
     * @param key The name of the key/section. Uses Enzyme.name if Empty or null.
     */
    static void add(Ini ini,Enzyme enzyme, String key) {
        //key should be upper case and no leading or trailing whitespace
        if(key == null){
            key = enzyme.name.trim().toUpperCase();
        }else{
            key = key.trim().toUpperCase();
        }
        Profile.Section section = ini.containsKey(key) ? (Profile.Section) ini.get(key) : ini.add(key);
        section.put("name", enzyme.name);
        //initialCutSiteRemnant
        StringBuilder builder = new StringBuilder(enzyme.initialCutSiteRemnant[0]);
        for (int i = 1; i < enzyme.initialCutSiteRemnant.length; ++i) {
            builder.append(",");
            builder.append(enzyme.initialCutSiteRemnant[i]);
        }
        section.put("initialCutSiteRemnant", builder.toString());
        //likelyReadEnd
        builder = new StringBuilder(enzyme.likelyReadEnd[0]);
        for (int i = 1; i < enzyme.likelyReadEnd.length; ++i) {
            builder.append(",");
            builder.append(enzyme.likelyReadEnd[i]);
        }
        section.put("likelyReadEnd", builder.toString());
        section.put("readEndCutSiteRemnantLength", enzyme.readEndCutSiteRemnantLength);
    }
    
    /**
     * Prints the list of Enzymes loaded in this object
     * @param print 
     * @throws IOException 
     */
    public void printEnzymes(PrintStream print) throws IOException{
        Ini ini = new Ini(); //using the Ini object format the output
        ini.getConfig().setLineSeparator("\n");
        if(this.map==null){
            loadDefaults();
        }
        for(String key : this.map.keySet()){
            Enzyme enzyme = this.map.get(key);
            add(ini,enzyme,key);
        }
        ini.store(print);
    }
    
    /**
     * Prints default enzyme list and exists
     * @param args
     * @throws IOException 
     */
    public static void main(String[] args) throws IOException{
        EnzymeList.defaultCache.printEnzymes(System.out);
    }

    private LinkedHashMap<String, Enzyme> map;
    /**
     * Replaces the previous GBSEnzyme class
     */
    public static class Enzyme {
        private Enzyme(
                String name,
                String[] initialCutSiteRemnant,
                String[] likelyReadEnd,
                int readEndCutSiteRemnantLength) {
            this.name = name;
            this.initialCutSiteRemnant = initialCutSiteRemnant;
            this.likelyReadEnd = likelyReadEnd;
            this.readEndCutSiteRemnantLength = readEndCutSiteRemnantLength;
        }
        public final String name;
        public final String[] initialCutSiteRemnant;
        public final String[] likelyReadEnd;
        public final int readEndCutSiteRemnantLength;
        
        /** 
         * This methods is for convenience based on old code. Use Enzyme.name
         * @deprecated 
         */
        public String enzyme() {return name;}
        /** 
         * This methods is for convenience based on old code. Use Enzyme.initialCutSiteRemnant
         * @deprecated 
         */
        public String[] initialCutSiteRemnant() {return initialCutSiteRemnant;}
        /** 
         * This methods is for convenience based on old code. Use Enzyme.likelyReadEnd
         * @deprecated 
         */
        public String[] likelyReadEnd() {return likelyReadEnd;}
        /** 
         * This methods is for convenience based on old code. Use Enzyme.readEndCutSiteRemnantLength
         * @deprecated 
         */
        public int readEndCutSiteRemnantLength() {return readEndCutSiteRemnantLength;}
    }
}

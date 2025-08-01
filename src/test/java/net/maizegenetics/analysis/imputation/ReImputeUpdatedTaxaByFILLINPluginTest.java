package net.maizegenetics.analysis.imputation;

import ch.systemsx.cisd.hdf5.HDF5Factory;
import ch.systemsx.cisd.hdf5.IHDF5Reader;
import ch.systemsx.cisd.hdf5.IHDF5Writer;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.HashMap;
import java.util.Map;
import java.util.Map.Entry;
import java.util.TreeMap;
import net.maizegenetics.constants.GBSConstants;
import net.maizegenetics.dna.map.PositionListBuilder;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableBuilder;
import net.maizegenetics.dna.snp.genotypecall.BasicGenotypeMergeRule;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.util.HDF5Utils;
import net.maizegenetics.util.Utils;
import org.junit.*;

import static org.junit.Assert.*;

/**
 *
 * @author jcg233
 */
public class ReImputeUpdatedTaxaByFILLINPluginTest {
    
    Map<String, Double> currentImpCallRates = new HashMap();
    Map<String, String> additionalDepthTaxa = new TreeMap(); // key = targetTaxon; value = sourceTaxon
    Map<String, String[]> taxaToMerge = new TreeMap(); // key = new taxon name; value = array of two source taxon names
    String tempInitialImpFile, tempRawFile, modRawFile;
    
    public ReImputeUpdatedTaxaByFILLINPluginTest() {
    }
    
    @Before
    public void setUp() {  // note that this is run before EVERY test
        File outdir = new File(GBSConstants.GBS_TEMP_REIMPUTE_UPDATED_TAXA_BY_FILLIN_PLUGIN_DIR);
        outdir.mkdirs();
        tempInitialImpFile = Utils.getFilename(GBSConstants.GBS_TEMP_REIMPUTE_UPDATED_TAXA_BY_FILLIN_PLUGIN_IMP_FILE)+".h5";
        String tempInitialImpFile2 = Utils.getFilename(GBSConstants.GBS_TEMP_REIMPUTE_UPDATED_TAXA_BY_FILLIN_PLUGIN_IMP_FILE2)+".h5";
        tempRawFile = Utils.getFilename(GBSConstants.GBS_TEMP_REIMPUTE_UPDATED_TAXA_BY_FILLIN_PLUGIN_RAW_FILE)+".h5";
        modRawFile = Utils.getFilename(GBSConstants.GBS_TEMP_REIMPUTE_UPDATED_TAXA_BY_FILLIN_PLUGIN_MOD_FILE)+".h5";
        if(outdir.listFiles()!=null) {
            for (File file : outdir.listFiles()) {
                if(file.getName().equals(tempInitialImpFile) || file.getName().equals(tempInitialImpFile2) 
                        || file.getName().equals(tempRawFile) || file.getName().equals(modRawFile)) {
                    file.delete();
                }
            }
        }
        // reassign these to their full paths
        tempInitialImpFile = GBSConstants.GBS_TEMP_REIMPUTE_UPDATED_TAXA_BY_FILLIN_PLUGIN_IMP_FILE;
        tempInitialImpFile2 = GBSConstants.GBS_TEMP_REIMPUTE_UPDATED_TAXA_BY_FILLIN_PLUGIN_IMP_FILE2;
        tempRawFile = GBSConstants.GBS_TEMP_REIMPUTE_UPDATED_TAXA_BY_FILLIN_PLUGIN_RAW_FILE;
        modRawFile = GBSConstants.GBS_TEMP_REIMPUTE_UPDATED_TAXA_BY_FILLIN_PLUGIN_MOD_FILE;
        
        File rawGenosFile = new File(GBSConstants.GBS_EXPECTED_ADD_REFERENCE_ALLELE_TO_HDF5_PLUGIN_OUT_FILE);
        try {
            Files.copy(rawGenosFile.toPath(), new File(tempRawFile).toPath());
            Files.copy(rawGenosFile.toPath(), new File(modRawFile).toPath());
        } catch (IOException e) {
            Assert.fail("Test failed: " + e.getMessage());
        }
        
        // generate an initial imputed genotypes file
        FILLINImputationPlugin fip = new FILLINImputationPlugin()
            .targetFile(tempRawFile)
            .outputFilename(tempInitialImpFile)
            .donorDir(GBSConstants.GBS_FILLIN_DONOR_DIR)
            .preferredHaplotypeSize(4000)
        ;
        fip.performFunction(null);
        
        // create a copy of the initial imputed file in an unlocked, "closeUnfinished" state
        GenotypeTableBuilder gtb = GenotypeTableBuilder.getTaxaIncremental(PositionListBuilder.getInstance(tempInitialImpFile),tempInitialImpFile2);
        IHDF5Reader h5ImpReader = HDF5Factory.openForReading(tempInitialImpFile);
        for (String tax : HDF5Utils.getAllTaxaNames(h5ImpReader)) {
            gtb.addTaxon(HDF5Utils.getTaxon(h5ImpReader, tax), HDF5Utils.getHDF5GenotypesCalls(h5ImpReader, tax));
        }
        gtb.closeUnfinished();
        tempInitialImpFile = tempInitialImpFile2;
    }
        
    private void modifyTestRawGenosFile() {
        // unlock the raw file for taxa modifications
        IHDF5Writer h5Writer = HDF5Factory.open(modRawFile);
        HDF5Utils.unlockHDF5TaxaModule(h5Writer);
        HDF5Utils.unlockHDF5GenotypeModule(h5Writer);
        h5Writer.close();

        // Add depth to some existing taxa in the original raw genotypes file
        // The source taxa (of the additional depth) were chosen from a NJ tree based on their high degree of IBS with the target taxa
        additionalDepthTaxa.put("A680:250047931", "B73Htrhm:250047925"); // key = targetTaxon; value = sourceTaxon
        additionalDepthTaxa.put("M0061:250021010", "M0063:250021011");
        additionalDepthTaxa.put("NC318:250047989", "NC320:250047967");
        additionalDepthTaxa.put("CM105:250047997", "CM174:250047927");        
        IHDF5Reader h5Reader = HDF5Factory.openForReading(tempRawFile);
        IHDF5Reader h5ImpReader = HDF5Factory.openForReading(tempInitialImpFile);
        BasicGenotypeMergeRule genoMergeRule = new BasicGenotypeMergeRule(0.01);
        GenotypeTableBuilder modGenos = GenotypeTableBuilder.mergeTaxaIncremental(modRawFile, genoMergeRule);
        for (Entry<String,String> depthTaxaPair : additionalDepthTaxa.entrySet()) {
            String targetTaxonName = depthTaxaPair.getKey();
            String sourceTaxonName = depthTaxaPair.getValue();
            currentImpCallRates.put(targetTaxonName, calculateCallRate(HDF5Utils.getHDF5GenotypesCalls(h5ImpReader, targetTaxonName)));
            Taxon modTaxon = new Taxon.Builder(targetTaxonName).addAnno("Flowcell_Lane", "FAKEFLOW_9").build();
            byte[] additionalGenos = HDF5Utils.getHDF5GenotypesCalls(h5Reader, sourceTaxonName);
            byte[][] additionalDepth = HDF5Utils.getHDF5GenotypesDepth(h5Reader, sourceTaxonName);
            modGenos.addTaxon(modTaxon, additionalGenos, additionalDepth);
        }
        h5ImpReader.close();
        
        // Add new taxa by merging depths from some genetically similar taxa from original raw genotypes file, giving the new taxa new names
        // The taxa to merge to create a new taxon were chosen from a NJ tree based on their high degree of IBS
        taxaToMerge.put("B73-M0382hyb:1001", new String[] {"B73:250020986", "M0382:250021080"}); // key = new taxon name; value = array of two source taxon names
        taxaToMerge.put("Mo17-M0326hyb:1002", new String[] {"Mo17:250020994", "M0326:250021071"});
        taxaToMerge.put("M0368-M0369hyb:1003", new String[] {"M0368:250021078", "M0369:250021067"});
        taxaToMerge.put("NC262-NC290Ahyb:1004", new String[] {"NC262:250047948", "NC290A:250047937"});
        for (Entry<String,String[]> mergeSet: taxaToMerge.entrySet()) {
            int counter = 0;
            for (String sourceTaxon : mergeSet.getValue()) {
                String newTaxonName = mergeSet.getKey();
                Taxon newTaxon = new Taxon.Builder(newTaxonName).addAnno("Flowcell_Lane", "FAKEFLOW_"+counter).build();
                byte[] genos = HDF5Utils.getHDF5GenotypesCalls(h5Reader, sourceTaxon);
                byte[][] depths = HDF5Utils.getHDF5GenotypesDepth(h5Reader, sourceTaxon);
                modGenos.addTaxon(newTaxon, genos, depths);
                counter++;
            }
        }
        h5Reader.close();
        modGenos.closeUnfinished();
    }
    
    @After
    public void tearDown() {
    }

    /**
     * Test of processData method, of class ReImputeUpdatedTaxaByFILLINPlugin.
     */
    @Ignore
    @Test
    public void testProcessData() {

        modifyTestRawGenosFile();  // add 4 new taxa, and additional depth to 4 other taxa
        
        // run the ReImputeUpdatedTaxaByFILLINPlugin
        ReImputeUpdatedTaxaByFILLINPlugin ReImpPlugin = new ReImputeUpdatedTaxaByFILLINPlugin()
            .rawHDF5GenotypeFile(modRawFile)
            .imputedHDF5GenotypeFile(tempInitialImpFile)
            .donorDir(GBSConstants.GBS_FILLIN_DONOR_DIR)
            .preferredHaplotypeSize(4000)
        ;
        ReImpPlugin.performFunction(null);

        // open the updated imputed file
        IHDF5Reader impGenosReader=HDF5Factory.openForReading(tempInitialImpFile);
        
        // test that the added taxa (= merger of two existing taxa) exist and were imputed
        assertEquals("Wrong number of taxa in the updated imputed file: ", 196, HDF5Utils.getAllTaxaNames(impGenosReader).size());
        for (String addedTaxonName : taxaToMerge.keySet()) {
            Taxon addedTaxon = HDF5Utils.getTaxon(impGenosReader, addedTaxonName);
            assertTrue("The taxon "+addedTaxonName+" was not added", addedTaxon != null);
            String[] fcLanes = addedTaxon.getAnnotation().getTextAnnotation("Flowcell_Lane");
            assertTrue("The number of Flowcell_Lane annotations for added taxon "+addedTaxonName+" does not equal 2", fcLanes.length == 2);
            for (String fcLane : fcLanes) {
                assertTrue("Flowcell_Lane annotation for added taxon "+addedTaxonName+" does not contain FAKEFLOW_", fcLane.contains("FAKEFLOW_"));
            }
            assertTrue("Genotypes were not added for the taxon: "+addedTaxonName, HDF5Utils.doTaxonCallsExist(impGenosReader, addedTaxonName));
            assertTrue("Call rate is <65% for imputed added taxon: "+addedTaxonName, 
                    calculateCallRate(HDF5Utils.getHDF5GenotypesCalls(impGenosReader, addedTaxonName)) >= 0.65);
        }
        
        // test that the taxa with additional depth ("augmented taxa") were reimputed
        for (String augmentedTaxonName : additionalDepthTaxa.keySet()) {
            Taxon augTaxon = HDF5Utils.getTaxon(impGenosReader, augmentedTaxonName);
            assertTrue("The augmented taxon "+augmentedTaxonName+" is not present", augTaxon != null);
            String[] fcLanes = augTaxon.getAnnotation().getTextAnnotation("Flowcell_Lane");
            assertTrue("The number of Flowcell_Lane annotations for augmented taxon "+augmentedTaxonName+" does not equal 2", fcLanes.length == 2);
            boolean found = false;
            for (String fcLane : fcLanes) {
                if (fcLane.equals("FAKEFLOW_9")) found = true;
            }
            assertTrue("Flowcell_Lane annotations for added taxon "+augmentedTaxonName+" does not contain FAKEFLOW_9", found);
            double newCallRate = calculateCallRate(HDF5Utils.getHDF5GenotypesCalls(impGenosReader, augmentedTaxonName));
            double oldCallRate = currentImpCallRates.get(augmentedTaxonName);
            assertTrue("Imputed SNP call rate < 97.5% of original imputed call rate for augmented taxon "
                    +augmentedTaxonName+": newCallRate:"+newCallRate+"  oldCallRate:"+oldCallRate, (newCallRate/oldCallRate) >= 0.975);
        }
        
    }
    
    private double calculateCallRate(byte[] genos) {
        if (genos.length == 0) {
            return Double.NaN;
        }
        int nonMissing = 0, total = 0;
        for (byte geno : genos) {
            if(geno != GenotypeTable.UNKNOWN_DIPLOID_ALLELE) {
                nonMissing++;
            }
            total++;
        }
        return (double) nonMissing/total;
    }
    
}

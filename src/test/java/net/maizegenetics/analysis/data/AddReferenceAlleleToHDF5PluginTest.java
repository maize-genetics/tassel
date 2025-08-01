/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.analysis.data;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;

import net.maizegenetics.constants.GBSConstants;
import net.maizegenetics.dna.WHICH_ALLELE;
import net.maizegenetics.dna.map.Position;
import net.maizegenetics.dna.map.PositionList;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.ImportUtils;

import org.junit.Before;
import org.junit.Test;

import static org.junit.Assert.*;

/**
 *
 * @author jcg233
 */
public class AddReferenceAlleleToHDF5PluginTest {
    
    public AddReferenceAlleleToHDF5PluginTest() {
    }
    
    @Before
    public void setUp() {
    }

    /**
     * Test of performFunction method, of class AddReferenceAlleleToHDF5Plugin.
     */
    @Test
    public void testPerformFunction() {
        File theDir = new File(GBSConstants.GBS_TEMP_ADD_REFERENCE_ALLELE_TO_HDF5_PLUGIN_DIR);
        if (!theDir.exists()) { // mkdirs() below returns false if directory already exists
            if (!(new File(GBSConstants.GBS_TEMP_ADD_REFERENCE_ALLELE_TO_HDF5_PLUGIN_DIR)).mkdirs()) {
                String absoluteDir = new File(GBSConstants.GBS_TEMP_ADD_REFERENCE_ALLELE_TO_HDF5_PLUGIN_DIR).getAbsolutePath();
                throw new IllegalStateException("AddReferenceAlleleToHDF5PluginTest: testPerformFunction: Can't create output directory: " + absoluteDir);
            }
        }
        try{
            Files.deleteIfExists(Paths.get(GBSConstants.GBS_TEMP_ADD_REFERENCE_ALLELE_TO_HDF5_PLUGIN_OUT_FILE));
        } catch (IOException e) {
            e.printStackTrace();
        }
        AddReferenceAlleleToHDF5Plugin instance = new AddReferenceAlleleToHDF5Plugin()
            .inputHDF5GenotypeFile(GBSConstants.GBS_EXPECTED_PRODUCTION_SNP_CALLER_PLUGIN_HDF5_OUT_FILE)
            .referenceGenomeFile(GBSConstants.GBS_REFERENCE_GENOME)
            .referenceGenomeVersion("B73 RefGenV2")
            .outputHDF5GenotypeFile(GBSConstants.GBS_TEMP_ADD_REFERENCE_ALLELE_TO_HDF5_PLUGIN_OUT_FILE);
        instance.performFunction(null);
        GenotypeTable result = ImportUtils.readGuessFormat(GBSConstants.GBS_TEMP_ADD_REFERENCE_ALLELE_TO_HDF5_PLUGIN_OUT_FILE);
        PositionList posits = result.positions();
        int nSitesRefAlleleDifferent = 0;
        int nSites = 0;
        for (Position pos : posits) {
            byte ref = pos.getAllele(WHICH_ALLELE.Reference);
            byte maj = pos.getAllele(WHICH_ALLELE.GlobalMajor);
            byte min = pos.getAllele(WHICH_ALLELE.GlobalMinor);
            if (ref != maj && ref != min) {
                nSitesRefAlleleDifferent++;
            }
            nSites++;
        }
        double propRefDifferent = (double) nSitesRefAlleleDifferent / nSites;
        double epsilon = 0.006;
        assertEquals("Proportion of sites with reference alleles that differ from either the major or minor allele is greater than "+epsilon, 0.0, propRefDifferent, epsilon);
    }

}

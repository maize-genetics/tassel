/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.analysis.data;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Paths;
import static java.nio.file.StandardCopyOption.REPLACE_EXISTING;
import java.util.Iterator;
import net.maizegenetics.constants.GBSConstants;
import net.maizegenetics.constants.GeneralConstants;
import net.maizegenetics.dna.map.Position;
import net.maizegenetics.dna.map.PositionList;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.ImportUtils;
import net.maizegenetics.plugindef.DataSet;
import org.junit.Before;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author jcg233
 */
public class SplitHDF5ByChromosomePluginTest {
    String tempInFile;
    
    public SplitHDF5ByChromosomePluginTest() {
    }
    
    @Before
    public void setUp() {
        String source = GBSConstants.GBS_EXPECTED_ADD_REFERENCE_ALLELE_TO_HDF5_PLUGIN_OUT_FILE;
        tempInFile = GeneralConstants.TEMP_DIR+new File(GBSConstants.GBS_EXPECTED_ADD_REFERENCE_ALLELE_TO_HDF5_PLUGIN_OUT_FILE).getName();
        try {
            Files.copy(Paths.get(source), Paths.get(tempInFile), REPLACE_EXISTING);
        } catch (Exception e) {
            throw new IllegalStateException("SplitHDF5ByChromosomePluginTest: setUp: Can't copy input file: " + e);
        }
    }

    /**
     * Test of performFunction method, of class SplitHDF5ByChromosomePlugin.
     */
    @Test
    public void testPerformFunction() {
        SplitHDF5ByChromosomePlugin instance = new SplitHDF5ByChromosomePlugin()
            .inputHDF5GenotypeFile(tempInFile);
        instance.performFunction(null);

        GenotypeTable orig = ImportUtils.readGuessFormat(GBSConstants.GBS_EXPECTED_ADD_REFERENCE_ALLELE_TO_HDF5_PLUGIN_OUT_FILE);
        PositionList OrigPosits = orig.positions();
        Iterator<Position> OrigPositsIter = OrigPosits.iterator();
        
        int totalNSites = 0;
        for (int chr = 9; chr < 11; chr++) {
            String splitFile = tempInFile.replaceFirst("\\.h5$", "_chr"+chr+".h5");
            GenotypeTable split = ImportUtils.readGuessFormat(splitFile);
            totalNSites += split.numberOfSites();
            assertEquals("Output file "+new File(splitFile).getName()+" contains a different number of taxa than input", orig.numberOfTaxa(), split.numberOfTaxa());
            PositionList SplitPosits = split.positions();
            Iterator<Position> SplitPositsIter = SplitPosits.iterator();
            while (SplitPositsIter.hasNext()) {
                Position splitPos = SplitPositsIter.next();
                assertEquals("Wrong chromosome",chr,splitPos.getChromosome().getChromosomeNumber());
                Position origPos = OrigPositsIter.next();
                assertEquals("Original and split Positions are not the same for postion "+origPos.getSNPID()+" vs "+splitPos.getSNPID(), origPos.compareTo(splitPos), 0);
                assertEquals("Annotations different between original and split position "+origPos.getSNPID()+" vs "+splitPos.getSNPID(),origPos.toString(),splitPos.toString());
            }
            // TODO also check that the genos and depths are the same for the appropriate site range?
        }
        assertEquals("Not all sites were retained in the set of split genotype files",orig.numberOfSites(),totalNSites);
    }
    
}

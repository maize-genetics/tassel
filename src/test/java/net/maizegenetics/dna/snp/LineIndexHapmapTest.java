/*
 *  LineIndexHapmapTest
 * 
 *  Created on Oct 30, 2015
 */
package net.maizegenetics.dna.snp;

import net.maizegenetics.constants.TutorialConstants;
import org.junit.*;
import java.io.IOException;
import net.maizegenetics.constants.GeneralConstants;
import net.maizegenetics.dna.snp.io.BuilderFromHapMapLIX;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListBuilder;
import net.maizegenetics.taxa.Taxon;

/**
 *
 * @author Terry Casstevens
 */
public class LineIndexHapmapTest {

    public LineIndexHapmapTest() {
    }

    @Test
    public void testLineIndexHapmap() throws IOException {

        System.out.println("Testing Line Index Hapmap...");

        GenotypeTable expectedGenotype = ImportUtils.readFromHapmap(TutorialConstants.HAPMAP_FILENAME);

        GenotypeTable lixGenotype = BuilderFromHapMapLIX.build(TutorialConstants.HAPMAP_BGZIP_FILENAME, TutorialConstants.HAPMAP_LIX_FILENAME);

        AlignmentTestingUtils.alignmentsEqual(expectedGenotype, lixGenotype);

        GenotypeTable expectedFilter = ImportUtils.readFromHapmap(GeneralConstants.EXPECTED_RESULTS_DIR + "FilterAlignmentTest.hmp.txt");

        TaxaList origTaxa = lixGenotype.taxa();
        Taxon[] newIds = new Taxon[origTaxa.numberOfTaxa() / 3 + 1];
        int count = 0;
        for (int i = 0; i < origTaxa.numberOfTaxa(); i = i + 3) {
            newIds[count++] = origTaxa.get(i);
        }
        TaxaList taxa = new TaxaListBuilder().addAll(newIds).build();
        GenotypeTable filterTaxa = FilterGenotypeTable.getInstance(lixGenotype, taxa);

        double minFreq = 0.1;
        double maxFreq = 0.8;
        int minCount = 0;
        GenotypeTable filterTaxaSites = GenotypeTableUtils.removeSitesBasedOnFreqIgnoreMissing(filterTaxa, minFreq, maxFreq, minCount);

        AlignmentTestingUtils.alignmentsEqual(expectedFilter, filterTaxaSites);

    }
}

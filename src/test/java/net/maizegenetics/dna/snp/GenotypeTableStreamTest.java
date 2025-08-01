/*
 *  GenotypeTableStreamTest
 * 
 *  Created on Mar 13, 2015
 */
package net.maizegenetics.dna.snp;

import java.io.IOException;
import java.util.List;
import java.util.stream.Collectors;
import net.maizegenetics.constants.GeneralConstants;
import net.maizegenetics.constants.TutorialConstants;
import net.maizegenetics.dna.map.PositionListBuilder;
import net.maizegenetics.dna.snp.genotypecall.GenotypeCallTableBuilder;
import net.maizegenetics.dna.snp.genotypecall.GenotypeTest;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListBuilder;
import net.maizegenetics.taxa.Taxon;
import org.junit.Ignore;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author Terry Casstevens
 */
public class GenotypeTableStreamTest {

    public GenotypeTableStreamTest() {
    }

    @Test
    public void testGenotypeStream() throws IOException {
        System.out.println("Testing SuperByteMatrixSingle Genotype Stream...");
        GenotypeTable genotypes = ImportUtils.readFromHapmap(TutorialConstants.HAPMAP_FILENAME, null);
        testStream(genotypes);
        testTaxaStream(genotypes);
    }

    @Test
    public void testFilterGenotypeStream() throws IOException {
        System.out.println("Testing Filter Genotype Stream...");
        GenotypeTable genotypes = ImportUtils.readFromHapmap(TutorialConstants.HAPMAP_FILENAME, null);
        TaxaList origTaxa = genotypes.taxa();
        Taxon[] newIds = new Taxon[origTaxa.numberOfTaxa() / 3 + 1];
        int count = 0;
        for (int i = 0; i < origTaxa.numberOfTaxa(); i += 3) {
            newIds[count++] = origTaxa.get(i);
        }
        TaxaList taxa = new TaxaListBuilder().addAll(newIds).build();
        GenotypeTable filterTaxa = FilterGenotypeTable.getInstance(genotypes, taxa);
        testStream(filterTaxa);
        testTaxaStream(filterTaxa);

        double minFreq = 0.1;
        double maxFreq = 0.8;
        int minCount = 0;
        GenotypeTable filterTaxaSites = GenotypeTableUtils.removeSitesBasedOnFreqIgnoreMissing(filterTaxa, minFreq, maxFreq, minCount);
        testStream(filterTaxaSites);
        testTaxaStream(filterTaxaSites);
    }

    @Test
//    @Ignore("HDF5 exception")
    public void testHDF5GenotypeStream() throws IOException {
        System.out.println("Testing HDF5 Genotype Stream...");
        String filename = GeneralConstants.DATA_DIR + "CandidateTests/regHmp.hmp.t5.h5";
        GenotypeTable genotypes = ImportUtils.readGuessFormat(filename);
        testStream(genotypes);
        testTaxaStream(genotypes);
    }

    @Test
    public void testBigGenotypeStream() throws IOException {
        System.out.println("Testing SuperByteMatrixMultiple Genotype Stream...");
        int numTaxa = 30000;
        int numSites = 74000;
        GenotypeCallTableBuilder builder = GenotypeCallTableBuilder.getUnphasedNucleotideGenotypeBuilder(numTaxa, numSites);
        GenotypeTest.setRandomNucleotideData(builder);
        TaxaListBuilder taxaBuilder = new TaxaListBuilder();
        for (int t = 0; t < numTaxa; t++) {
            taxaBuilder.add(new Taxon("taxon" + t));
        }
        GenotypeTable genotypes = GenotypeTableBuilder.getInstance(builder.build(), PositionListBuilder.getInstance(numSites), taxaBuilder.build());
        testTaxaStream(genotypes);
    }

    private void testStream(GenotypeTable genotypes) {
        List<Byte> stream = genotypes.streamGenotype().collect(Collectors.toList());
        int count = 0;
        for (int t = 0; t < genotypes.numberOfTaxa(); t++) {
            for (int s = 0; s < genotypes.numberOfSites(); s++) {
                assertEquals("Genotypes don't match stream at taxon: " + t + "  site: " + s, genotypes.genotype(t, s), stream.get(count++).byteValue());
            }
        }
    }

    private void testTaxaStream(GenotypeTable genotypes) {
        for (int t = 0; t < genotypes.numberOfTaxa(); t++) {
            List<Byte> stream = genotypes.streamGenotype(t).collect(Collectors.toList());
            int count = 0;
            for (int s = 0; s < genotypes.numberOfSites(); s++) {
                assertEquals("Genotypes don't match stream at taxon: " + t + "  site: " + s, genotypes.genotype(t, s), stream.get(count++).byteValue());
            }
        }
    }
}

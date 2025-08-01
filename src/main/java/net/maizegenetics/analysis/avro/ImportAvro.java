/*
 *  ImportAvro
 * 
 *  Created on Oct 11, 2016
 */
package net.maizegenetics.analysis.avro;

import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import net.maizegenetics.dna.map.Chromosome;
import net.maizegenetics.dna.map.GeneralPosition;
import net.maizegenetics.dna.map.PositionList;
import net.maizegenetics.dna.map.PositionListBuilder;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableBuilder;
import net.maizegenetics.dna.snp.genotypecall.GOBIIAvroGenotypeCallTable;
import net.maizegenetics.dna.snp.genotypecall.GenotypeCallTable;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListBuilder;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.util.GeneralAnnotationStorage;
import org.apache.avro.generic.GenericArray;
import org.apache.avro.generic.GenericRecord;
import org.apache.avro.util.Utf8;

/**
 *
 * @author Terry Casstevens
 */
public class ImportAvro {

    private ImportAvro() {
        // utility
    }

    public static GenotypeTable genotypeTable(GenericRecord genotypeTable) {

        TaxaList taxa = taxa((GenericArray<GenericRecord>) genotypeTable.get(AvroConstants.GENOTYPE_TABLE_COMPONENTS.taxa.name()));
        PositionList positions = positions((GenericArray<GenericRecord>) genotypeTable.get(AvroConstants.GENOTYPE_TABLE_COMPONENTS.positions.name()));
        GenotypeCallTable genotypes = genotypeCallTable(taxa.numberOfTaxa(), positions.numberOfSites(), false, (GenericRecord) genotypeTable.get(AvroConstants.GENOTYPE_TABLE_COMPONENTS.genotype.name()));
        return GenotypeTableBuilder.getInstance(genotypes, positions, taxa, null, null, null, null, null);

    }

    private static GenotypeCallTable genotypeCallTable(int numTaxa, int numSites, boolean phased, GenericRecord genotypes) {
        return GOBIIAvroGenotypeCallTable.getInstance(numTaxa, numSites, phased, genotypes);
    }

    private static TaxaList taxa(GenericArray<GenericRecord> taxa) {

        TaxaListBuilder builder = new TaxaListBuilder();
        Iterator<GenericRecord> itr = taxa.iterator();
        while (itr.hasNext()) {
            GenericRecord current = itr.next();
            builder.add(new Taxon(current.get(AvroConstants.TAXON_INDICES.name.name()).toString(),
                    annotations((Map<Utf8, Utf8>) current.get(AvroConstants.TAXON_INDICES.annotations.name()))));
        }
        return builder.build();

    }

    private static PositionList positions(GenericArray<GenericRecord> positions) {

        PositionListBuilder builder = new PositionListBuilder();
        Iterator<GenericRecord> itr = positions.iterator();
        while (itr.hasNext()) {
            GenericRecord record = itr.next();
            GeneralPosition.Builder current = new GeneralPosition.Builder(chromosome((GenericRecord) record.get(AvroConstants.POSITION_INDICES.chromosome.name())),
                    (int) record.get(AvroConstants.POSITION_INDICES.position.name()))
                    .snpName(record.get(AvroConstants.POSITION_INDICES.snp_id.name()).toString());
            annotations((Map<Utf8, Utf8>) record.get(AvroConstants.POSITION_INDICES.annotations.name()), current);
            builder.add(current.build());
        }
        return builder.build();

    }

    private static final Map<String, Chromosome> CHROMOSOMES = new HashMap<>();

    private static Chromosome chromosome(GenericRecord chromosome) {
        String name = chromosome.get(AvroConstants.CHROMOSOME_INDICES.name.name()).toString();
        Chromosome result = CHROMOSOMES.get(name);
        if (result == null) {
            result = Chromosome.instance(name,
                    -1,
                    annotations((Map<Utf8, Utf8>) chromosome.get(AvroConstants.CHROMOSOME_INDICES.annotations.name())));
            CHROMOSOMES.put(name, result);
        }
        return result;
    }

    private static void annotations(Map<Utf8, Utf8> annotations, GeneralPosition.Builder positionBuilder) {

        if (annotations.isEmpty()) {
            return;
        }

        for (Map.Entry<Utf8, Utf8> current : annotations.entrySet()) {
            String key = current.getKey().toString();
            if (key.equals(AvroConstants.POSITION_STRAND)) {
                positionBuilder.strand(current.getValue().toString());
            } else {
                positionBuilder.addAnno(key, current.getValue());
            }
        }

    }

    private static GeneralAnnotationStorage annotations(Map<Utf8, Utf8> annotations) {

        if (annotations.isEmpty()) {
            return null;
        }

        GeneralAnnotationStorage.Builder builder = GeneralAnnotationStorage.getBuilder();
        for (Map.Entry<Utf8, Utf8> current : annotations.entrySet()) {
            builder.addAnnotation(current.getKey().toString(), current.getValue().toString());
        }
        return builder.build();

    }

}

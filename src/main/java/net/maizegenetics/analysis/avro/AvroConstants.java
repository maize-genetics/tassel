/*
 *  AvroConstants
 * 
 *  Created on Oct 10, 2016
 */
package net.maizegenetics.analysis.avro;

import net.maizegenetics.dna.snp.genotypecall.GOBIIAvroGenotypeCallTable;
import org.apache.avro.Schema;
import org.apache.avro.SchemaBuilder;

/**
 *
 * @author Terry Casstevens
 */
public class AvroConstants {

    public static final int GENOTYPE_BLOCK_SIZE = GOBIIAvroGenotypeCallTable.GENOTYPE_BLOCK_SIZE;

    public static final Schema BYTE_SCHEMA = SchemaBuilder
            .builder("net.maizegenetics")
            .bytesType();

    public static final Schema BYTE_BLOCK_SCHEMA = SchemaBuilder
            .builder("net.maizegenetics")
            .array().items(BYTE_SCHEMA);

    public static final Schema ANNOTATION_SCHEMA = SchemaBuilder
            .builder("net.maizegenetics")
            .map().values().stringType();

    public static enum TAXON_INDICES {
        name, annotations
    };

    public static final Schema TAXON_SCHEMA = SchemaBuilder
            .builder("net.maizegenetics")
            .record("taxon")
            .fields()
            .name(TAXON_INDICES.name.name()).type().stringType().noDefault()
            .name(TAXON_INDICES.annotations.name()).type(ANNOTATION_SCHEMA).noDefault()
            .endRecord();

    public static final Schema TAXA_SCHEMA = SchemaBuilder
            .builder("net.maizegenetics")
            .array().items(TAXON_SCHEMA);

    public static enum CHROMOSOME_INDICES {
        name, annotations
    };

    public static final Schema CHROMOSOME_SCHEMA = SchemaBuilder
            .builder("net.maizegenetics")
            .record("chromosome")
            .fields()
            .name(CHROMOSOME_INDICES.name.name()).type().stringType().noDefault()
            .name(CHROMOSOME_INDICES.annotations.name()).type(ANNOTATION_SCHEMA).noDefault()
            .endRecord();

    public static enum POSITION_INDICES {
        chromosome, position, snp_id, annotations
    };

    public static final Schema POSITION_SCHEMA = SchemaBuilder
            .builder("net.maizegenetics")
            .record("position")
            .fields()
            .name(POSITION_INDICES.chromosome.name()).type(CHROMOSOME_SCHEMA).noDefault()
            .name(POSITION_INDICES.position.name()).type().intType().noDefault()
            .name(POSITION_INDICES.snp_id.name()).type().stringType().noDefault()
            .name(POSITION_INDICES.annotations.name()).type(ANNOTATION_SCHEMA).noDefault()
            .endRecord();

    public static final Schema POSITIONS_SCHEMA = SchemaBuilder
            .builder("net.maizegenetics")
            .array().items(POSITION_SCHEMA);

    public static final String POSITION_STRAND = "strand";
    public static final String POSITION_CM = "cm";
    public static final String POSITION_IS_NUCLEOTIDE = "isNucleotide";
    public static final String POSITION_IS_INDEL = "isIndel";
    public static final String POSITION_MAF = "maf";
    public static final String POSITION_SITE_COVERAGE = "siteCoverage";
    public static final String POSITION_ALLELE_VALUE = "alleleValue";

    public static enum GENOTYPE_TABLE_COMPONENTS {
        taxa, positions, genotype
    };

    public static long getCacheKey(int taxon, int site) {
        return GOBIIAvroGenotypeCallTable.getCacheKey(taxon, site);
    }

    public static String getKey(int taxon, int site) {
        return GOBIIAvroGenotypeCallTable.getKey(taxon, site);
    }

    public static int[] getTaxonSiteFromKey(String keyStr) {
        long key = Long.valueOf(keyStr.substring(1));
        int[] result = new int[2];
        result[0] = (int) (key >>> 32) * GENOTYPE_BLOCK_SIZE;
        result[1] = (int) (key & 0xFFFFFFFF) * GENOTYPE_BLOCK_SIZE;
        return result;
    }

    private AvroConstants() {
        // utility
    }

}

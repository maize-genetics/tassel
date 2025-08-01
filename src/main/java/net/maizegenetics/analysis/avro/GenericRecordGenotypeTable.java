/*
 *  GenericRecordGenotypeTable
 * 
 *  Created on Oct 10, 2016
 */
package net.maizegenetics.analysis.avro;

import net.maizegenetics.dna.snp.GenotypeTable;
import org.apache.avro.Schema;
import org.apache.avro.generic.GenericRecord;

/**
 *
 * @author Terry Casstevens
 */
public class GenericRecordGenotypeTable implements GenericRecord {

    private final Schema mySchema;
    private final GenotypeTable myTable;

    public GenericRecordGenotypeTable(Schema tasselSchema, GenotypeTable table) {
        mySchema = tasselSchema;
        myTable = table;
    }

    @Override
    public void put(String arg0, Object arg1) {
        throw new UnsupportedOperationException("Not Mutable.");
    }

    @Override
    public Object get(String key) {
        switch (key) {
            case "taxa":
                return new GenericArrayTaxa(myTable.taxa());
            case "positions":
                return new GenericArrayPositions(myTable.positions());
            case "genotype":
                return new GenericRecordGenotype(mySchema.getField("genotype").schema(), myTable);
            default:
                throw new IllegalStateException("GenericRecordGenotypeTable: get: Unknown key: " + key);
        }
    }

    @Override
    public void put(int arg0, Object arg1) {
        throw new UnsupportedOperationException("Not Mutable.");
    }

    @Override
    public Object get(int i) {
        return get(AvroConstants.GENOTYPE_TABLE_COMPONENTS.values()[i].name());
    }

    @Override
    public Schema getSchema() {
        return mySchema;
    }

}

/*
 *  GenericRecordTaxon
 * 
 *  Created on Oct 11, 2016
 */
package net.maizegenetics.analysis.avro;

import net.maizegenetics.taxa.Taxon;
import org.apache.avro.Schema;
import org.apache.avro.generic.GenericRecord;

/**
 *
 * @author Terry Casstevens
 */
public class GenericRecordTaxon implements GenericRecord {

    private final Taxon myTaxon;

    public GenericRecordTaxon(Taxon taxon) {
        myTaxon = taxon;
    }

    @Override
    public void put(String arg0, Object arg1) {
        throw new UnsupportedOperationException("Not Mutable.");
    }

    @Override
    public Object get(String key) {

        if (key.equals(AvroConstants.TAXON_INDICES.name.name())) {
            return myTaxon.getName();
        } else if (key.equals(AvroConstants.TAXON_INDICES.annotations.name())) {
            return new GenericMapAnnotations(myTaxon.getAnnotation());
        } else {
            throw new IllegalArgumentException("GenericRecordTaxon: get: unknown key: " + key);
        }

    }

    @Override
    public void put(int arg0, Object arg1) {
        throw new UnsupportedOperationException("Not Mutable.");
    }

    @Override
    public Object get(int i) {
        return get(AvroConstants.TAXON_INDICES.values()[i].name());
    }

    @Override
    public Schema getSchema() {
        return AvroConstants.TAXON_SCHEMA;
    }

}

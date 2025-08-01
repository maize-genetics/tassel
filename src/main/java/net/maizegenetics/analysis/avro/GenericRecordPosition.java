/*
 *  GenericRecordPosition
 * 
 *  Created on Oct 31, 2016
 */

package net.maizegenetics.analysis.avro;

import net.maizegenetics.dna.map.Position;
import org.apache.avro.Schema;
import org.apache.avro.generic.GenericRecord;

/**
 *
 * @author Terry Casstevens
 */
public class GenericRecordPosition implements GenericRecord {

    private final Position myPosition;

    public GenericRecordPosition(Position position) {
        myPosition = position;
    }

    @Override
    public void put(String arg0, Object arg1) {
        throw new UnsupportedOperationException("Not Mutable.");
    }

    @Override
    public Object get(String key) {

        if (key.equals(AvroConstants.POSITION_INDICES.chromosome.name())) {
            return GenericRecordChromosome.getInstance(myPosition.getChromosome());
        } else if (key.equals(AvroConstants.POSITION_INDICES.position.name())) {
            return myPosition.getPosition();
        } else if (key.equals(AvroConstants.POSITION_INDICES.snp_id.name())) {
            return myPosition.getSNPID();
        } else if (key.equals(AvroConstants.POSITION_INDICES.annotations.name())) {
            return new GenericMapAnnotations(myPosition.getAnnotation());
        } else {
            throw new IllegalArgumentException("GenericRecordPosition: get: unknown key: " + key);
        }

    }

    @Override
    public void put(int arg0, Object arg1) {
        throw new UnsupportedOperationException("Not Mutable.");
    }

    @Override
    public Object get(int i) {
        return get(AvroConstants.POSITION_INDICES.values()[i].name());
    }

    @Override
    public Schema getSchema() {
        return AvroConstants.POSITION_SCHEMA;
    }

}


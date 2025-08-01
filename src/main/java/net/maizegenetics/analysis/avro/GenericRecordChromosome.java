/*
 *  GenericRecordChromosome
 * 
 *  Created on Oct 31, 2016
 */
package net.maizegenetics.analysis.avro;

import java.util.HashMap;
import java.util.Map;
import net.maizegenetics.dna.map.Chromosome;
import org.apache.avro.Schema;
import org.apache.avro.generic.GenericRecord;

/**
 *
 * @author Terry Casstevens
 */
public class GenericRecordChromosome implements GenericRecord {
    
    private static Map<Chromosome, GenericRecordChromosome> INSTANCES = new HashMap<>();

    private final Chromosome myChromosome;

    private GenericRecordChromosome(Chromosome chromosome) {
        myChromosome = chromosome;
    }
    
    public static GenericRecordChromosome getInstance(Chromosome chromosome) {
        GenericRecordChromosome result = INSTANCES.get(chromosome);
        if (result == null) {
            result = new GenericRecordChromosome(chromosome);
            INSTANCES.put(chromosome, result);
        }
        return result;
    }

    @Override
    public void put(String arg0, Object arg1) {
        throw new UnsupportedOperationException("Not Mutable.");
    }

    @Override
    public Object get(String key) {

        if (key.equals(AvroConstants.CHROMOSOME_INDICES.name.name())) {
            return myChromosome.getName();
        } else if (key.equals(AvroConstants.CHROMOSOME_INDICES.annotations.name())) {
            return new GenericMapAnnotations(myChromosome.getAnnotation());
        } else {
            throw new IllegalArgumentException("GenericRecordChromosome: get: unknown key: " + key);
        }

    }

    @Override
    public void put(int arg0, Object arg1) {
        throw new UnsupportedOperationException("Not Mutable.");
    }

    @Override
    public Object get(int i) {
        return get(AvroConstants.CHROMOSOME_INDICES.values()[i].name());
    }

    @Override
    public Schema getSchema() {
        return AvroConstants.CHROMOSOME_SCHEMA;
    }

}

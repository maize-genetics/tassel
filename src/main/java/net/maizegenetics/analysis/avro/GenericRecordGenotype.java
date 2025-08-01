/*
 *  GenericRecordGenotype
 * 
 *  Created on Sep 10, 2016
 */
package net.maizegenetics.analysis.avro;

import java.io.File;
import java.nio.ByteBuffer;
import net.maizegenetics.dna.snp.GenotypeTable;
import org.apache.avro.Schema;
import org.apache.avro.file.DataFileReader;
import org.apache.avro.generic.GenericData;
import org.apache.avro.generic.GenericDatumReader;
import org.apache.avro.generic.GenericRecord;
import org.apache.avro.io.DatumReader;

/**
 *
 * @author Terry Casstevens
 */
public class GenericRecordGenotype implements GenericRecord {

    private final Schema myGenotypeSchema;
    private final GenotypeTable myGenotype;
    private final int myNumTaxa;
    private final int myNumSites;
    private final int myNumTaxaBlocks;

    public GenericRecordGenotype(Schema genotypeSchema, GenotypeTable genotype) {
        myGenotypeSchema = genotypeSchema;
        myGenotype = genotype;
        myNumTaxa = myGenotype.numberOfTaxa();
        myNumSites = myGenotype.numberOfSites();
        int temp = myNumTaxa / AvroConstants.GENOTYPE_BLOCK_SIZE;
        if (myNumTaxa % AvroConstants.GENOTYPE_BLOCK_SIZE == 0) {
            myNumTaxaBlocks = temp;
        } else {
            myNumTaxaBlocks = temp + 1;
        }
    }

    @Override
    public void put(String arg0, Object arg1) {
        throw new UnsupportedOperationException("Not Mutable.");
    }

    @Override
    public Object get(String key) {
        System.out.println("GenericRecordGenotype: get: key: " + key);
        int[] taxonSite = AvroConstants.getTaxonSiteFromKey(key);
        int tSize = Math.min(AvroConstants.GENOTYPE_BLOCK_SIZE, myNumTaxa - taxonSite[0]);
        int sSize = Math.min(AvroConstants.GENOTYPE_BLOCK_SIZE, myNumSites - taxonSite[1]);
        System.out.println("tSize: " + tSize + "  sSize: " + sSize + "   taxon: " + taxonSite[0] + "  site: " + taxonSite[1]);
        GenericData.Array<ByteBuffer> result = new GenericData.Array<>(sSize, AvroConstants.BYTE_BLOCK_SCHEMA);
        for (int s = 0; s < sSize; s++) {
            byte[] genotype = myGenotype.genotypeAllTaxa(taxonSite[1] + s);
            ByteBuffer temp = ByteBuffer.allocateDirect(tSize);
            temp.rewind();
            for (int t = 0; t < tSize; t++) {
                temp.put(genotype[taxonSite[0] + t]);
            }
            temp.rewind();
            result.add(temp);
        }
        return result;
    }

    @Override
    public void put(int arg0, Object arg1) {
        throw new UnsupportedOperationException("Not Mutable.");
    }

    @Override
    public Object get(int i) {
        System.out.println("GenericRecordGenotype: get: " + i);
        int siteBlock = i / myNumTaxaBlocks;
        int taxaBlock = i % myNumTaxaBlocks;
        long key = ((long) taxaBlock << 32) + siteBlock;
        return get("B" + key);
    }

    @Override
    public Schema getSchema() {
        return myGenotypeSchema;
    }

    public static void main(String[] args) {
        try {
            DatumReader<GenericRecord> datumReader = new GenericDatumReader<>();
            DataFileReader<GenericRecord> reader = new DataFileReader<>(new File("tassel.avro"), datumReader);
            GenericRecord temp = reader.next();
            GenericData.Array<ByteBuffer> temp1 = (GenericData.Array<ByteBuffer>) temp.get("B0");
            System.out.println("temp1: " + temp1.getClass().getName());
            System.out.println("temp1 get 0: " + temp1.get(0).getClass().getName());
            System.out.println("temp1 size: " + temp1.size());
            System.out.println(temp1.get(1).array().length);
            ByteBuffer temp3 = temp1.get(2);
            while (temp3.hasRemaining()) {
                System.out.print(temp3.get() + " ");
            }
            System.out.println("");
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

}

/*
 *  ExportAvro
 * 
 *  Created on Nov 3, 2016
 */
package net.maizegenetics.analysis.avro;

import java.io.File;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.ImportUtils;
import net.maizegenetics.util.LoggingUtils;
import net.maizegenetics.util.Utils;
import org.apache.avro.Schema;
import org.apache.avro.SchemaBuilder;
import org.apache.avro.file.CodecFactory;
import org.apache.avro.file.DataFileWriter;
import org.apache.avro.generic.GenericDatumWriter;
import org.apache.avro.generic.GenericRecord;
import org.apache.avro.io.DatumWriter;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 *
 * @author Terry Casstevens
 */
public class ExportAvro {

    private static final Logger myLogger = LogManager.getLogger(ExportAvro.class);

    private ExportAvro() {
        // utility
    }

    public static String write(GenotypeTable genotype, String filename) {

        filename = Utils.addSuffixIfNeeded(filename, ".avro");

        try {

            int numTaxa = genotype.numberOfTaxa();
            int numSites = genotype.numberOfSites();

            SchemaBuilder.FieldAssembler<Schema> genotypeSchemaBuilder = SchemaBuilder
                    .builder("net.maizegenetics")
                    .record("genotype")
                    .fields();
            for (int s = 0; s < numSites; s += AvroConstants.GENOTYPE_BLOCK_SIZE) {
                for (int t = 0; t < numTaxa; t += AvroConstants.GENOTYPE_BLOCK_SIZE) {
                    genotypeSchemaBuilder = genotypeSchemaBuilder
                            .name(AvroConstants.getKey(t, s))
                            .type(AvroConstants.BYTE_BLOCK_SCHEMA)
                            .noDefault();
                }
            }
            Schema genotypeSchema = genotypeSchemaBuilder.endRecord();

            Schema tasselSchema = SchemaBuilder
                    .builder("net.maizegenetics")
                    .record("tassel")
                    .fields()
                    .name("taxa").type(AvroConstants.TAXA_SCHEMA).noDefault()
                    .name("positions").type(AvroConstants.POSITIONS_SCHEMA).noDefault()
                    .name("genotype").type(genotypeSchema).noDefault()
                    .endRecord();

            DatumWriter<GenericRecord> datumWriter = new GenericDatumWriter<>(tasselSchema);
            try (DataFileWriter<GenericRecord> dataFileWriter = new DataFileWriter<>(datumWriter)) {
                dataFileWriter.setCodec(CodecFactory.snappyCodec());
                dataFileWriter.create(tasselSchema, new File(filename));
                dataFileWriter.append(new GenericRecordGenotypeTable(tasselSchema, genotype));
            }

        } catch (Exception e) {
            myLogger.debug(e.getMessage(), e);
            throw new IllegalStateException("ExportAvro: write: problem writing file: " + filename + ". " + e.getMessage());
        }
        
        return filename;

    }
    
    public static void main(String[] args) {
        LoggingUtils.setupDebugLogging();
        GenotypeTable genotype = ImportUtils.read("mdp_genotype.hmp.txt");
        write(genotype, "test");
    }

}

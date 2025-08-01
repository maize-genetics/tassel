/*
 *  GenerateRCode
 */
package net.maizegenetics.plugindef;

import com.google.common.base.CaseFormat;
import net.maizegenetics.analysis.association.FastMultithreadedAssociationPlugin;
import net.maizegenetics.analysis.association.FixedEffectLMPlugin;
import net.maizegenetics.analysis.association.GenomicSelectionPlugin;
import net.maizegenetics.analysis.association.MLMPlugin;
import net.maizegenetics.analysis.distance.KinshipPlugin;
import net.maizegenetics.analysis.filter.FilterSiteBuilderPlugin;
import net.maizegenetics.analysis.filter.FilterTaxaBuilderPlugin;
import net.maizegenetics.analysis.popgen.LinkageDisequilibrium;
import net.maizegenetics.dna.map.*;
import net.maizegenetics.dna.snp.*;
import net.maizegenetics.dna.snp.genotypecall.GenotypeCallTableBuilder;
import net.maizegenetics.dna.snp.io.FlapjackUtils;
import net.maizegenetics.phenotype.*;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListBuilder;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.taxa.distance.DistanceMatrix;
import net.maizegenetics.taxa.distance.DistanceMatrixBuilder;
import net.maizegenetics.util.BitSet;
import net.maizegenetics.util.OpenBitSet;
import net.maizegenetics.util.TableReport;
import net.maizegenetics.util.Utils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.awt.*;
import java.lang.reflect.Constructor;
import java.lang.reflect.Field;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;
import java.util.List;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import static net.maizegenetics.dna.map.Position.STRAND_MINUS;
import static net.maizegenetics.dna.map.Position.STRAND_PLUS;

/**
 * @author Terry Casstevens
 * @author Ed Buckler
 * @author Brandon Monier
 */
public class GenerateRCode {

    private static final Logger myLogger = LogManager.getLogger(GenerateRCode.class);

    private GenerateRCode() {

    }

    public static void main(String[] args) {
        printHeader();
        generate(FilterSiteBuilderPlugin.class, "genotypeTable", "genotypeTable");
        generate(FilterTaxaBuilderPlugin.class, "genotypeTable", "genotypeTable");
        generate(KinshipPlugin.class, "genotypeTable", "distanceMatrix");
        generate(FixedEffectLMPlugin.class, "phenotypeGenotypeTable", "tableReport");

    }

    private static void printHeader() {

        System.out.println("#!/usr/bin/env Rscript");

        System.out.println("\n#--------------------------------------------------------------------");
        System.out.println("# Script Name:   TasselPluginWrappers.R");
        System.out.println("# Description:   Generated R interface to TASSEL 5");
        System.out.println("# Author:        Brandon Monier, Ed Buckler, Terry Casstevens");
        System.out.print("# Created:       ");
        System.out.println(new Date());
        System.out.println("#--------------------------------------------------------------------");

        System.out.println("# Preamble\n");

        System.out.println("\n## Load packages");
        System.out.println("if (!requireNamespace(\"BiocManager\")) {");
        System.out.println("    install.packages(\"BiocManager\")");
        System.out.println("}");

        System.out.println("\npackages <- c(");
        System.out.println("\"rJava\"");
        System.out.println(")");
        System.out.println("BiocManager::install(packages)");
        System.out.println("library(rJava)");

        System.out.println("\n## Init JVM");
        System.out.println("rJava::.jinit()");

        System.out.println("\n## Add TASSEL 5 class path");
        System.out.println("rJava::.jaddClassPath(\"/tassel-5-standalone/lib\")");
        System.out.println("rJava::.jaddClassPath(\"/tassel-5-standalone/sTASSEL.jar\")\n");

        System.out.println("source(\"R/AllClasses.R\")\n");

    }

    public static void generate(Class currentMatch, String inputObject, String outputObject) {
        try {
            Constructor constructor = currentMatch.getConstructor(Frame.class);
            generate((AbstractPlugin) constructor.newInstance((Frame) null), inputObject, outputObject);
        } catch (Exception ex) {
            try {
                Constructor constructor = currentMatch.getConstructor(Frame.class, boolean.class);
                generate((AbstractPlugin) constructor.newInstance(null, false), inputObject, outputObject);
            } catch (NoSuchMethodException nsme) {
                myLogger.warn("Self-describing Plugins should implement this constructor: " + currentMatch.getClass().getName());
                myLogger.warn("public Plugin(Frame parentFrame, boolean isInteractive) {");
                myLogger.warn("   super(parentFrame, isInteractive);");
                myLogger.warn("}");
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
    }

    private static void generate(AbstractPlugin plugin, String inputObject, String outputObject) {

        String clazz = Utils.getBasename(plugin.getClass().getName());
        //StringBuilder sb = new StringBuilder("rTASSEL::" + stringToCamelCase(clazz));
        StringBuilder sb = new StringBuilder(stringToCamelCase(clazz));
        sb.append(" <- function(");
        sb.append(inputObject + ",\n");
        for (Field field : plugin.getParameterFields()) {
            PluginParameter<?> current = null;
            try {
                current = (PluginParameter) field.get(plugin);
            } catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }
            String methodName = current.cmdLineName();
            if ((current.defaultValue() instanceof Number) || (current.defaultValue() instanceof Boolean) || (current.defaultValue() instanceof String)) {
                sb.append("            " + methodName + "=" + current.defaultValue() + ",\n");
            } else if ((current.defaultValue() instanceof Enum)) {
                sb.append("            " + methodName + "=\"" + current.defaultValue() + "\",\n");
            } else if ((current.defaultValue() instanceof Boolean)) {
                sb.append("            " + methodName + "=\"" + current.defaultValue() + "\",\n");
            }

        }
        sb.deleteCharAt(sb.lastIndexOf(",")); //remove last comma

        sb.append(") {\n");

        if (inputObject.equals("genotypeTable")) {
            sb.append("    if(is(genotypeTable, \"GenotypeTable\") == TRUE) {\n");
            sb.append("        genotypeTable <- genotypeTable@jtsGenotypeTable\n");
            sb.append("    }\n");
        }

        sb.append("    plugin <- new(J(\"" + plugin.getClass().getCanonicalName() + "\"), .jnull(), FALSE)\n");
        //sb.append("    plugin <- rJava::.jnew(\"" + plugin.getClass().getCanonicalName() + "\")\n");
        for (Field field : plugin.getParameterFields()) {
            PluginParameter<?> current = null;
            try {
                current = (PluginParameter) field.get(plugin);
            } catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }
            String methodName = current.cmdLineName();
            if ((current.defaultValue() instanceof Number) || (current.defaultValue() instanceof Boolean) || (current.defaultValue() instanceof String)) {
                sb.append("    plugin$setParameter(\"" + methodName + "\",toString(" + methodName + "))\n");
            } else if ((current.defaultValue() instanceof Enum)) {
                sb.append("    plugin$setParameter(\"" + methodName + "\",toString(" + methodName + "))\n");
            } else if ((current.defaultValue() instanceof Boolean)) {
                sb.append("    plugin$setParameter(\"" + methodName + "\",toString(" + methodName + "))\n");
            }
        }

        if (outputObject.equals("genotypeTable")) {
            sb.append("    filteredGT <- plugin$runPlugin(" + inputObject + ")\n");
            sb.append("    new(\n");
            sb.append("        Class = \"GenotypeTable\",\n");
            sb.append("        name = paste0(\"Filtered:\"),\n");
            sb.append("        jtsGenotypeTable = filteredGT\n");
            sb.append("    )\n");
        } else {
            sb.append("    plugin$runPlugin(" + inputObject + ")\n");
        }

        sb.append("}\n");
        System.out.println(sb.toString());

    }

    private static final int DEFAULT_DESCRIPTION_LINE_LENGTH = 50;

    private static String createDescription(String description) {
        int count = 0;
        StringBuilder builder = new StringBuilder();
        builder.append("     * ");
        for (int i = 0, n = description.length(); i < n; i++) {
            count++;
            if (description.charAt(i) == '\n') {
                builder.append("\n");
                builder.append("     * ");
                count = 0;
            } else if ((count > DEFAULT_DESCRIPTION_LINE_LENGTH) && (description.charAt(i) == ' ')) {
                builder.append("\n");
                builder.append("     * ");
                count = 0;
            } else {
                builder.append(description.charAt(i));
            }
        }
        return builder.toString();
    }

    private static String stringToCamelCase(String str) {
        StringBuilder builder = new StringBuilder();
        builder.append(Character.toLowerCase(str.charAt(0)));
        boolean makeUpper = false;
        for (int i = 1; i < str.length(); i++) {
            char current = str.charAt(i);
            if (current == ' ') {
                makeUpper = true;
            } else if (makeUpper) {
                builder.append(Character.toUpperCase(current));
                makeUpper = false;
            } else {
                builder.append(current);
            }
        }
        return builder.toString();
    }

    private static String removeMyFromString(String str) {
        String lower = str.toLowerCase();
        if (lower.startsWith("my")) {
            str = str.substring(2);
        }
        return CaseFormat.UPPER_CAMEL.to(CaseFormat.LOWER_CAMEL, str);
    }

    /**
     * This converts a genotype table to a double dimension dosage byte array.
     *
     * @param genotype genotype table
     *
     * @return byte[][] of dosage values.  First index is taxa and second index is site.
     * Byte.MIN_VALUE is used for NA
     */
    public static byte[][] genotypeTableToDosageByteArray(GenotypeTable genotype) {

        byte[][] result = new byte[genotype.numberOfTaxa()][genotype.numberOfSites()];

        PositionVectors posVectors = genotypeTableToPositionListOfArrays(genotype.positions());
        for (int site = 0; site < genotype.numberOfSites(); site++) {
            byte[] siteGenotypes = genotype.genotypeAllTaxa(site);

            byte refAllele = NucleotideAlignmentConstants.getNucleotideAlleleByte(posVectors.refAllele[site]);

            // value assigned to site / taxon is the number of alleles
            // that do not match the ref allele / major allele.
            for (int taxon = 0; taxon < genotype.numberOfTaxa(); taxon++) {
                byte value = 0;
                byte[] alleles = GenotypeTableUtils.getDiploidValues(siteGenotypes[taxon]);
                if (alleles[0] == GenotypeTable.UNKNOWN_ALLELE || alleles[1] == GenotypeTable.UNKNOWN_ALLELE) {
                    value = Byte.MIN_VALUE;
                } else {
                    if (alleles[0] != refAllele) value++;
                    if (alleles[1] != refAllele) value++;
                }
                result[taxon][site] = value;
            }
        }
        return result;
    }

    /**
     * This converts the taxa of the specified genotype table to
     * an array of taxa name strings.
     *
     * @param genotype genotype table
     *
     * @return string array of taxa names
     */
    public static String[] genotypeTableToSampleNameArray(GenotypeTable genotype) {
        return genotypeTableToSampleNameArray(genotype.taxa());
    }

    /**
     * This converts the taxa list to an array of taxa name strings.
     *
     * @param taxa taxa list
     *
     * @return string array of taxa names
     */
    public static String[] genotypeTableToSampleNameArray(TaxaList taxa) {
        return taxa.stream()
                .map(Taxon::getName)
                .toArray(String[]::new);
    }

    /**
     * This converts the given position list to a @{@link PositionVectors}
     *
     * @param positions
     *
     * @return @{@link PositionVectors}
     */
    public static PositionVectors genotypeTableToPositionListOfArrays(PositionList positions) {

        String[] chromosomes = new String[positions.numberOfSites()];
        int[] startPos = new int[positions.numberOfSites()];
        int[] strand = new int[positions.numberOfSites()];
        String[] refAllele = new String[positions.numberOfSites()];
        String[] altAllele = new String[positions.numberOfSites()];

        for (int site = 0; site < positions.numberOfSites(); site++) {
            Position p = positions.get(site);
            chromosomes[site] = p.getChromosome().getName();
            startPos[site] = p.getPosition();
            strand[site] = (p.getStrand() == STRAND_PLUS) ? (1) : ((p.getStrand() == STRAND_MINUS) ? (-1) : (Integer.MIN_VALUE));
            String[] variants = p.getKnownVariants();
            refAllele[site] = (variants.length > 0) ? variants[0] : "";
            altAllele[site] = (variants.length > 1) ? variants[1] : "";
        }

        return new PositionVectors(chromosomes, startPos, strand, refAllele, altAllele);
    }

    /**
     * This converts the position list of the given genotype table to
     * a @{@link PositionVectors}
     *
     * @param genotype genotype table
     *
     * @return @{@link PositionVectors}
     */
    public static PositionVectors genotypeTableToPositionListOfArrays(GenotypeTable genotype) {
        return genotypeTableToPositionListOfArrays(genotype.positions());
    }

    public static class PositionVectors {
        public String[] chromosomes;
        public int[] startPos;
        public int[] strand;
        public String[] refAllele;
        public String[] altAllele;

        public PositionVectors(String[] chromosomes, int[] startPos, int[] strand, String[] refAllele, String[] altAllele) {
            this.chromosomes = chromosomes;
            this.startPos = startPos;
            this.strand = strand;
            this.refAllele = refAllele;
            this.altAllele = altAllele;
        }
    }

    /**
     * Temporary place for this experimental method.
     *
     * @param tableReport
     *
     * @return int[] in column order with NA set to R approach
     */
    public static TableReportVectors tableReportToVectors(TableReport tableReport) {
        if (tableReport.getRowCount() > Integer.MAX_VALUE)
            throw new IndexOutOfBoundsException("R cannot handle more than " + Integer.MAX_VALUE + " rows");
        int rows = (int) tableReport.getRowCount();
        String[] columnNames = Stream.of(tableReport.getTableColumnNames()).map(Object::toString).toArray(String[]::new);
        List dataVector = new ArrayList();
        String[] columnType = Stream.of(tableReport.getRow(0)).map(t -> t.getClass().toString()).toArray(String[]::new);
        for (int column = 0; column < tableReport.getColumnCount(); column++) {
            Object o = tableReport.getValueAt(0, column);
            if (o instanceof Float || o instanceof Double) {
                double[] result = new double[rows];
                for (int row = 0; row < tableReport.getRowCount(); row++) {
                    Object value = tableReport.getValueAt(row, column);
                    if (value == null) {
                        result[row] = Double.NaN;
                    } else {
                        result[row] = ((Number) value).doubleValue();
                    }
                }
                dataVector.add(result);
            } else if (o instanceof Byte || o instanceof Short || o instanceof Integer || o instanceof Long) {
                int[] result = new int[rows];
                for (int row = 0; row < tableReport.getRowCount(); row++) {
                    try {
                        result[row] = ((Number) tableReport.getValueAt(row, column)).intValue();
                    } catch (ClassCastException cce) {
                        result[row] = Integer.MIN_VALUE;
                    }
                }
                dataVector.add(result);
            } else {
                String[] result = new String[rows];
                for (int row = 0; row < tableReport.getRowCount(); row++) {
                    Object value = tableReport.getValueAt(row, column);
                    if (value == null) {
                        result[row] = "";
                    } else {
                        result[row] = value.toString();
                    }
                }
                dataVector.add(result);
            }
        }
        String[] annotation = new String[tableReport.getColumnCount()];
        return new TableReportVectors(columnNames, columnType, annotation, dataVector);
    }

    public static class TableReportVectors {
        public String[] columnNames;
        public String[] columnType;
        public String[] annotation;
        public List dataVector = new ArrayList();

        public TableReportVectors(String[] columnNames, String[] columnType, String[] annotation, List dataVector) {
            this.columnNames = columnNames;
            this.columnType = columnType;
            this.annotation = annotation;
            this.dataVector = dataVector;
        }
    }

    /**
     * This converts a GIGWA R Dataframe to a GenotypeTable
     */
    public static GenotypeTable createGenotypeFromRDataFrameElements(
            String[] taxa,
            String[] chromosomes,
            int[] snpPos,
            String[] snpId,
            String[] alleles,
            int[][] markerMatrix
    ) {

        if (chromosomes.length != snpPos.length || chromosomes.length != snpId.length || chromosomes.length != alleles.length)
            throw new IllegalArgumentException("createGenotypeFromRDataFrameElements: length of chrom, snpPos, snpId, and alleles must be equal");

        TaxaListBuilder taxaListBuilder = new TaxaListBuilder();
        for (String taxon : taxa) {
            taxaListBuilder.add(new Taxon(taxon));
        }
        TaxaList taxaList = taxaListBuilder.build();

        int numberOfTaxa = taxaList.numberOfTaxa();
        int numberOfSites = snpPos.length;

        if (markerMatrix.length != numberOfSites)
            throw new IllegalArgumentException("createGenotypeFromRDataFrameElements: number of rows in markerMatrix must equal number of positions");

        if (markerMatrix[0].length != numberOfTaxa)
            throw new IllegalArgumentException("createGenotypeFromRDataFrameElements: number of columns in markerMatrix must equal number of taxa");

        String delimiter = "/";
        boolean phased = false;
        if (alleles[0].contains("/")) {
            delimiter = "/";
            phased = false;
        } else if (alleles[0].contains("|")) {
            delimiter = "|";
            phased = true;
        } else {
            throw new IllegalArgumentException("createGenotypeFromRDataFrameElements: alleles must be in the form A/C or A|C");
        }

        GenotypeCallTableBuilder genotypeCallTableBuilder = null;
        if (phased) {
            genotypeCallTableBuilder = GenotypeCallTableBuilder.getUnphasedNucleotideGenotypeBuilder(numberOfTaxa, numberOfSites);
            genotypeCallTableBuilder.isPhased(true);
        } else {
            genotypeCallTableBuilder = GenotypeCallTableBuilder.getUnphasedNucleotideGenotypeBuilder(numberOfTaxa, numberOfSites);
        }

        PositionListBuilder positionListBuilder = new PositionListBuilder();

        for (int sIdx = 0; sIdx < snpPos.length; sIdx++) {

            String[] variants = alleles[sIdx].split(delimiter);
            if (variants.length != 2)
                throw new IllegalArgumentException("createGenotypeFromRDataFrameElements: there must be exactly two alleles per position: " + alleles[sIdx]);

            Position position = new GeneralPosition.Builder(Chromosome.instance(chromosomes[sIdx]), snpPos[sIdx])
                    .knownVariants(variants)
                    .snpName(snpId[sIdx])
                    .build();
            positionListBuilder.add(position);

            for (int tIdx = 0; tIdx < numberOfTaxa; tIdx++) {

                // 0 = homozygous reference, 1 = heterozygous, 2 = homozygous alternate, MIN_VALUE = missing
                int genotypeIndex = markerMatrix[sIdx][tIdx];
                byte homoRef = NucleotideAlignmentConstants.getNucleotideDiploidByte(variants[0]);
                byte heterozygous = NucleotideAlignmentConstants.getNucleotideDiploidByte(variants[0] + variants[1]);
                byte homoAlt = NucleotideAlignmentConstants.getNucleotideDiploidByte(variants[1]);

                if (genotypeIndex == 0) {
                    genotypeCallTableBuilder.setBase(tIdx, sIdx, homoRef);
                } else if (genotypeIndex == 1) {
                    genotypeCallTableBuilder.setBase(tIdx, sIdx, heterozygous);
                } else if (genotypeIndex == 2) {
                    genotypeCallTableBuilder.setBase(tIdx, sIdx, homoAlt);
                } else if (genotypeIndex == Integer.MIN_VALUE) {
                    genotypeCallTableBuilder.setBase(tIdx, sIdx, GenotypeTable.UNKNOWN_DIPLOID_ALLELE);
                } else {
                    throw new IllegalArgumentException("createGenotypeFromRDataFrameElements: genotype index must be 0, 1, 2, or NA (MIN_VALUE)");
                }

            }

        }

        PositionList positionList = positionListBuilder.build();

        return GenotypeTableBuilder.getInstance(genotypeCallTableBuilder.build(), positionList, taxaList);

    }

    public static Phenotype createPhenotypeFromRDataFrameElements(String[] taxaArray, String[] colNames, String[] attributeType, List dataVectors) {
        if (colNames.length != attributeType.length || colNames.length != dataVectors.size())
            throw new IllegalArgumentException("ColNames, attributeType, and dataVectors need to be same size");

        List<Taxon> taxaList = Stream.of(taxaArray).map(Taxon::new).collect(Collectors.toList());

        List<PhenotypeAttribute> attributes = new ArrayList<PhenotypeAttribute>();
        List<Phenotype.ATTRIBUTE_TYPE> types = new ArrayList<Phenotype.ATTRIBUTE_TYPE>();
        TaxaAttribute ta = new TaxaAttribute(taxaList);
        attributes.add(ta);
        types.add(Phenotype.ATTRIBUTE_TYPE.taxa);

        for (int i = 0; i < colNames.length; i++) {
            Object o = dataVectors.get(i);
            if (o instanceof double[]) {
                attributes.add(new NumericAttribute(colNames[i], (double[]) o));
                Phenotype.ATTRIBUTE_TYPE attribute_type = Phenotype.ATTRIBUTE_TYPE.valueOf(attributeType[i]);
                types.add(attribute_type);
            } else if (o instanceof int[]) {
                //INT min is NA
                int[] initialTraits = (int[]) o;
                Phenotype.ATTRIBUTE_TYPE attribute_type = Phenotype.ATTRIBUTE_TYPE.valueOf(attributeType[i]);
                if (attribute_type == Phenotype.ATTRIBUTE_TYPE.data || attribute_type == Phenotype.ATTRIBUTE_TYPE.covariate) {
                    float[] traitsValueWithNaN = new float[initialTraits.length];
                    for (int j = 0; j < initialTraits.length; j++)
                        traitsValueWithNaN[j] = (initialTraits[j] == Integer.MIN_VALUE) ? Float.NaN : initialTraits[j];
                    attributes.add(new NumericAttribute(colNames[i], traitsValueWithNaN));
                    types.add(attribute_type);
                } else if (attribute_type == Phenotype.ATTRIBUTE_TYPE.factor) {
                    String[] categories = Arrays.stream(initialTraits).mapToObj(intval -> Integer.toString(intval)).toArray(String[]::new);
                    attributes.add(new CategoricalAttribute(colNames[i], categories));
                    types.add(attribute_type);
                } else {
                    throw new IllegalArgumentException("attribute type of " + attributeType[i] + "inconsistent with data type of int for " + colNames[i]);
                }
            } else if (o instanceof String[]) {
                attributes.add(new CategoricalAttribute(colNames[i], (String[]) o));
                Phenotype.ATTRIBUTE_TYPE attribute_type = Phenotype.ATTRIBUTE_TYPE.valueOf(attributeType[i]);
                types.add(attribute_type);
            } else {
                throw new IllegalArgumentException("Unsupported type for phenotype table");
            }
            System.out.println("phenotype column added: " + colNames[i] + "," + attributeType[i]);
        }

        return new PhenotypeBuilder().fromAttributeList(attributes, types).build().get(0);
    }

    //    traitName traitType       traitAttribute isReponse isPredictor   newType effect
    // 1      Taxa      taxa        TaxaAttribute     FALSE       FALSE      <NA>   <NA>
    // 2  location    factor CategoricalAttribute     FALSE        TRUE    factor  fixed
    // 3     EarHT      data     NumericAttribute      TRUE       FALSE      data   <NA>
    // 4     dpoll      data     NumericAttribute     FALSE       FALSE      <NA>   <NA>
    // 5    EarDia      data     NumericAttribute      TRUE       FALSE      data   <NA>
    // 6        Q1 covariate     NumericAttribute     FALSE        TRUE covariate  fixed
    // 7        Q2 covariate     NumericAttribute     FALSE        TRUE covariate  fixed
    // 8        Q3 covariate     NumericAttribute     FALSE        TRUE covariate  fixed
    // 9         G  genotype             Genotype     FALSE        TRUE  genotype  fixed

    public static Map<String, Object> association(DistanceMatrix kinship, GenotypeTable genotype, Phenotype phenotype, GenotypePhenotype genoPheno, int minClassSize, boolean biallelicOnly, boolean appendAddDom, boolean saveToFile, String outputFile, double maxP) {

        String timeStr = LocalDateTime.now().format(DateTimeFormatter.ofPattern("MMM d, uuuu H:mm:s"));
        myLogger.info("Starting association: time: " + timeStr);

        try {

            if (genotype == null && phenotype == null) {
                myLogger.warn("association: genotype and phenotype are null.  Nothing calculated");
            }

            if (kinship == null) {

                myLogger.info("association: running GLM");

                FixedEffectLMPlugin plugin = new FixedEffectLMPlugin(null, false);
                plugin.biallelicOnly(biallelicOnly);
                plugin.minClassSize(minClassSize);
                plugin.appendAddDom(appendAddDom);
                plugin.saveAsFile(saveToFile);
                plugin.siteReportFilename(outputFile + "_site");
                plugin.alleleReportFilename(outputFile + "_allele");
                plugin.bluesReportFilename(outputFile + "_blues");
                plugin.anovaReportFilename(outputFile + "_anova");
                plugin.maxPvalue(maxP);

                DataSet input;

                if (genotype == null) {
                    plugin.phenotypeOnly(true);
                    input = DataSet.getDataSet(phenotype);
                } else {
                    plugin.phenotypeOnly(false);
                    input = DataSet.getDataSet(genoPheno);
                }

                DataSet output = plugin.performFunction(input);

                return tableReportsMap(output);

            } else {

                myLogger.info("association: running MLM");

                MLMPlugin plugin = new MLMPlugin(null, false);
                plugin.setWriteOutputToFile(saveToFile);
                if (outputFile != null && !outputFile.isEmpty()) plugin.setOutputName(outputFile);

                Datum genoDatum = new Datum("GenotypePhenotype", genoPheno, null);
                Datum kinshipDatum = new Datum("Kinship", kinship, null);
                DataSet input = new DataSet(new Datum[]{genoDatum, kinshipDatum}, null);

                DataSet output = plugin.performFunction(input);

                return tableReportsMap(output);

            }

        } finally {
            timeStr = LocalDateTime.now().format(DateTimeFormatter.ofPattern("MMM d, uuuu H:mm:s"));
            myLogger.info("Finished association: time: " + timeStr);
        }

    }

    public static Map<String, Object> fastAssociation(GenotypePhenotype genoPheno, double maxp, Integer maxThreads, boolean writeToFile, String outputFile) {

        String timeStr = LocalDateTime.now().format(DateTimeFormatter.ofPattern("MMM d, uuuu H:mm:s"));
        myLogger.info("Starting fastAssociation: time: " + timeStr);

        try {

            if (genoPheno == null) {
                myLogger.warn("fastAssociation: GenotypePhenotype is null.  Nothing calculated");
            }

            FastMultithreadedAssociationPlugin plugin = new FastMultithreadedAssociationPlugin(null, false);

            plugin.maxp(maxp);

            if (maxThreads != null) plugin.maxThreads(maxThreads);

            plugin.saveAsFile(writeToFile);
            plugin.reportFilename(outputFile);

            DataSet input = DataSet.getDataSet(genoPheno);

            DataSet output = plugin.performFunction(input);

            return tableReportsMap(output);

        } finally {
            timeStr = LocalDateTime.now().format(DateTimeFormatter.ofPattern("MMM d, uuuu H:mm:s"));
            myLogger.info("Finished fastAssociation: time: " + timeStr);
        }

    }

    private static Map<String, Object> tableReportsMap(DataSet output) {

        Map<String, Object> result = new HashMap<>();

        if (output == null) return result;

        for (Datum temp : output.getDataOfType(TableReport.class)) {
            String name = temp.getName();

            if (name.startsWith("GLM_Genotypes")) {
                result.put("GLM_Genotypes", temp.getData());
            } else if (name.startsWith("GLM_Stats")) {
                result.put("GLM_Stats", temp.getData());
            } else if (name.startsWith("MLM_statistics")) {
                result.put("MLM_Stats", temp.getData());
            } else if (name.startsWith("MLM_effects")) {
                result.put("MLM_Effects", temp.getData());
            } else if (name.startsWith("Residuals for")) {
                String[] tokens = name.split(" ");
                if (tokens[2].endsWith(".")) {
                    result.put("MLM_Residuals_" + tokens[2].substring(0, tokens[2].length() - 1), temp.getData());
                } else {
                    result.put("MLM_Residuals_" + tokens[2], temp.getData());
                }
            } else if (name.startsWith("MLM_compression")) {
                result.put("MLM_Compression", temp.getData());
            } else if (name.startsWith("BLUEs")) {
                result.put("BLUE", temp.getData());
            } else if (name.startsWith("Phenotype_ANOVA")) {
                result.put("BLUE_ANOVA", temp.getData());
            } else if (name.startsWith("Fast Association")) {
                result.put("FastAssociation", temp.getData());
            }
        }

        return result;

    }

    public static TableReport linkageDiseq(GenotypeTable genotype, String ldType, int windowSize, String hetTreatment) {

        LinkageDisequilibrium.testDesign testDesign = null;
        try {
            testDesign = LinkageDisequilibrium.testDesign.valueOf(ldType);
        } catch (Exception e) {
            myLogger.debug(e.getMessage(), e);
            myLogger.error("linkageDiseq: ldType: " + ldType + " is unknown");
            throw new IllegalArgumentException("GenerateRCode: linkageDiseq: ldType: " + ldType + " is unknown");
        }

        LinkageDisequilibrium.HetTreatment treatment = null;
        if (hetTreatment.equalsIgnoreCase("ignore")) {
            treatment = LinkageDisequilibrium.HetTreatment.Haplotype;
        } else if (hetTreatment.equalsIgnoreCase("missing")) {
            treatment = LinkageDisequilibrium.HetTreatment.Homozygous;
        } else if (hetTreatment.equalsIgnoreCase("third")) {
            treatment = LinkageDisequilibrium.HetTreatment.Genotype;
        } else {
            myLogger.error("linkageDiseq: unknown LD Type: " + hetTreatment);
            throw new IllegalArgumentException("GenerateRCode: linkageDiseq: unknown LD Type: " + hetTreatment);
        }

        LinkageDisequilibrium result = new LinkageDisequilibrium(genotype, windowSize, testDesign, -1, null, false, 100, null, treatment);
        result.run();

        return result;

    }

    public static void exportToFlapjack(GenotypeTable genotype, String filename) {
        FlapjackUtils.writeToFlapjack(genotype, filename, '\t');
    }

    public static GenotypeTable read(String filename, boolean keepDepth, boolean sortPositions) {
        return ImportUtils.read(filename, keepDepth, sortPositions);
    }

    public static TableReport genomicSelection(Phenotype phenotype, DistanceMatrix matrix, boolean doCV, int kFolds, int nIter) {

        GenomicSelectionPlugin plugin = new GenomicSelectionPlugin(null, false);
        plugin.performCrossValidation(doCV);
        plugin.kFolds(kFolds);
        plugin.nIterations(nIter);

        Datum phenoInput = new Datum("phenotype", phenotype, null);
        Datum matrixInput = new Datum("matrix", matrix, null);
        DataSet result = plugin.performFunction(new DataSet(new Datum[]{phenoInput, matrixInput}, null));

        return (TableReport) result.getDataOfType(TableReport.class).get(0).getData();

    }

    /**
     * Filters given genotypes to only the ranges specified.
     *
     * @param input genotype table to be filtered
     * @param seqName sequence names (i.e. chromosomes) for ranges
     * @param start start positions for each range (inclusive)
     * @param end end positions for each range (inclusive)
     *
     * @return filtered genotype table
     */
    public static GenotypeTable filterSitesByGRanges(GenotypeTable input, String[] seqName, int[] start, int[] end) {

        int numOfRanges = seqName.length;

        if (numOfRanges != start.length || numOfRanges != end.length) {
            throw new IllegalArgumentException("GenerateRCode: filterSitesByGRanges: arrays seqName, start, and end must be same length");
        }

        int numSites = input.numberOfSites();
        BitSet sitesToInclude = new OpenBitSet(numSites);
        for (int index = 0; index < numOfRanges; index++) {

            if (start[index] > end[index]) {
                throw new IllegalArgumentException("GenerateRCode: filterSitesByGRanges: index: " + index + " start position: " + start[index] + " is greater than end position: " + end[index]);
            }

            int startSite = input.siteOfPhysicalPosition(start[index], Chromosome.instance(seqName[index]));
            // if startSite is negative, the position wasn't found and
            // (-startSite - 1) is the insertion location.
            if (startSite < 0) {
                startSite = -startSite - 1;
            }

            int endSite = input.siteOfPhysicalPosition(end[index], Chromosome.instance(seqName[index]));
            // if endSite is negative, the positon wasn't found and
            // (-endSite - 1) is the insertion location. Another 1 is subtracted
            // to move back to the last site included in the range.
            if (endSite < 0) {
                endSite = -endSite - 2;
            }
            
            for (int i = startSite; i <= endSite; i++) {
                sitesToInclude.fastSet(i);
            }
        }

        int numNewSites = (int) sitesToInclude.cardinality();
        int[] result = new int[numNewSites];
        int count = 0;
        for (int s = 0; s < numSites; s++) {
            if (sitesToInclude.fastGet(s)) {
                result[count++] = s;
            }
        }

        return FilterGenotypeTable.getInstance(input, result);

    }

    /**
     * Convert an R matrix to a DistanceMatrix object
     *
     * @param m A numeric R matrix passed as a 2d array of doubles
     * @param taxa A String array of taxa IDs
     *
     * @return DistanceMatrix object
     */
    public static DistanceMatrix asTasselDistanceMatrix(double[][] m, String[] taxa) {

        if (taxa.length != m.length) {
            throw new IllegalArgumentException("GenerateRCode: Number of taxa does not equal dimensions of matrix.");
        }

        DistanceMatrixBuilder distBuilder = DistanceMatrixBuilder.getInstance(m.length);
        for (int i = 0; i < m.length; i++) {
            for (int j = 0; j < m[0].length; j++) {
                distBuilder.set(j, i, m[i][j]);
            }
            distBuilder.addTaxon(new Taxon(taxa[i]));
        }
        return distBuilder.build();
    }

}

package net.maizegenetics.analysis.association;

import java.awt.Frame;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.TimeUnit;
import java.util.stream.Collectors;

import javax.swing.ImageIcon;

import net.maizegenetics.stats.linearmodels.*;
import org.apache.commons.math3.distribution.FDistribution;
import org.apache.commons.math3.exception.OutOfRangeException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import net.maizegenetics.dna.map.Position;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableUtils;
import net.maizegenetics.dna.snp.GenotypeTable.GENOTYPE_TABLE_COMPONENT;
import net.maizegenetics.phenotype.CategoricalAttribute;
import net.maizegenetics.phenotype.GenotypePhenotype;
import net.maizegenetics.phenotype.NumericAttribute;
import net.maizegenetics.phenotype.Phenotype;
import net.maizegenetics.phenotype.PhenotypeAttribute;
import net.maizegenetics.phenotype.Phenotype.ATTRIBUTE_TYPE;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.GeneratePluginCode;
import net.maizegenetics.plugindef.PluginParameter;
import net.maizegenetics.prefs.TasselPrefs;
import net.maizegenetics.util.TableReport;
import net.maizegenetics.util.TableReportBuilder;

public class FastMultithreadedAssociationPlugin extends AbstractPlugin {
    private static Logger myLogger = LogManager.getLogger(FastMultithreadedAssociationPlugin.class);
    private GENOTYPE_TABLE_COMPONENT[] GENOTYPE_COMP = new GENOTYPE_TABLE_COMPONENT[] {
            GENOTYPE_TABLE_COMPONENT.Genotype, GENOTYPE_TABLE_COMPONENT.ReferenceProbability,
            GENOTYPE_TABLE_COMPONENT.AlleleProbability };

    private final byte NN = GenotypeTable.UNKNOWN_DIPLOID_ALLELE;
    private Phenotype myPhenotype;
    private GenotypeTable myGenotype;
    List<String >phenotypeNames;
    double minR2;
    private FDistribution Fdist;
    GenotypePhenotype myGenoPheno;
    
    //plugin parameter definitions
    private PluginParameter<Double> maxp =
            new PluginParameter.Builder<>("MaxPValue", .001, Double.class)
                    .guiName("MaxPValue")
                    .description("The maximum p-value that will be output by the analysis.")
                    .build();
    private PluginParameter<GENOTYPE_TABLE_COMPONENT> myGenotypeTable =
            new PluginParameter.Builder<>("genotypeComponent", GENOTYPE_TABLE_COMPONENT.Genotype, GENOTYPE_TABLE_COMPONENT.class)
                    .genotypeTable()
                    .range(GENOTYPE_COMP)
                    .description("If the genotype table contains more than one type of genotype data, choose the type to use for the analysis.")
                    .build();
    private PluginParameter<Boolean> saveAsFile =
            new PluginParameter.Builder<>("writeToFile", false, Boolean.class)
                    .description("Should the results be saved to a file rather than stored in memory? It true, the results will be written to a file as each SNP is analyzed in order to reduce memory requirements"
                            + "and the results will NOT be saved to the data tree. Default = false.")
                    .guiName("Write to file")
                    .build();
    private PluginParameter<String> reportFilename =
            new PluginParameter.Builder<>("outputFile", null, String.class)
                    .outFile()
                    .dependentOnParameter(saveAsFile)
                    .description("The name of the file to which these results will be saved.")
                    .guiName("Output File")
                    .build();
    private PluginParameter<Integer> maxThreads = new PluginParameter.Builder<>("maxThreads", TasselPrefs.getMaxThreads(), Integer.class)
    		.description("the maximum number of threads to be used by this plugin.")
    		.guiName("Max Threads")
    		.build();
    
    public FastMultithreadedAssociationPlugin() {
        this(null, false);
    }

    public FastMultithreadedAssociationPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    @Override
    protected void preProcessParameters(DataSet input) {
        List<Datum> inData = input.getDataOfType(GenotypePhenotype.class);
        if (inData.size() != 1) throw new IllegalArgumentException("Fast Association requires exactly one joined genotype-phenotype data set.");
    }

    @Override
    public DataSet processData(DataSet input) {
        long start = System.currentTimeMillis();
        
        int maxSitesInQueue = 2000;
        int maxObjectsInQueue = 1000;
        
        Datum inDatum = input.getDataOfType(GenotypePhenotype.class).get(0);
        myGenoPheno = (GenotypePhenotype) inDatum.getData();
        myGenotype = myGenoPheno.genotypeTable();
        myPhenotype = myGenoPheno.phenotype();
        int numberOfObservations = myPhenotype.numberOfObservations();
        testMissingDataInTheBaseModel();

        //calculate orthogonal phenotypes
        SolveByOrthogonalizing sbo = initializeOrthogonalizer();
        
        //determine errdf
        double errdf = numberOfObservations - sbo.baseDf() - 1;
        Fdist = new FDistribution(1, errdf);
        
        //calculate minR2
        calculateR2Fromp(errdf);
        
        //initialize report builder
        TableReportBuilder myReport = initializeOutput(inDatum);
        
        //create a thread pool
        int nthreads = maxThreads.value();
        nthreads = Math.max(nthreads, 2);
        int siteTesterThreads = nthreads - 1;
        ExecutorService myExecutor = Executors.newFixedThreadPool(nthreads);
                
        //start report thread
        BlockingQueue<Object[]> reportQueue = new LinkedBlockingQueue<>();
        
        //start processing and output threads
        BlockingQueue<Marker> siteQueue = new LinkedBlockingQueue<>(maxSitesInQueue);
        List<double[]> dataList = sbo.getOrthogonalizedData();
        List<double[]> uList = sbo.getUColumns();
        
        for (int i = 0; i < siteTesterThreads; i++) {
            myExecutor.execute(new SiteTester(dataList, phenotypeNames, uList, siteQueue, reportQueue, minR2, errdf, numberOfObservations));
            //the next line can be used to test whether each thread should have its own copy of the phenotype data
//            myExecutor.execute(new SiteTester(sbo.copyOrthogonalizedData(), phenotypeNames, sbo.copyUColumns(), siteQueue, reportQueue, minR2, errdf, numberOfObservations));
        }
        
        //start the reporter
        myExecutor.execute(new ReportWriter(myReport, reportQueue, siteTesterThreads));
        
        System.out.printf("Time to set up threads = %d ms.\n", System.currentTimeMillis() - start);
        start = System.currentTimeMillis();
        
        //add sites to the siteQueue
        int nsites = myGenotype.numberOfSites();
        System.out.printf("myGenotype has %d sites\n", nsites);
        for (int s = 0; s < nsites; s++) {
            if (s % 1000000 == 0) myLogger.info("Adding site " + s + " to the site queue.");
            try {
                byte major = myGenotype.majorAllele(s);
                double freq = myGenotype.majorAlleleFrequency(s);
                byte[] geno = myGenoPheno.genotypeAllTaxa(s);
                siteQueue.put(new Marker(geno, major, freq, myGenotype.positions().get(s)));
            } catch (Exception e) {
                throw new RuntimeException("Site thread interrupted at site " + s, e);
            }
        }
        
        //add end signal times number of threads to queue
        for (int i = 0; i < nthreads; i++) {
            byte zerobyte = (byte) 0;
            try {
                siteQueue.put(new Marker(new byte[0], zerobyte, 0.0, null));
            } catch (InterruptedException e) {
                throw new RuntimeException("siteQueue interrupted", e);
            }
        }
        
        //wait here until all threads finish
        myExecutor.shutdown();
        try {
            myExecutor.awaitTermination(1, TimeUnit.HOURS);
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
        
        System.out.printf("Time to process sites = %d ms.\n", System.currentTimeMillis() - start);
        if (saveAsFile.value()) {
            myReport.build();
            return null;
        }
        
        String name = String.format("Fast Association_%s", inDatum.getName());
        String comment = String.format("Fast Association Test Results\n Source = %s", inDatum.getName());
        return new DataSet(new Datum(name, myReport.build(), comment), this);
    }

    private void testMissingDataInTheBaseModel() {
        for (PhenotypeAttribute attr : myPhenotype.attributeListOfType(ATTRIBUTE_TYPE.factor)) {
            if (attr.missing().cardinality() > 0) {
                String msg = "There is missing data in the factor " + attr.name();
                throw new IllegalArgumentException(msg);
            }
        }
        for (PhenotypeAttribute attr : myPhenotype.attributeListOfType(ATTRIBUTE_TYPE.covariate)) {
            if (attr.missing().cardinality() > 0) {
                String msg = "There is missing data in the covariate " + attr.name();
                throw new IllegalArgumentException(msg);
            }
        }
        for (PhenotypeAttribute attr : myPhenotype.attributeListOfType(ATTRIBUTE_TYPE.data)) {
            if (attr.missing().cardinality() > 0) {
                String msg = "There is missing data in the phenotype " + attr.name();
                throw new IllegalArgumentException(msg);
            }
        }
    }
    
    private SolveByOrthogonalizing initializeOrthogonalizer() {
        List<PhenotypeAttribute> phenotypeList =
                myPhenotype.attributeListOfType(ATTRIBUTE_TYPE.data);
        List<PhenotypeAttribute> covariateList =
                myPhenotype.attributeListOfType(ATTRIBUTE_TYPE.covariate);
        List<PhenotypeAttribute> factorList =
                myPhenotype.attributeListOfType(ATTRIBUTE_TYPE.factor);

        //build the model, no mean necessary because it will not be used
        List<ModelEffect> baseModel = new ArrayList<>();
        for (PhenotypeAttribute pa : factorList) {
            CategoricalAttribute ca = (CategoricalAttribute) pa;
            baseModel.add(new FactorModelEffect(ca.allIntValues(), true, ca.name()));
        }
        for (PhenotypeAttribute pa : covariateList) {
            NumericAttribute na = (NumericAttribute) pa;
            CovariateModelEffect cme =
                    new CovariateModelEffect(AssociationUtils.convertFloatArrayToDouble(na.floatValues()), na.name());
            baseModel.add(cme);
        }

        List<double[]> dataList = phenotypeList.stream()
                .map(pa -> (float[]) pa.allValues())
                .map(a -> AssociationUtils.convertFloatArrayToDouble(a))
                .collect(Collectors.toList());

        phenotypeNames =
                phenotypeList.stream().map(PhenotypeAttribute::name).collect(Collectors.toList());

        return SolveByOrthogonalizing.getInstanceFromModel(baseModel, dataList);
    }

    private TableReportBuilder initializeOutput(Datum myDatum) {
        //output is a TableReport with p-value; site position information: chr, position, id; trait name
        //add separate values for additive test and dominant test later
        String[] columnNames = new String[] { AssociationConstants.STATS_HEADER_TRAIT,
                AssociationConstants.STATS_HEADER_MARKER,
                AssociationConstants.STATS_HEADER_CHR,
                AssociationConstants.STATS_HEADER_POSITION,
                "df",
                "r2",
                AssociationConstants.STATS_HEADER_P_VALUE };
        String name = "EqtlReport_" + myDatum.getName();
        if (saveAsFile.value())
            return TableReportBuilder.getInstance(name, columnNames, reportFilename.value());
        else
            return TableReportBuilder.getInstance(name, columnNames);
    }

    private void calculateR2Fromp(double errdf) {
        //returns the value of R^2 corresponding to the value of F, f for which P(F>f) = alpha
        double p = 1 - maxp.value();
        try {
            double F = Fdist.inverseCumulativeProbability(p);
            minR2 = F / (errdf + F);
        } catch (OutOfRangeException e) {
            e.printStackTrace();
            minR2 = Double.NaN;
        }
    }

    @Override
    public ImageIcon getIcon() {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public String getButtonName() {
        return "Fast-MT Association";
    }

    @Override
    public String getToolTipText() {
        return "Multi-threaded version of Fast Association";
    }

    class SiteTester extends Thread {
        
        final List<double[]> orthogonalPhenotypes;
        final List<String> phenotypeNames;
        final List<double[]> Ucolumns;
        final BlockingQueue<Marker> siteQueue;
        final BlockingQueue<Object[]> outQueue;
        final double minR2;
        final int nphenotypes;
        final int numberOfObservations;
        final double errdf;
        final private FDistribution Fdist;
        
        SiteTester(List<double[]> orthogonalPhenotypes, List<String> phenotypeNames, List<double[]> Ucol, BlockingQueue<Marker> siteQueue, BlockingQueue<Object[]> outQueue, double minRSquare, double errdf, int ntaxa) {
            this.orthogonalPhenotypes = orthogonalPhenotypes;
            this.phenotypeNames = phenotypeNames;
            Ucolumns = Ucol;
            this.siteQueue = siteQueue;
            this.outQueue = outQueue;
            minR2 = minRSquare;
            numberOfObservations = ntaxa;
            this.errdf = errdf;
            nphenotypes = orthogonalPhenotypes.size();
            Fdist = new FDistribution(1, errdf);
        }
        
        public void run() {
            try {
                Marker thisMarker = siteQueue.poll(30, TimeUnit.SECONDS);
                if (thisMarker == null) {
                    //send end signal to reporter
                    outQueue.put(new Object[0]);
                    throw new IllegalStateException("ERROR: The site tester timeout was exceeded.");
                }
                byte[] geno = thisMarker.geno;
                while (geno.length > 0) {
                    byte major = thisMarker.major;
                    double genoMean = thisMarker.majorFrequency;
                    
                    //convert genotypes to centered values
                    double[] siteValues = new double[numberOfObservations];
                    for (int t = 0; t < numberOfObservations; t++) {
                        if (geno[t] == NN) siteValues[t] = 0;
                        else {
                            siteValues[t] = -genoMean;
                            byte[] alleles = GenotypeTableUtils.getDiploidValues(geno[t]);
                            if (alleles[0] == major) siteValues[t] += 0.5;
                            if (alleles[1] == major) siteValues[t] += 0.5;
                        }
                    }

                    //debug
                    double sum = 0;
                    for (double d : siteValues) sum += d;
                    sum /= numberOfObservations;
                    
                    siteValues = orthogonalizeByBase(siteValues);
                    siteValues = SolveByOrthogonalizing.centerAndScale(siteValues);
                    
                    if (siteValues == null) {
                        System.err.printf("siteValues null at position %d, probably invariant\n", thisMarker.myPosition.getPosition());
                    } else {
                        double[] r2values = new double[nphenotypes];
                        for (int p = 0; p < nphenotypes; p++) {
                            double sumprod = 0;
                            double[] pheno = orthogonalPhenotypes.get(p);
                            for (int t = 0; t < numberOfObservations; t++) sumprod += siteValues[t] * pheno[t];
                            r2values[p] = sumprod * sumprod;
                        }
                        
                        outputResult(r2values, thisMarker.myPosition);
                    }

                    thisMarker = siteQueue.poll(2, TimeUnit.SECONDS);
                    if (thisMarker == null) {
                        //send end signal to reporter
                        outQueue.put(new Object[0]);
                        throw new IllegalStateException("Error: The site tester timeout was exceeded.");
                    }

                    geno = thisMarker.geno;
                }
                //send end signal to reporter
                outQueue.put(new Object[0]);
                
            } catch (InterruptedException e) {
                throw new RuntimeException("InterruptedException occurred in SiteTester thread", e);
            }
        }
        
        private double[] orthogonalizeByBase(double[] vector) {
            if (Ucolumns == null || Ucolumns.size() == 0) return vector;
            int nrows = vector.length;
            double[] result = Arrays.copyOf(vector, nrows);

            for (double[] u : Ucolumns) {
                double ip = SolveByOrthogonalizing.innerProduct(vector, u);
                for (int j = 0; j < nrows; j++)
                    result[j] -= ip * u[j];
            }

            return result;
        }

        private void outputResult(double[] rvalues, Position pos) throws InterruptedException {
            //for r2 values >= minR2, create an output record and add it to the output queue
        	for (int p = 0; p < nphenotypes; p++) {
        		if (rvalues[p] >= minR2) {
        		    Object[] result = new Object[]{ phenotypeNames.get(p), pos.getSNPID(),
                                pos.getChromosome().getName(), pos.getPosition(),
                                1, rvalues[p],
                                pvalue(rvalues[p])};
        		    outQueue.put(result);
        		}
        	}
        }
        
        private double pvalue(double rvalue) {
            double F = rvalue / (1 - rvalue) * errdf;
            double p;
            try {
                p = LinearModelUtils.Ftest(F, 1, errdf);
            } catch (Exception e) {
                p = Double.NaN;
            }
            return p;
        }

    }
    
    class ReportWriter extends Thread {
        TableReportBuilder myReportBuilder;
        BlockingQueue<Object[]> myReportQueue;
        int numberOfSources;
        
        ReportWriter(TableReportBuilder reportBuilder, BlockingQueue<Object[]> reportQueue, int numberOfSourceThreads) {
            myReportBuilder = reportBuilder;
            myReportQueue = reportQueue;
            numberOfSources = numberOfSourceThreads;
        }
        
        @Override
        public void run() {
            int numberOfFinishedThreads = 0;
            try {
                do {
                    Object[] reportRow = myReportQueue.poll(30, TimeUnit.MINUTES);
                    if (reportRow == null) {
                        throw new IllegalStateException("ERROR: report queue timed out.");
                    }
                    if (reportRow.length > 0) {
                        myReportBuilder.add(reportRow);
                    }
                    else {
                        numberOfFinishedThreads++;
                        System.out.printf("number of threads finished = %d\n", numberOfFinishedThreads);
                    }
                    
                } while (numberOfFinishedThreads < numberOfSources);
                
            } catch (InterruptedException e) {
                throw new RuntimeException("Report thread was interrupted.", e);
            }
            System.out.println("report thread finished");
        }
    }
    
    class Marker {
    	byte[] geno;
    	byte major;
    	double majorFrequency;
    	Position myPosition;
    	
    	Marker(byte[] geno, byte major, double majorFreq, Position pos) {
    		this.geno = geno;
    		this.major = major;
    		majorFrequency = majorFreq;
    		myPosition = pos;
    	}
    	
    }
    
    // The following getters and setters were auto-generated.
    // Please use this method to re-generate.
    //
    // public static void main(String[] args) {
    //     GeneratePluginCode.generate(FastMultithreadedAssociationPlugin.class);
    // }

    /**
     * Convenience method to run plugin with one return object.
     */
    public TableReport runPlugin(DataSet input) {
        return (TableReport) performFunction(input).getData(0).getData();
    }

    /**
     * The maximum p-value that will be output by the analysis.
     *
     * @return MaxPValue
     */
    public Double maxp() {
        return maxp.value();
    }

    /**
     * Set MaxPValue. The maximum p-value that will be output
     * by the analysis.
     *
     * @param value MaxPValue
     *
     * @return this plugin
     */
    public FastMultithreadedAssociationPlugin maxp(Double value) {
        maxp = new PluginParameter<>(maxp, value);
        return this;
    }

    /**
     * If the genotype table contains more than one type of
     * genotype data, choose the type to use for the analysis.
     *
     * @return Genotype Component
     */
    public GENOTYPE_TABLE_COMPONENT genotypeTable() {
        return myGenotypeTable.value();
    }

    /**
     * Set Genotype Component. If the genotype table contains
     * more than one type of genotype data, choose the type
     * to use for the analysis.
     *
     * @param value Genotype Component
     *
     * @return this plugin
     */
    public FastMultithreadedAssociationPlugin genotypeTable(GENOTYPE_TABLE_COMPONENT value) {
        myGenotypeTable = new PluginParameter<>(myGenotypeTable, value);
        return this;
    }

    /**
     * Should the results be saved to a file rather than stored
     * in memory? It true, the results will be written to
     * a file as each SNP is analyzed in order to reduce memory
     * requirementsand the results will NOT be saved to the
     * data tree. Default = false.
     *
     * @return Write to file
     */
    public Boolean saveAsFile() {
        return saveAsFile.value();
    }

    /**
     * Set Write to file. Should the results be saved to a
     * file rather than stored in memory? It true, the results
     * will be written to a file as each SNP is analyzed in
     * order to reduce memory requirementsand the results
     * will NOT be saved to the data tree. Default = false.
     *
     * @param value Write to file
     *
     * @return this plugin
     */
    public FastMultithreadedAssociationPlugin saveAsFile(Boolean value) {
        saveAsFile = new PluginParameter<>(saveAsFile, value);
        return this;
    }

    /**
     * The name of the file to which these results will be
     * saved.
     *
     * @return Output File
     */
    public String reportFilename() {
        return reportFilename.value();
    }

    /**
     * Set Output File. The name of the file to which these
     * results will be saved.
     *
     * @param value Output File
     *
     * @return this plugin
     */
    public FastMultithreadedAssociationPlugin reportFilename(String value) {
        reportFilename = new PluginParameter<>(reportFilename, value);
        return this;
    }

    /**
     * the maximum number of threads to be used by this plugin.
     *
     * @return Max Threads
     */
    public Integer maxThreads() {
        return maxThreads.value();
    }

    /**
     * Set Max Threads. the maximum number of threads to be
     * used by this plugin.
     *
     * @param value Max Threads
     *
     * @return this plugin
     */
    public FastMultithreadedAssociationPlugin maxThreads(Integer value) {
        maxThreads = new PluginParameter<>(maxThreads, value);
        return this;
    }

}

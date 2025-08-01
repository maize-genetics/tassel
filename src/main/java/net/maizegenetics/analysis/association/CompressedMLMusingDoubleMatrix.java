package net.maizegenetics.analysis.association;

import net.maizegenetics.analysis.data.ExportPlugin;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableUtils;
import net.maizegenetics.dna.snp.NucleotideAlignmentConstants;
import net.maizegenetics.dna.snp.score.SiteScore.SITE_SCORE_TYPE;
import net.maizegenetics.phenotype.GenotypePhenotype;
import net.maizegenetics.phenotype.NumericAttribute;
import net.maizegenetics.phenotype.Phenotype;
import net.maizegenetics.phenotype.PhenotypeAttribute;
import net.maizegenetics.phenotype.PhenotypeBuilder;
import net.maizegenetics.phenotype.TaxaAttribute;
import net.maizegenetics.phenotype.Phenotype.ATTRIBUTE_TYPE;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.matrixalgebra.Matrix.DoubleMatrix;
import net.maizegenetics.matrixalgebra.Matrix.DoubleMatrixFactory;
import net.maizegenetics.matrixalgebra.Matrix.DoubleMatrixFactory.FactoryType;
import net.maizegenetics.taxa.distance.DistanceMatrix;
import net.maizegenetics.util.BitSet;
import net.maizegenetics.util.OpenBitSet;
import net.maizegenetics.util.TableReport;
import net.maizegenetics.util.TableReportBuilder;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListBuilder;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.taxa.tree.UPGMATree;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.stats.EMMA.EMMAforDoubleMatrix;
import net.maizegenetics.stats.linearmodels.FactorModelEffect;
import net.maizegenetics.stats.linearmodels.LinearModelUtils;
import net.maizegenetics.stats.linearmodels.ModelEffectUtils;
import net.maizegenetics.stats.linearmodels.SweepFast;
import net.maizegenetics.stats.linearmodels.SymmetricMatrixInverterDM;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;
import net.maizegenetics.taxa.distance.DistanceMatrixBuilder;

public class CompressedMLMusingDoubleMatrix {

    private static final Logger myLogger = LogManager.getLogger(CompressedMLMusingDoubleMatrix.class);
    private static final List<String> homGenotypes = Arrays.asList("A","C","G","T","Z");
    private static final List<String> hetGenotypes = Arrays.asList("R","W","K","Y","S","M","0");
    private final boolean useCompression;
    private final boolean useP3D;
    private final double compression;
    private boolean outputResiduals = true;
    
    private final GenotypePhenotype myGenoPheno;
    private final Phenotype myPhenotype;
    private final GenotypeTable myGenotype;
    private final boolean hasGenotype;
    private final MLMPlugin parentPlugin;
    private final DistanceMatrix kinshipMatrix;
    private double resvar, genvar, lnlk;
    private boolean testMarkers = true;
    private SymmetricMatrixInverterDM Vminus = null;
    private String datasetName;
	private List<PhenotypeAttribute> factorAttributeList;
	private List<PhenotypeAttribute> covariateAttributeList;

    private final Phenotype myWeightMatrix;
    
    
    private final TableReportBuilder siteReportBuilder;
    private final TableReportBuilder alleleReportBuilder;
    private final TableReportBuilder compressionReportBuilder;
    
    private boolean useGenotypeCalls = true;
    private boolean useReferenceProbability = false;
    private boolean useAlleleProbabilities = false;
    
    public CompressedMLMusingDoubleMatrix(MLMPlugin parentPlugin, Datum dataset, DistanceMatrix kinshipMatrix, boolean useCompression, boolean useP3D, double compression) {
        this.parentPlugin = parentPlugin;
        this.kinshipMatrix = kinshipMatrix;
        this.useCompression = useCompression;
        this.useP3D = useP3D;
        this.compression = compression;
        this.myWeightMatrix = null;
        datasetName = dataset.getName();
        
        if (dataset.getData().getClass().equals(GenotypePhenotype.class)) {
            myGenoPheno = (GenotypePhenotype) dataset.getData();
            myPhenotype = myGenoPheno.phenotype();
            myGenotype = myGenoPheno.genotypeTable();
            hasGenotype = true;
        } else if (dataset.getData() instanceof Phenotype) {
            myGenoPheno = null;
            myPhenotype = (Phenotype) dataset.getData();
            myGenotype = null;
            hasGenotype = false;
        } else {
            myGenoPheno = null;
            myPhenotype = null;
            myGenotype = null;
            hasGenotype = false;
        }
        
        // String[] headerMain = new String[]{"Trait", "Marker", "Locus", "Site", "df", "F", "p", "errordf", "markerR2", "Genetic Var", "Residual Var", "-2LnLikelihood"};
        String[] headerMain = new String[]{AssociationConstants.STATS_HEADER_TRAIT, AssociationConstants.STATS_HEADER_MARKER,
            AssociationConstants.STATS_HEADER_CHR, AssociationConstants.STATS_HEADER_POSITION,
            "df", "F", AssociationConstants.STATS_HEADER_P_VALUE,
            "add_effect", "add_F", "add_p", "dom_effect", "dom_F", "dom_p", "errordf",
            "MarkerR2", "Genetic Var", "Residual Var", "-2LnLikelihood"};
        String[] headerAlleles = new String[]{"Trait", "Marker", "Locus", "Site", "Allele", "Effect", "Obs"};
        String[] headerCompression = new String[]{"Trait", "# groups", "Compression", "-2LnLk", "Var_genetic", "Var_error"};
        
        if (parentPlugin.isWriteOutputToFile()) {
        	String outputbase = parentPlugin.getOutputName();
        	String datasetNameNoSpace = datasetName.trim().replaceAll("\\ ", "_");
        	
        	StringBuilder sb = new StringBuilder();
        	sb.append(outputbase).append("_").append(datasetNameNoSpace).append("_stats.txt");
            siteReportBuilder = TableReportBuilder.getInstance("Marker Statistics - " + datasetName, headerMain, sb.toString());
            
        	sb = new StringBuilder();
        	sb.append(outputbase).append("_").append(datasetNameNoSpace).append("_effects.txt");
            alleleReportBuilder = TableReportBuilder.getInstance("Allele Estimates - " + datasetName, headerAlleles, sb.toString());
            
        	sb = new StringBuilder();
        	sb.append(outputbase).append("_").append(datasetNameNoSpace).append("_compression.txt");
            if (useCompression) compressionReportBuilder = TableReportBuilder.getInstance("Compression - " + datasetName, headerCompression, sb.toString());
            else compressionReportBuilder = null;
        } else {
            siteReportBuilder = TableReportBuilder.getInstance("Marker Statistics - " + datasetName, headerMain);
            alleleReportBuilder = TableReportBuilder.getInstance("Allele Estimates - " + datasetName, headerAlleles);
            if (useCompression) compressionReportBuilder = TableReportBuilder.getInstance("Compression - " + datasetName, headerCompression);
            else compressionReportBuilder = null;
        }

//        solve();
    }

    public CompressedMLMusingDoubleMatrix(WeightedMLMPlugin parentPlugin, Datum dataset, DistanceMatrix kinshipMatrix, Datum weights, boolean useCompression, boolean useP3D, double compression) {
       
        this.parentPlugin = parentPlugin;
        //this.parentPlugin = null;
        this.kinshipMatrix = kinshipMatrix;
        this.useCompression = useCompression;
        this.useP3D = useP3D;
        this.compression = compression;
        this.myWeightMatrix = (Phenotype)weights.getData();
        datasetName = dataset.getName();
        
        if (dataset.getData().getClass().equals(GenotypePhenotype.class)) {
            myGenoPheno = (GenotypePhenotype) dataset.getData();
            myPhenotype = myGenoPheno.phenotype();
            myGenotype = myGenoPheno.genotypeTable();
            hasGenotype = true;
        } else if (dataset.getData() instanceof Phenotype) {
            myGenoPheno = null;
            myPhenotype = (Phenotype) dataset.getData();
            myGenotype = null;
            hasGenotype = false;
        } else {
            myGenoPheno = null;
            myPhenotype = null;
            myGenotype = null;
            hasGenotype = false;
        }
        
        // String[] headerMain = new String[]{"Trait", "Marker", "Locus", "Site", "df", "F", "p", "errordf", "markerR2", "Genetic Var", "Residual Var", "-2LnLikelihood"};
        String[] headerMain = new String[]{AssociationConstants.STATS_HEADER_TRAIT, AssociationConstants.STATS_HEADER_MARKER,
            AssociationConstants.STATS_HEADER_CHR, AssociationConstants.STATS_HEADER_POSITION,
            "df", "F", AssociationConstants.STATS_HEADER_P_VALUE,
            "add_effect", "add_F", "add_p", "dom_effect", "dom_F", "dom_p", "errordf",
            "MarkerR2", "Genetic Var", "Residual Var", "-2LnLikelihood"};
        String[] headerAlleles = new String[]{"Trait", "Marker", "Locus", "Site", "Allele", "Effect", "Obs"};
        String[] headerCompression = new String[]{"Trait", "# groups", "Compression", "-2LnLk", "Var_genetic", "Var_error"};
        
        if (parentPlugin.isWriteOutputToFile()) {
                String outputbase = parentPlugin.getOutputName();
                String datasetNameNoSpace = datasetName.trim().replaceAll("\\ ", "_");
                
                StringBuilder sb = new StringBuilder();
                sb.append(outputbase).append("_").append(datasetNameNoSpace).append("_stats.txt");
            siteReportBuilder = TableReportBuilder.getInstance("Marker Statistics - " + datasetName, headerMain, sb.toString());
            
                sb = new StringBuilder();
                sb.append(outputbase).append("_").append(datasetNameNoSpace).append("_effects.txt");
            alleleReportBuilder = TableReportBuilder.getInstance("Allele Estimates - " + datasetName, headerAlleles, sb.toString());
            
                sb = new StringBuilder();
                sb.append(outputbase).append("_").append(datasetNameNoSpace).append("_compression.txt");
            if (useCompression) compressionReportBuilder = TableReportBuilder.getInstance("Compression - " + datasetName, headerCompression, sb.toString());
            else compressionReportBuilder = null;
        } else {
            siteReportBuilder = TableReportBuilder.getInstance("Marker Statistics - " + datasetName, headerMain);
            alleleReportBuilder = TableReportBuilder.getInstance("Allele Estimates - " + datasetName, headerAlleles);
            if (useCompression) compressionReportBuilder = TableReportBuilder.getInstance("Compression - " + datasetName, headerCompression);
            else compressionReportBuilder = null;
        }

//        solve();
    }
    public void useGenotypeCalls(boolean use) {
    	useGenotypeCalls = use;
    }

    public void useReferenceProbability(boolean use) {
    	useReferenceProbability = use;
    }

    public void useAlleleProbabilities(boolean use) {
    	useAlleleProbabilities = use;
    }
    
    public List<Datum> solve() {
    	List<Datum> results = new LinkedList<Datum>();

    	int numberOfMarkers = 0;
    	if (hasGenotype) numberOfMarkers = myGenotype.numberOfSites();
    	List<PhenotypeAttribute> dataAttributeList =  myPhenotype.attributeListOfType(ATTRIBUTE_TYPE.data);
    	factorAttributeList =  myPhenotype.attributeListOfType(ATTRIBUTE_TYPE.factor);
    	covariateAttributeList =  myPhenotype.attributeListOfType(ATTRIBUTE_TYPE.covariate);
    	TaxaAttribute myTaxaAttribute = myPhenotype.taxaAttribute();
    	
        int numberOfPhenotypes = dataAttributeList.size();

        //calculate total iterations
        int expectedIterations = numberOfPhenotypes * numberOfMarkers;
        int iterationsSofar = 0;

        //cycle through the phenotypes
        for (PhenotypeAttribute attr : dataAttributeList) {
            //get phenotype data
        	myLogger.debug("Running MLM for " + attr.name());
            double[] phenotypeData = doubleDataFromAttribute(attr);

            //get the taxa
            Taxon[] theTaxa = myTaxaAttribute.allTaxa();

            //keep track of missing rows
            OpenBitSet missing = new OpenBitSet(attr.missing());
            for (PhenotypeAttribute factorAttribute : factorAttributeList) missing.or(factorAttribute.missing());
            for (PhenotypeAttribute covariateAttribute : covariateAttributeList) missing.or(covariateAttribute.missing());

            //update missing for taxa not in the kinship matrix or the distance matrix.
            //Create kinship and distance matrices with taxa in phenotype
            TaxaList nonmissingIds = updateMissingWithKinship(missing, theTaxa);
            DistanceMatrix kin = new DistanceMatrix(kinshipMatrix, nonmissingIds);

            //calculate the number of nonmissing observations
            int totalObs = attr.size();
            int nonMissingObs = totalObs - (int) missing.cardinality();

            //create phenotype matrix
            double[] nonMissingData = AssociationUtils.getNonMissingDoubles(phenotypeData, missing);
            Taxon[] nonMissingTaxa = AssociationUtils.getNonMissingValues(theTaxa, missing);
            DoubleMatrix y = DoubleMatrixFactory.DEFAULT.make(nonMissingObs,1, nonMissingData);

            //create the Z matrix
            DoubleMatrix Z = DoubleMatrixFactory.DEFAULT.make(nonMissingObs, kin.numberOfTaxa());
            for (int i = 0; i < nonMissingObs; i++) {
                   Z.set(i, kin.whichIdNumber(nonMissingTaxa[i]), 1);
            }

            //fixed effects matrix
            DoubleMatrix fixed = AssociationUtils.createFixedEffectsArray(factorAttributeList, covariateAttributeList, missing, nonMissingObs);
            DoubleMatrix W = null;
            DoubleMatrix WOriginal = null;
            EMMAforDoubleMatrix emlm = null;
            DoubleMatrix[] zk = null;
            
            //Check to see if weightedMLM is used
            if(myWeightMatrix!=null) {
                //Grab Attribute from myWeightMatrix which matches attr.name()
                int weightAttrIndex = myWeightMatrix.attributeIndexForName(attr.name());
                //Check to see if Attribute index = -1(missing)
                if(weightAttrIndex!=-1) {
                    //Grab PhenotypeAttribute Object
                    PhenotypeAttribute weightAttribute = myWeightMatrix.attribute(weightAttrIndex);
                    //Remove Missing values from Weight List
                    double[] weightValues = doubleDataFromAttribute(weightAttribute);
                    double[] nonMissingWeights = AssociationUtils.getNonMissingDoubles(weightValues, missing);
                    
                    //Check to make sure numRows(y) == numRows(weight)
                    if(nonMissingObs == nonMissingWeights.length) { 
                        //Set Original Weight for use in Testing
                        WOriginal = DoubleMatrixFactory.DEFAULT.diagonal(nonMissingWeights);
                        
                        //Calculate W = W^{-1/2}
                        for(int i = 0;i<nonMissingWeights.length;i++) {
                            nonMissingWeights[i] = Math.pow(nonMissingWeights[i],-.5);
                        }
                              
                        //Build Weight Matrix W
                        W = DoubleMatrixFactory.DEFAULT.diagonal(nonMissingWeights);
                        
                        //run emlm with weighted values
                        zk = computeZKZ(W.mult(y), W.mult(fixed), W.mult(Z), kin, attr.name());
                        emlm = new EMMAforDoubleMatrix(W.mult(y), W.mult(fixed), zk[1], zk[0], 0, Double.NaN);
                        emlm.solve();
                        //recompute zk with original values for testing
                        zk = computeZKZ(y, fixed, Z, kin, attr.name());
                       
                    }
                    else {
                       //run Standard EMMA without Weights
                        zk = computeZKZ(y, fixed, Z, kin, attr.name());
                        emlm = new EMMAforDoubleMatrix(y, fixed, zk[1], zk[0], 0, Double.NaN);
                        emlm.solve();
                    }
                }
                else {
                    //Maybe send message to user saying We are using normal MLM for this attribute
                    //run Standard EMMA without Weights
                    zk = computeZKZ(y, fixed, Z, kin, attr.name());
                    emlm = new EMMAforDoubleMatrix(y, fixed, zk[1], zk[0], 0, Double.NaN);
                    emlm.solve();
                }
            }
            else {
                //run Standard EMMA without Weights
                zk = computeZKZ(y, fixed, Z, kin, attr.name());
                emlm = new EMMAforDoubleMatrix(y, fixed, zk[1], zk[0], 0, Double.NaN);
                emlm.solve();
            }
            
            //fit data without markers
            genvar = emlm.getVarRan();
            resvar = emlm.getVarRes();
            lnlk = emlm.getLnLikelihood();
            int baseModeldf = emlm.getDfModel();

            //record the results
            if (outputResiduals) {
                List<Taxon> phenotypeTaxa = Arrays.asList(nonMissingTaxa);
                Datum residuals = createResPhenotype(emlm, phenotypeTaxa, attr.name());
                results.add(residuals);
                if(parentPlugin.isWriteOutputToFile()){
                    ExportPlugin exporter = new ExportPlugin(null, false);
                    String outfile = parentPlugin.getOutputName() + "_" + residuals.getName() +  "_residuals.txt";
                    exporter.setSaveFile(outfile);
                    exporter.performFunction(new DataSet(residuals, parentPlugin));
                }
            }
            
            Object[] tableRow;
            //{"Trait", "Marker", "Locus", "Site", "df", "F", "p", "errordf", "markerR2", "Genetic Var", "Residual Var", "-2LnLikelihood"}
            //{"Trait","Marker","Chr","Pos","Locus","Site","df","F","p","add_effect","add_F","add_p","dom_effect","dom_F","dom_p","errordf","MarkerR2","Genetic Var","Residual Var", "-2LnLikelihood"}
            tableRow = new Object[]{
            		attr.name(),
            		"None",
            		"",
            		"",
            		new Integer(0),
            		new Double(Double.NaN),
            		new Double(Double.NaN),
            		new Double(Double.NaN),
            		new Double(Double.NaN),
            		new Double(Double.NaN),
            		new Double(Double.NaN),
            		new Double(Double.NaN),
            		new Double(Double.NaN),
            		new Integer(nonMissingObs - baseModeldf),
            		new Double(Double.NaN),
            		new Double(genvar),
            		new Double(resvar),
            		new Double(-2 * lnlk)};

            siteReportBuilder.add(tableRow);

            //the BLUPs
            //not implemented

            if (useP3D) {
                DoubleMatrix ZKZ = zk[0].mult(zk[1]).tcrossproduct(zk[0]);
                //Check to see if W exists. Alter V calculation accordingly
            	if(W==null) {
            	    Vminus = new SymmetricMatrixInverterDM(calculateV(ZKZ, genvar, resvar));
            	}
            	else {
            	    Vminus = new SymmetricMatrixInverterDM(calculateV(ZKZ, WOriginal, genvar, resvar));
            	}
            }

            //iterate markers
            if (testMarkers) {
                for (int m = 0; m < numberOfMarkers; m++) {
                	OpenBitSet missingObsForSite = new OpenBitSet(missing);
                	missingObsForSite.or(missingForSite(m));
                	
                    //only data for which missing=false are in the Z matrix
                    //the block below finds the rows of Z that have no marker data.
                    //Those rows/columns will need to be removed from ZKZ or from V, depending on the analysis method.
                	OpenBitSet missingFromZ = new OpenBitSet(nonMissingObs);
                	
                    int nonMissingCount = 0;
                    for (int i = 0; i < totalObs; i++) {
                        if (!missing.fastGet(i)) {
                            if (missingObsForSite.fastGet(i)) {
                            	missingFromZ.fastSet(nonMissingCount);
                            }
                            nonMissingCount++;
                        }
                    }

                    //test for all missing
                    if (nonMissingCount == 0) {
                        myLogger.error("All values missing for " + myGenotype.siteName(m));
                        continue;
                    }

                    //adjust y for missing data
                    DoubleMatrix ymarker = AssociationUtils.getNonMissingValues(y, missingFromZ);

                    //adjust the fixed effects
                    DoubleMatrix fixed2 = AssociationUtils.getNonMissingValues(fixed, missingFromZ);

                    //add marker data to fixed effects
                    ArrayList<Byte> markerIds = new ArrayList<>();
                    int nAlleles = 0;
                    int markerdf = 0;
                    DoubleMatrix X;
                    int[] alleleCounts = null;
                    
                    if (useGenotypeCalls) {
                    	byte[] genotypes = ModelEffectUtils.genotypesToUnphasedSorted(AssociationUtils.getNonMissingBytes(myGenoPheno.genotypeAllTaxa(m), missingObsForSite));
                        FactorModelEffect markerEffect = new FactorModelEffect(ModelEffectUtils.getIntegerLevels(genotypes, markerIds), true);
                        X = fixed2.concatenate(markerEffect.getX(), false);
                        nAlleles = markerEffect.getNumberOfLevels();
                        alleleCounts = markerEffect.getLevelCounts();
                        markerdf = nAlleles - 1;
                    } else if (useReferenceProbability) { 
                        double[] genotypes = AssociationUtils.getNonMissingDoubles(myGenoPheno.referenceProb(m), missingObsForSite);
                        int nrows = genotypes.length;
                        X = fixed2.concatenate(DoubleMatrixFactory.DEFAULT.make(nrows, 1, genotypes), false);
                        nAlleles = 1;
                        alleleCounts = new int[]{nrows};
                        markerdf = 1;
                    } else {
                    	X = null;
                    }
                    
                    CompressedMLMResult result = new CompressedMLMResult();
                    //need to add marker information to result once Alignment is stable

                    if (useP3D) {
                        testMarkerUsingP3D(result, ymarker, X, Vminus.getInverse(missingFromZ, nonMissingObs), markerdf, markerIds);
                    } else {
                    	DoubleMatrix Zsel = AssociationUtils.getNonMissingValues(zk[0], missingFromZ);
                        if (!testMarkerUsingEMMA(result, ymarker, X, zk[1], Zsel, nAlleles, markerIds)) continue;
                        markerdf = result.modeldf - baseModeldf;
                    }

                    //if the results are to be filtered on pmax check for that condition
                    boolean recordTheseResults = true;
                    if (parentPlugin.isFilterOutput() && result.p > parentPlugin.getMaxp()) {
                        recordTheseResults = false;
                    }

                    if (recordTheseResults) {
                        //add result to main
                        //{"Trait","Marker","Chr","Pos","Locus","Site","df","F","p","errordf","MarkerR2","Genetic Var","Residual Var", "-2LnLikelihood"};
                    	//results with additive and dominance effects
                    	//{"Trait","Marker","Chr","Pos","Locus","Site","df","F","p","add_effect","add_F","add_p","dom_effect","dom_F","dom_p","errordf","MarkerR2","Genetic Var","Residual Var", "-2LnLikelihood"}
                    	
                        String markername = myGenotype.siteName(m);
                        String chr = "";
                        String pos = "";
                        String locus = myGenotype.chromosomeName(m);
                        String site = Integer.toString(myGenotype.chromosomalPosition(m));
                        double errordf = (double) (ymarker.numberOfRows() - result.modeldf);

                        tableRow = new Object[]{attr.name(),
                        		markername,
                        		locus,
                        		site,
                        		new Integer(markerdf),
                        		new Double(result.F),
                        		new Double(result.p),
                        		new Double(result.addEffect),
                        		new Double(result.Fadd),
                        		new Double(result.padd),
                        		new Double(result.domEffect),
                        		new Double(result.Fdom),
                        		new Double(result.pdom),
                        		new Double(errordf),
                        		new Double(result.r2),
                        		new Double(genvar),
                        		new Double(resvar),
                        		new Double(-2 * lnlk)};
                        siteReportBuilder.add(tableRow);

                        //add result to alleles
                        //"Trait","Marker","Chr","Pos","Allele","Effect", obs
                        int numberOfRowsKept = totalObs - (int) missingObsForSite.cardinality();
                        if (useReferenceProbability) {
                        	tableRow = new Object[]{attr.name(),
                        			markername,
                        			locus,
                        			site,
                        			"",
                        			result.beta.get(result.beta.numberOfRows() - 1, 0),
                        			numberOfRowsKept
                        	};

                            //record the results
                        	alleleReportBuilder.add(tableRow);
                        } else if (nAlleles > 1) {
                            for (int a = 0; a < nAlleles; a++) {
                                Double estimate;
                                if (a < nAlleles - 1) {
                                    estimate = result.beta.get(result.beta.numberOfRows() - nAlleles + 1 + a, 0);
                                } else {
                                    estimate = 0.0;
                                }
                                tableRow = new Object[]{attr.name(),
                                		markername,
                                		locus,
                                		site,
//                                		markerIds.get(a),
                                		NucleotideAlignmentConstants.getNucleotideIUPAC(markerIds.get(a)),
                                		estimate,
                                		alleleCounts[a]
                                };

                                //record the results
                            	alleleReportBuilder.add(tableRow);
                            }
                        }

                    }
                    iterationsSofar++;
                    int progress = (int) ((double) iterationsSofar / (double) expectedIterations * 100);
                    progress = Math.min(99, progress);
                    parentPlugin.updateProgress(progress);
                }
            }

        }

        parentPlugin.updateProgress(100);

        results.addAll(formatResults());
        
        return results;
    }

//    private BitSet missingForSiteX(int site) {
//    	int ntaxa = myGenotype.numberOfTaxa();
//    	OpenBitSet missing = new OpenBitSet(ntaxa);
//    	if (useGenotypeCalls) {
//        	byte[] siteGenotype = myGenotype.genotypeAllTaxa(site);
//        	byte missingByte = GenotypeTable.UNKNOWN_DIPLOID_ALLELE;
//        	for (int i = 0; i < ntaxa; i++) if (siteGenotype[i] == missingByte) missing.fastSet(i);
//    	} else if (useReferenceProbability) {
//    		for (int t = 0; t < ntaxa; t++) {
//    			if (myGenotype.referenceProbability(t, site) == Float.NaN) missing.fastSet(t);
//    		}
//    	} else {
//    		for (int t = 0; t < ntaxa; t++) {
//    			if (myGenotype.alleleProbability(t, site, SITE_SCORE_TYPE.DepthA) == Float.NaN) missing.fastSet(t);
//    		}
//    	}
//    	return missing;
//    }
    
    private BitSet missingForSite(int site) {
    	//returns BitSet with missing set for each observation with a missing genotype value
    	int nobs = myGenoPheno.phenotype().numberOfObservations();
    	OpenBitSet missing = new OpenBitSet(nobs);
    	if (useGenotypeCalls) {
        	byte[] siteGenotype = myGenoPheno.genotypeAllTaxa(site);
        	byte missingByte = GenotypeTable.UNKNOWN_DIPLOID_ALLELE;
        	for (int i = 0; i < nobs; i++) {
        		if (siteGenotype[i] == missingByte) missing.fastSet(i);
        	}
    	} else if (useReferenceProbability) {
    		float[] probs = myGenoPheno.referenceProb(site);
    		for (int t = 0; t < nobs; t++) {
    			if (probs[t] == Float.NaN) missing.fastSet(t);
    		}
    	} else {
			float[] probs = myGenoPheno.alleleProbsOfType(SITE_SCORE_TYPE.DepthA, site);
    		for (int t = 0; t < nobs; t++) {
    			if (probs[t] == Float.NaN) missing.fastSet(t);
    		}
    	}
    	return missing;
    }
    
    private double[] doubleDataFromAttribute(PhenotypeAttribute attribute) {
    	float[] floatData = (float[]) ((NumericAttribute) attribute).allValues();
    	int n = floatData.length;
    	double[] doubleData = new double[n];
    	for (int i = 0; i < n; i++) doubleData[i] = floatData[i];
    	return doubleData;
    }
    
    private String getTabbedStringFromArray(Object[] array) {
        StringBuffer sb = new StringBuffer();
        sb.append(array[0]);
        int n = array.length;
        for (int i = 1; i < n; i++) {
            sb.append("\t").append(array[i]);
        }
        return sb.toString();
    }

    public List<Datum> formatResults() {
        LinkedList<Datum> output = new LinkedList<Datum>();

        //generate comments
        StringBuilder options = new StringBuilder();
        options.append("Use compression = ").append(useCompression).append("\n");
        options.append("Use P3D = ").append(useP3D).append("\n");
        if (useCompression) {
            options.append(", compression level = ").append(compression).append("\n");
        }
        if (useP3D) {
            options.append("P3D = ").append(useP3D).append(". Variance components were estimated only for the model without any markers.\n");
        } else {
            options.append("P3D = ").append(useP3D).append(". Variance components were estimated for each marker.\n");
        }

        StringBuilder model = new StringBuilder();
        model.append("Model: trait = mean + ");
        int nFactors = factorAttributeList.size();
        for (PhenotypeAttribute factor:factorAttributeList) {
            model.append(factor.name()).append(" + ");
        }
        int nCovar = covariateAttributeList.size();
        for (PhenotypeAttribute cov : covariateAttributeList) {
            model.append(cov.name()).append(" + ");
        }
        model.append("marker\n");

        String reportName = "MLM_statistics_for_" + datasetName;
        StringBuilder comment = new StringBuilder();
        comment.append("MLM statistics for compressed MLM\n");
        comment.append("Dataset: ").append(datasetName).append("\n");
        comment.append(options).append(model);
        TableReport myTableReport = siteReportBuilder.build();
        if (myTableReport != null) output.add(new Datum(reportName, myTableReport, comment.toString()));

        reportName = "MLM_effects_for_" + datasetName;
        comment = new StringBuilder();
        comment.append("MLM SNP effect estimates\n");
        comment.append("Dataset: ").append(datasetName).append("\n");
        comment.append(options).append(model);
        myTableReport = alleleReportBuilder.build();
        if (myTableReport != null) output.add(new Datum(reportName, myTableReport, comment.toString()));
        
        if (useCompression) {
        	reportName = "MLM_compression_for_" + datasetName;
        	comment = new StringBuilder();
        	comment.append("MLM compression report\n");
        	comment.append("Dataset: ").append(datasetName).append("\n");
        	comment.append(options).append(model);
        	myTableReport = compressionReportBuilder.build();
            if (myTableReport != null) output.add(new Datum(reportName, myTableReport, comment.toString()));        }
        
        return output;
    }

    /**
     * Computes ZKZ. If compression is specified then the compressed ZKZ is calculated along with compressed versions of Z and K.
     * @param data	the phenotype data. Needed for optimizing compression.
     * @param X	the incidence matrix specifying all fixed effects other than markers
     * @param Z the kinship incidence matrix
     * @param kin	the genetic similarity matrix
     * @param traitname	the name of the phenotype passed in data
     * @return	an array containing the Z matrix as its first element and the K matrix as its second element. If compression is specified, then both are the compressed versions.
     */
    public DoubleMatrix[] computeZKZ(DoubleMatrix data, DoubleMatrix X, DoubleMatrix Z, DistanceMatrix kin, String traitname) {
    	myLogger.debug("Running compression for " + traitname);
    	DoubleMatrix[] zkMatrices = new DoubleMatrix[2]; 
    	CompressedDoubleMatrix.kinshipMethod kinmethod = CompressedDoubleMatrix.kinshipMethod.avg;
    	
        //Kmatrix
        int nkin = kin.getSize();
        int nrow = nkin;
        int ncol = nrow;

        DoubleMatrix K = DoubleMatrixFactory.DEFAULT.make(nrow, ncol);
        for (int r = 0; r < nrow; r++) {
            for (int c = 0; c < ncol; c++) {
                K.set(r, c, kin.getDistance(r, c));
            }
        }

        if (!useCompression) {
        	zkMatrices[0] = Z;
        	zkMatrices[1] = K;
        } else if (Double.isNaN(compression)) {
            //are taxa replicated?
            //sum columns of Z. If any sum > 1, then yes
            int n = Z.numberOfColumns();
            int count = 0;
            boolean taxaReplicated = false;
            while (count < n && !taxaReplicated) {
                if (Z.columnSum(count++) > 1.5) {
                    taxaReplicated = true;
                }
            }

            DistanceMatrix distance = calculateDistanceFromKin(kin);
            CompressedDoubleMatrix cm = new CompressedDoubleMatrix(kin, new UPGMATree(distance));
            EMMAforDoubleMatrix emlm = new EMMAforDoubleMatrix(data, X, K, Z, 0, Double.NaN);
            
            emlm.solve();
            double bestlnlk = emlm.getLnLikelihood();
            int bestCompression = nkin;

            double exponent = 1;
            double base = 0.98;
            double maxexponent = Math.log(1 / ((double) nkin)) / Math.log(base);
            parentPlugin.updateProgress((int) (exponent * 100 / maxexponent));
            //int g = (int) (nkin * Math.pow(base, exponent));
            int g = (int) (nkin);
            while (g > 3) {
                cm.setNumberOfGroups(g);
                DoubleMatrix compressedZ = cm.getCompressedZ(Z);
                DoubleMatrix compressedK = cm.getCompressedMatrix(kinmethod);
                try {
                    emlm = new EMMAforDoubleMatrix(data, X, compressedK, compressedZ, 0, Double.NaN);
                    emlm.solve();

                    //output number of groups, compression level (= number of taxa / number of groups), -2L, genvar, resvar
                    compressionReportBuilder.add(new Object[]{traitname, g,
                                ((double) nkin) / ((double) g),
                                -2 * emlm.getLnLikelihood(),
                                emlm.getVarRan(),
                                emlm.getVarRes()});

                    if (Double.isNaN(bestlnlk) || emlm.getLnLikelihood() > bestlnlk) {
                        bestlnlk = emlm.getLnLikelihood();
                        bestCompression = g;
                        resvar = emlm.getVarRes();
                        genvar = emlm.getVarRan();
                    }
                } catch (Exception e) {
                    System.out.println("Compression failed for g = " + g);
                }

                int prev = g;
                while (g == prev) {
                    exponent++;
                    int prog = (int) (exponent * 100 / maxexponent);
                    prog = Math.min(prog, 99);
                    parentPlugin.updateProgress(prog);
                    g = (int) (nkin * Math.pow(base, exponent));
                }
            }

            //for g = 1 use GLM to estimate beta and errvar
//            if (!taxaReplicated) {
//                SweepFast sweep = new SweepFast(X, data);
//                sweep.XTXSweepSetDmin();
//                n = X.numberOfColumns();
//                double ssres = sweep.getResidualSS();
//                double errordf = (double) (data.numberOfRows() - n);
//                double errvar = ssres / errordf;
//                double lnlk = (errordf * Math.log(2 * Math.PI * errvar) + errordf);
//
//                compressionReportBuilder.add(new Object[]{traitname, g,
//                            ((double) nkin) / ((double) g),
//                            lnlk,
//                            new Double(0.0),
//                            errvar});
//
//                if (Double.isNaN(bestlnlk) || emlm.getLnLikelihood() > bestlnlk) {
//                    bestlnlk = emlm.getLnLikelihood();
//                    bestCompression = g;
//                    resvar = emlm.getVarRes();
//                    genvar = 0;
//                }
//
//            }

            cm.setNumberOfGroups(bestCompression);
            zkMatrices[0] = cm.getCompressedZ(Z);
            zkMatrices[1] = cm.getCompressedMatrix(kinmethod);
            parentPlugin.updateProgress(0);

        } else {
            DistanceMatrix distance = calculateDistanceFromKin(kin);
            CompressedDoubleMatrix cm = new CompressedDoubleMatrix(kin, new UPGMATree(distance));
            int g = (int) Math.round(nkin / compression);
            cm.setNumberOfGroups(g);
            zkMatrices[0] = cm.getCompressedZ(Z);
            zkMatrices[1] = cm.getCompressedMatrix(kinmethod);
        }
        
        return zkMatrices;
    }

    public boolean testMarkerUsingEMMA(CompressedMLMResult result, DoubleMatrix y, DoubleMatrix X, DoubleMatrix K, DoubleMatrix Z, int nAlleles, ArrayList<Byte> markerIds) {
        EMMAforDoubleMatrix emlm;
        try {
            emlm = new EMMAforDoubleMatrix(y, X, K, Z, nAlleles, Double.NaN);
        } catch(Exception except) {
            myLogger.error(except.getMessage());
            except.printStackTrace();
            return false;
        }

        emlm.solve();
        result.beta = emlm.getBeta();
        double[] Fp = emlm.getMarkerFp();
        result.F = Fp[0];
        result.p = Fp[1];
        result.modeldf = emlm.getDfModel();
        genvar = emlm.getVarRan();
        resvar = emlm.getVarRes();
        lnlk = emlm.getLnLikelihood();
        
        calculateRsquare(X, y, emlm.getInvH(), result, nAlleles - 1);
        
        boolean markerTest = markerIds.size() == 3;
        if (markerTest) {
        	markerTest = markerTest && !GenotypeTableUtils.isHeterozygous(markerIds.get(0));
        	markerTest = markerTest && !GenotypeTableUtils.isHeterozygous(markerIds.get(1));
        	markerTest = markerTest && GenotypeTableUtils.isHeterozygous(markerIds.get(2));
        }
        if (markerTest && Fp.length == 8) { //calculate additive and dominance tests and effects
        	//from EMMA,  return new double[]{F,p,addEffect,Fadd,padd,domEffect,Fdom,pdom}
        	result.addEffect = Fp[2];
        	result.Fadd = Fp[3];
        	result.padd = Fp[4];
        	result.domEffect = Fp[5];
        	result.Fdom = Fp[6];
        	result.pdom = Fp[7];
        }
        return true;
    }

    public void testMarkerUsingP3D(CompressedMLMResult result, DoubleMatrix y, DoubleMatrix X, DoubleMatrix invV, int markerdf, ArrayList<Byte> markerIds) {
        //calculate beta
        DoubleMatrix invXVX = X.crossproduct(invV).mult(X);
        invXVX.invert();
        result.beta = invXVX.mult(X.crossproduct(invV.mult(y)));

        //test for markerdf = 0
        if (markerdf == 0) {
            result.F = Double.NaN;
            result.p = Double.NaN;
            result.r2 = 0.0;
        } else {  //full model
        	
            //calculate F test, p-value of F test
            int nparm = result.beta.numberOfRows();
            DoubleMatrix M = DoubleMatrixFactory.DEFAULT.make(markerdf, nparm, 0);
            for (int i = 0; i < markerdf; i++) {
                M.set(i, nparm - markerdf + i, 1);
            }
            DoubleMatrix Mb = M.mult(result.beta);
            DoubleMatrix invMiM = M.mult(invXVX.tcrossproduct(M));
            try {
                invMiM.invert();
                result.F = Mb.crossproduct(invMiM.mult(Mb)).get(0, 0) / markerdf;
            } catch (Exception ex) {
                result.F = Double.NaN;
            }
            try {
                result.p = LinearModelUtils.Ftest(result.F, markerdf, y.numberOfRows() - nparm);
            } catch (Exception e) {
                result.p = Double.NaN;
            }

            calculateRsquare(X, y, invV, result, markerdf);
            
            boolean markerTest = markerIds.size() == 3;
            if (markerTest) {
            	markerTest = markerTest && !GenotypeTableUtils.isHeterozygous(markerIds.get(0));
            	markerTest = markerTest && !GenotypeTableUtils.isHeterozygous(markerIds.get(1));
            	markerTest = markerTest && GenotypeTableUtils.isHeterozygous(markerIds.get(2));
            }
            if (markerdf == 2 && markerTest) { //calculate additive and dominance tests and effects
            	
            	//additive test 
                M = DoubleMatrixFactory.DEFAULT.make(1, nparm, 0);
                M.set(0, nparm - 2, 0.5);
                M.set(0, nparm - 1, -0.5);
                    
                Mb = M.mult(result.beta);
                result.addEffect = Mb.get(0, 0);
                try {
                   result.Fadd = Mb.get(0, 0) * Mb.get(0,0) / (M.mult(invXVX.tcrossproduct(M))).get(0,0);
                } catch (Exception ex) {
                    result.Fadd = Double.NaN;
                }
                try {
                    result.padd = LinearModelUtils.Ftest(result.Fadd, 1, y.numberOfRows() - nparm);
                } catch (Exception e) {
                    result.padd = Double.NaN;
                }

                //dominance test
                M = DoubleMatrixFactory.DEFAULT.make(1, nparm, 0);
                M.set(0, nparm - 2, -0.5);
                M.set(0, nparm - 1, -0.5);
                    
                Mb = M.mult(result.beta);
                result.domEffect = Mb.get(0, 0);
                try {
                    result.Fdom = Mb.get(0, 0) * Mb.get(0,0) / (M.mult(invXVX.tcrossproduct(M))).get(0,0);
                } catch (Exception ex) {
                    result.Fdom = Double.NaN;
                }
                try {
                    result.pdom = LinearModelUtils.Ftest(result.Fdom, 1, y.numberOfRows() - nparm);
                } catch (Exception e) {
                    result.pdom = Double.NaN;
                }
                
            }
        }

    }
    
    private void calculateRsquare(DoubleMatrix X, DoubleMatrix y, DoubleMatrix invV, CompressedMLMResult result, int markerdf) {
        //calculate R2
        //from Buse(1973) Am. Stat. 27:106-108.
        //R^2 = ymarker'*inverseV*ymarker / (y-mean)'*inverseV*(y-mean)
        //where ymarker = yhat(full model) - yhat(model without marker)
        //as Xm*betam where Xm is the columns of X due to the marker and betam is the portion of beta due to the markers adjusted for the marker mean

        int dimX = X.numberOfColumns();
        int dimXreduced = dimX - markerdf;
        int[] colsToKeep = new int[dimXreduced];
        for (int i = 0; i < dimXreduced; i++) {
            colsToKeep[i] = i;
        }
        DoubleMatrix Xreduced = X.getSelection(null, colsToKeep);

        //calculate reduced beta
        DoubleMatrix invXVX = Xreduced.crossproduct(invV).mult(Xreduced);
        invXVX.invert();
        DoubleMatrix betaReduced = invXVX.mult(Xreduced.crossproduct(invV.mult(y)));

        //calculate yhat = yhatFull - yhatReduced
        DoubleMatrix yhat = X.mult(result.beta);
        DoubleMatrix yhatReduced = Xreduced.mult(betaReduced);
        yhat.minusEquals(yhatReduced);

        //calculate ydev = y - mean
        double sum = 0;
        int n = y.numberOfRows();
        for (int i = 0; i < n; i++) {
            sum += y.get(i, 0);
        }
        double mean = sum / n;

        DoubleMatrix ydev = y.scalarAdd(-mean);
        double numerator = yhat.crossproduct(invV).mult(yhat).get(0, 0);
        double denominator = ydev.crossproduct(invV).mult(ydev).get(0, 0);
        result.r2 = numerator / denominator;
    }

    public DoubleMatrix calculateV(DoubleMatrix ZKZ, double genvar, double resvar) {
        DoubleMatrix V = ZKZ.scalarMult(genvar);
        int n = V.numberOfRows();
        for (int i = 0; i < n; i++) {
            V.set(i, i, V.get(i, i) + resvar);
        }
        return V;
    }
    
    public DoubleMatrix calculateV(DoubleMatrix ZKZ, DoubleMatrix W ,double genvar, double resvar) {
        DoubleMatrix V = ZKZ.scalarMult(genvar);
        int n = V.numberOfRows();
        for(int i = 0; i < n; i++) {
            V.set(i, i, V.get(i, i)+W.get(i,i) * resvar);
        }
        return V;
    }

    public DistanceMatrix calculateDistanceFromKin(DistanceMatrix kin) {
        int n = kin.getSize();
        double max = kin.getDistance(0, 0);
        for (int i = 0; i < n; i++) {
            max = Math.max(max, kin.getDistance(i, i));
        }

        double constant;
        if (max > 2) {
            constant = max;
        } else if (max > 1) {
            constant = 2;
        } else {
            constant = 1;
        }

        DistanceMatrixBuilder builder = DistanceMatrixBuilder.getInstance(kin.getTaxaList());
        for (int r = 0; r < n; r++) {
            builder.set(r, r, constant - kin.getDistance(r, r));
            for (int c = r + 1; c < n; c++) {
                double newval = constant - kin.getDistance(r, c);
                builder.set(r, c, newval);
            }
        }
        return builder.build();
    }

    /**
     * @param missing			a BitSet with bits equal set when a value is missing in that row
     * @param phenotypeTaxa 	the taxa
     * @return 					a TaxaList with the taxa that are in both the kinship matrix and the phenotype.
     * Sets indices of taxa in missing that are not in the kinship matrix.
     */
    public TaxaList updateMissingWithKinship(BitSet missing, Taxon[] phenotypeTaxa) {
    	int n = phenotypeTaxa.length;
    	for (int i = 0; i < n; i++) {
    		int ndx = kinshipMatrix.whichIdNumber(phenotypeTaxa[i]);
    		if (ndx < 0) missing.fastSet(i);
    	}
    	Taxon[] nonMissingTaxa = AssociationUtils.getNonMissingValues(phenotypeTaxa, missing);
    	Set<Taxon> taxaSet = Arrays.stream(nonMissingTaxa).collect(Collectors.toCollection(HashSet::new));
    	
        return new TaxaListBuilder().addAll(taxaSet).build();
    }

    public Datum createResPhenotype(EMMAforDoubleMatrix emma, List<Taxon> taxa, String traitName) {
 
    	emma.calculateBlupsPredictedResiduals();
        DoubleMatrix res = emma.getRes();
        int nres = res.numberOfRows();
        
        float[] resarray = new float[nres];
        
    	for (int i = 0; i < nres; i++) resarray[i] = (float) res.get(i,0);
        
    	List<PhenotypeAttribute> attrList = new ArrayList<PhenotypeAttribute>();
    	List<ATTRIBUTE_TYPE> typeList = new ArrayList<Phenotype.ATTRIBUTE_TYPE>();
    	attrList.add(new TaxaAttribute(taxa));
    	typeList.add(ATTRIBUTE_TYPE.taxa);
    	attrList.add(new NumericAttribute(traitName, resarray, new OpenBitSet(nres)));
    	typeList.add(ATTRIBUTE_TYPE.data);
        Phenotype residualPheno = new PhenotypeBuilder().fromAttributeList(attrList, typeList).build().get(0);
        
        String name = String.format("Residuals for %s.", traitName);
        String comment = String.format("Residuals for %s calculated by MLM, no markers fit\nDataset: %s\n", traitName, datasetName);
        Datum output = new Datum(name, residualPheno, comment);
        return output;
     }
    
    class CompressedMLMResult {

        DoubleMatrix beta = null;
        double F = Double.NaN;
        double p = Double.NaN;
        double Fadd = Double.NaN;
        double padd = Double.NaN;
        double Fdom = Double.NaN;
        double pdom = Double.NaN;
        double r2 = Double.NaN;
        double addEffect = Double.NaN;
        double domEffect = Double.NaN;
        int modeldf;
        int markerdf;
        int ngroups;
    }

    public void setTestMarkers(boolean testMarkers) {
        this.testMarkers = testMarkers;
    }

}

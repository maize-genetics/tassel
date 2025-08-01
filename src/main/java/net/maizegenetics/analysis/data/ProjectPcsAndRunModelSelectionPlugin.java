/*
 * ProjectPcsAndRunModelSelectionPlugin
 */
package net.maizegenetics.analysis.data;

import java.awt.Frame;
import java.util.ArrayList;
import java.util.List;
import javax.swing.ImageIcon;
import net.maizegenetics.analysis.association.AssociationUtils;
import net.maizegenetics.analysis.numericaltransform.ImputationPlugin;
import net.maizegenetics.analysis.numericaltransform.NumericalGenotypePlugin;
import net.maizegenetics.dna.map.Chromosome;
import net.maizegenetics.dna.snp.FilterGenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.genotypecall.ProjectionGenotypeCallTable;
import net.maizegenetics.dna.snp.io.ProjectionGenotypeIO;
import net.maizegenetics.matrixalgebra.Matrix.DoubleMatrix;
import net.maizegenetics.matrixalgebra.Matrix.DoubleMatrixFactory;
import net.maizegenetics.phenotype.NumericAttribute;
import net.maizegenetics.phenotype.Phenotype;
import net.maizegenetics.phenotype.Phenotype.ATTRIBUTE_TYPE;
import net.maizegenetics.phenotype.PhenotypeAttribute;
import net.maizegenetics.phenotype.PhenotypeBuilder;
import net.maizegenetics.phenotype.TaxaAttribute;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.PluginEvent;
import net.maizegenetics.plugindef.PluginParameter;
import net.maizegenetics.stats.PCA.PrinComp;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.util.OpenBitSet;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 *
 * @author Alex Lipka
 *
 * This should enable users to read in a projection alignment, run a PCA within
 * a given window, and then conduct model selection
 */
public class ProjectPcsAndRunModelSelectionPlugin extends AbstractPlugin {

    private static final Logger myLogger = LogManager.getLogger(ProjectPcsAndRunModelSelectionPlugin.class);

    private PluginParameter<String> myRecombinationBreakpoints = new PluginParameter.Builder<>("recombinationBreakpoints", null, String.class).required(true).inFile()
            .description("").build();

    private GenotypeTable myHighDensityMarkersGenotypeTable = null;

    private GenotypeTable myCharacterAlignment;
    private double minRequiredData = 0.00;

    /**
     * Creates a new instance of ProjectPcsAndRunModelSelectionLPlugin
     */
    public ProjectPcsAndRunModelSelectionPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    @Override
    protected void preProcessParameters(DataSet input) {
        if (input == null) {
            throw new IllegalArgumentException("ProjectPcsAndRunModelSelectionPlugin: preProcessParameters: Please select one Genotype Table.");
        }
        List<Datum> genotypeTables = input.getDataOfType(GenotypeTable.class);
        if (genotypeTables.size() == 1) {
            myHighDensityMarkersGenotypeTable = (GenotypeTable) genotypeTables.get(0).getData();
        } else {
            throw new IllegalArgumentException("ProjectPcsAndRunModelSelectionPlugin: preProcessParameters: Please select one Genotype Table.");
        }
    }

    @Override
    public DataSet processData(DataSet input) {
        try {
            return loadFile(myRecombinationBreakpoints.value(), myHighDensityMarkersGenotypeTable);
        } catch (Exception e) {
            throw new IllegalStateException("ProjectPcsAndRunModelSelectionPlugin: processData: Problem loading: " + myRecombinationBreakpoints.value() + "\n" + e.getMessage());
        } finally {
            fireProgress(100);
        }

    }

    // The following getters and setters were auto-generated.
    // Please use this method to re-generate.
    //
    // public static void main(String[] args) {
    //     GeneratePluginCode.generate(ProjectionLoadPlugin.class);
    // }
    /**
     * Convenience method to run plugin with one return object.
     */
    public GenotypeTable runPlugin(DataSet input) {
        return (GenotypeTable) performFunction(input).getData(0).getData();
    }

    /**
     * Recombination Breakpoints
     *
     * @return Recombination Breakpoints
     */
    public String recombinationBreakpoints() {
        return myRecombinationBreakpoints.value();
    }

    /**
     * Set Recombination Breakpoints. Recombination Breakpoints
     *
     * @param value Recombination Breakpoints
     *
     * @return this plugin
     */
    public ProjectPcsAndRunModelSelectionPlugin recombinationBreakpoints(String value) {
        myRecombinationBreakpoints = new PluginParameter<>(myRecombinationBreakpoints, value);
        return this;
    }

    /**
     * Icon for this plugin to be used in buttons, etc.
     *
     * @return ImageIcon
     */
    public ImageIcon getIcon() {
        return null;
    }

    /**
     * Button name for this plugin to be used in buttons, etc.
     *
     * @return String
     */
    public String getButtonName() {
        return "Load Projection Alignment";
    }

    /**
     * Tool Tip Text for this plugin
     *
     * @return String
     */
    public String getToolTipText() {
        return "Load Projection Alignments";
    }

    public DataSet loadFile(String theRecombinationBreakpoints, GenotypeTable theHighDensityMarkers) {
        Datum test = new Datum("Full", theHighDensityMarkers, null);
        DataSet tests = new DataSet(test, this);
        fireDataSetReturned(new PluginEvent(tests, ProjectPcsAndRunModelSelectionPlugin.class));

        //Calcualte PCs across the NAM founders
        System.out.println("------------------------Calculating the PCs among the NAM founders--------------");
        Chromosome[] chr = theHighDensityMarkers.chromosomes();
        ArrayList<String> chrVector = new ArrayList<String>();
        ArrayList<Double> posVector = new ArrayList<Double>();//You can also use an ArrayList. it has an "add()" and "get()" method
        ArrayList<Double> startPosVector = new ArrayList<Double>();
        ArrayList<Double> endPosVector = new ArrayList<Double>();
        int increment = 10000;
        int[] selectedColumns = new int[]{0, 1, 2, 3, 4};
        DoubleMatrix PCResults = calculatePCsAcrossNAMFounders(chr, theHighDensityMarkers,
                chrVector, posVector, startPosVector, endPosVector, increment, selectedColumns);

        DataSet tdr = displayNamPCsOnTASSELGUI(PCResults, chrVector, posVector, startPosVector, endPosVector,
                theHighDensityMarkers);

        System.out.println("------------------------Done:- Calculating the PCs among the NAM founders--------------");

        //Create the projeciton alignment
        System.out.println("------------------------Creating the projection alignment--------------");
        GenotypeTable theAlignmentForGenotype = null;
        try {
            theAlignmentForGenotype = ProjectionGenotypeIO.getInstance(theRecombinationBreakpoints, theHighDensityMarkers);
        } catch (Exception e) {
            e.printStackTrace();
            return null;
        }
        System.out.println("------------------------Done:- Creating the projection alignment--------------");

        //Project the PCs onto the NAM
        System.out.println("------------------------Projecting PCs onto the NAM popluation--------------");
        DoubleMatrix ProjectedPCs = projectPCsOntoNAMFounders(theAlignmentForGenotype, PCResults,
                chrVector, posVector, theHighDensityMarkers, chr);
        System.out.println("------------------------Done:- Projecting PCs onto the NAM popluation--------------");

        System.out.println("------------------------Displaying Results on TASSEL GUI--------------");
        DataSet tds = displayProjectedPCsOnTASSELGUI(ProjectedPCs, chrVector, posVector, startPosVector,
                endPosVector, theAlignmentForGenotype);
        System.out.println("------------------------Done: Displaying Results on TASSEL GUI--------------");

        //fireDataSetReturned(new PluginEvent(tds, ProjectPcsAndRunModelSelectionPlugin.class));
        return tds;

    }

    public DoubleMatrix calculatePCsAcrossNAMFounders(Chromosome[] chr,
            GenotypeTable theGenotypesForCalculatingPCs, ArrayList<String> chrVector,
            ArrayList<Double> posVector, ArrayList<Double> startPosVector,
            ArrayList<Double> endPosVector, int increment, int[] selectedColumns) {

        GenotypeTable theGenotypesForCalculatingPCsOneChr = theGenotypesForCalculatingPCs;
        DoubleMatrix PCResults = null;
        int[] chrStartAndStop = new int[2];

        for (int i = 0; i < chr.length; i++) {
            chrStartAndStop = theGenotypesForCalculatingPCs.firstLastSiteOfChromosome(chr[i]);
            theGenotypesForCalculatingPCsOneChr = FilterGenotypeTable.getInstance(theGenotypesForCalculatingPCs, chrStartAndStop[0], chrStartAndStop[1]);

            int[] positions = theGenotypesForCalculatingPCsOneChr.physicalPositions();
            for (int j = 0; j < positions.length; j += increment) {
                int diffBetweenIncrementAndIndexj = positions.length - j;

                int startPos = positions[j];
                int endPos;
                if (diffBetweenIncrementAndIndexj >= increment) {
                    endPos = positions[j + increment];
                } else {
                    endPos = positions[(positions.length - 1)];
                }

                int myStart = theGenotypesForCalculatingPCsOneChr.siteOfPhysicalPosition(startPos, chr[i]);
                int myEnd = theGenotypesForCalculatingPCsOneChr.siteOfPhysicalPosition(endPos, chr[i]);

                GenotypeTable theGenotypesForCalculatingPCsReduced = theGenotypesForCalculatingPCs;
                theGenotypesForCalculatingPCsReduced = FilterGenotypeTable.getInstance(theGenotypesForCalculatingPCsReduced, myStart, myEnd);

                Datum test1 = new Datum("Reduced", theGenotypesForCalculatingPCsReduced, null);
                DataSet test1s = new DataSet(test1, this);
                //fireDataSetReturned(new PluginEvent(test1s, ProjectPcsAndRunModelSelectionPlugin.class));     

                //Create a numeric data set
                //SimplePhenotype numericalGenotypesForCalculatingPCs = NumericalGenotypePlugin.collapseTransform(theGenotypesForCalculatingPCsReduced);
                NumericalGenotypePlugin NGPConverter = new NumericalGenotypePlugin();
                GenotypeTable theGenotypesForCalculatingPCsReducedPartTwo = NGPConverter.setAlternateMinorAllelesToMinor(theGenotypesForCalculatingPCsReduced);
                ImputationPlugin imputor = new ImputationPlugin(null, false);
                imputor.by_mean(true);
                DataSet genoData = new DataSet(new Datum("name", theGenotypesForCalculatingPCsReducedPartTwo, "no comment"), null);
                DataSet numericalData = imputor.processData(genoData);

                myCharacterAlignment = (GenotypeTable) numericalData.getData(0).getData();

                int ntaxa = myCharacterAlignment.numberOfTaxa();
                int nsites = myCharacterAlignment.numberOfSites();
                DoubleMatrix dataMatrix = DoubleMatrixFactory.DEFAULT.make(ntaxa, nsites);
                for (int t = 0; t < ntaxa; t++) {
                    for (int s = 0; s < nsites; s++) {
                        dataMatrix.set(t, s, myCharacterAlignment.referenceProbability(t, s));
                    }
                }

                PrinComp myPrinComp = new PrinComp(dataMatrix, PrinComp.PC_TYPE.cov);
                //Use KNN to impute missing values  NOTE: Wait until the new KNN imputation code is up and running   

                //Datum ImpNumGeno4CalcPCsAsDatum = createImputedData();//You need to go into createImputedData() and fix things
                //DataSet ImpNumGeno4CalcPCs = new DataSet(ImpNumGeno4CalcPCsAsDatum, this);
                //fireDataSetReturned(new PluginEvent(ImpNumGeno4CalcPCs));
                //Obtain the PCs, which was ran in createImputedData()
                //Question: How do I get the first k PCs from myPCs?
                DoubleMatrix myPCs = myPrinComp.getPrincipalComponents();
                //System.out.println(myPCs.toString());

                if ((i == 0) & (j == 0)) {
                    PCResults = myPCs.getSelection(null, selectedColumns);
                } else {
                    PCResults = PCResults.concatenate(myPCs.getSelection(null, selectedColumns), false);
                    //DoubleFactory2D.dense.appendColumns() is concatenate() in DoubleMatrix
                    //.viewSelection(null,selectedColumns))is getSelection() in DoubleMatrix
                }
                //Append chrVector and posVector with the current chromosome and midpoint of the interval, respectively
                int posMidPoint = (startPos + endPos) / 2;
                for (int k = 0; k < selectedColumns.length; k++) {
                    chrVector.add(chr[i].toString());//TODO: Change this to the number of PCs per interval
                    posVector.add((double) posMidPoint);//TODO: Change this to the number of PCs per interval
                    startPosVector.add((double) startPos);//TODO: Change this to the number of PCs per interval
                    endPosVector.add((double) endPos);//TODO: Change this to the number of PCs per interval
                }

            }
        }
        return PCResults;
    }

    public DoubleMatrix projectPCsOntoNAMFounders(GenotypeTable theAlignmentForGenotype, DoubleMatrix PCResults,
            ArrayList<String> chrVector, ArrayList<Double> posVector, GenotypeTable theGenotypesForCalculatingPCs,
            Chromosome[] chr) {
        // theAlignmentForGenotype.chromosomalPosition(myEnd);
        ProjectionGenotypeCallTable pg = (ProjectionGenotypeCallTable) theAlignmentForGenotype.genotypeMatrix();

        // System.out.println("pg.numberOfTaxa(): "+ pg.numberOfTaxa());
        DoubleMatrix ProjectedPCs = null;
        for (int midpointPCSite = 0; midpointPCSite < chrVector.size(); midpointPCSite++) {
            double[] ProjectedPCColumn = new double[pg.numberOfTaxa()];
            DoubleMatrix ProjectedPCColumnAsDoubleMatrix = null;
            //Figure out the flanking sites of the midpoint of the interval of SNPs in which PCs were taken
            int[] leftAndRightFlankingMarkerSite = identifySitesOfFlankingMarkers(midpointPCSite, chrVector, posVector,
                    theGenotypesForCalculatingPCs, pg, chr);
            //IMPORTANT: this method parses out the sites on the given chromosome. Thus, the sites output are relative to
            // one chromosome at a time. This is why the code on lines 304-307 are there.
            for (int individual = 0; individual < pg.numberOfTaxa(); individual++) {

                int leftFlankingMarkerSite = leftAndRightFlankingMarkerSite[0];
                int rightFlankingMarkerSite = leftAndRightFlankingMarkerSite[1];

                //*****************Find out the donor parents for the two flanking sites, if such information is available
                double projectedPCElement;
                try {//If parental information is available at the sites
                    //Look at only the sites that are on the given chromosome

                    int[] theDonorsOnLeftFlank = pg.taxonDonors(individual, leftFlankingMarkerSite);

                    int[] theDonorsOnRightFlank = pg.taxonDonors(individual, rightFlankingMarkerSite);

                    //DoubleMatrix1D SpecificPCColumn = PCResults.viewColumn(midpointPCSite);
                    DoubleMatrix SpecificPCColumn = PCResults.column(midpointPCSite);

                    projectedPCElement = (0.25 * SpecificPCColumn.get(theDonorsOnLeftFlank[0], 0))
                            + (0.25 * SpecificPCColumn.get(theDonorsOnLeftFlank[1], 0))
                            + (0.25 * SpecificPCColumn.get(theDonorsOnRightFlank[0], 0))
                            + (0.25 * SpecificPCColumn.get(theDonorsOnRightFlank[1], 0));
                } catch (Exception e) {//If parental information is not available at the sites, indicate this by missing
                    projectedPCElement = Double.NaN;
                }
                //System.out.println("projectedPCElement "+ projectedPCElement);
                //System.out.println("SpecificPCColumn.get(theDonorsOnLeftFlank[0]) "+ SpecificPCColumn.get(theDonorsOnLeftFlank[0]));
                ProjectedPCColumn[individual] = projectedPCElement;
            }
            ProjectedPCColumnAsDoubleMatrix = DoubleMatrixFactory.DEFAULT.make(ProjectedPCColumn.length, 1, ProjectedPCColumn);
            //System.out.println(ProjectedPCColumnAsDoubleMatrix2D.toString());
            if (midpointPCSite == 0) {
                ProjectedPCs = ProjectedPCColumnAsDoubleMatrix;
            } else {
                ProjectedPCs = ProjectedPCs.concatenate(ProjectedPCColumnAsDoubleMatrix, false);
            }
        }
        return ProjectedPCs;
    }

    public DataSet displayProjectedPCsOnTASSELGUI(DoubleMatrix ProjectedPCs, ArrayList<String> chrVector,
            ArrayList<Double> posVector, ArrayList<Double> startPosVector,
            ArrayList<Double> endPosVector, GenotypeTable theAlignmentForGenotype) {

        TaxaList theTaxa = theAlignmentForGenotype.taxa();
        List<PhenotypeAttribute> myAttributes = new ArrayList<>();
        List<ATTRIBUTE_TYPE> types = new ArrayList<>();
        myAttributes.add(new TaxaAttribute(theTaxa));
        types.add(ATTRIBUTE_TYPE.taxa);
        Integer counter = 0;
        int ntaxa = theTaxa.numberOfTaxa();
        for (int i = 0; i < chrVector.size(); i++) {
            counter = counter + 1;
            if ((i > 0) && (!posVector.get(i).equals(posVector.get(i - 1)))) {
                counter = 1;
            }
            String name = "Chr_" + chrVector.get(i).toString() + "_Start_BP_"
                    + startPosVector.get(i).toString() + "_End_BP_"
                    + endPosVector.get(i).toString() + "_End_BP_" + "_PC_" + counter.toString();
            float[] data = AssociationUtils.convertDoubleArrayToFloat(ProjectedPCs.column(i).to1DArray());
            myAttributes.add(new NumericAttribute(name, data, new OpenBitSet(ntaxa)));
            types.add(ATTRIBUTE_TYPE.covariate);

        }
        double[][] ProjectedPCsAsDouble = new double[ProjectedPCs.numberOfRows()][ProjectedPCs.numberOfColumns()];
        for (int i = 0; i < ProjectedPCs.numberOfRows(); i++) {
            for (int j = 0; j < ProjectedPCs.numberOfColumns(); j++) {
                ProjectedPCsAsDouble[i][j] = ProjectedPCs.get(i, j);
            }
        }

        String ProjectedPCsReportName = "Projected PCs";
        String ProjectedPCsReportComments = "These are the projected PCs";
        Phenotype ProjectedPCsAsPhenotype = new PhenotypeBuilder().fromAttributeList(myAttributes, types).build().get(0);
        Datum ProjectedPCsDatum = new Datum(ProjectedPCsReportName, ProjectedPCsAsPhenotype, ProjectedPCsReportComments);
        DataSet ProjectedPCsDataSet = new DataSet(ProjectedPCsDatum, this);
        fireDataSetReturned(new PluginEvent(ProjectedPCsDataSet, ProjectPcsAndRunModelSelectionPlugin.class));
        return ProjectedPCsDataSet;
    }

    public DataSet displayNamPCsOnTASSELGUI(DoubleMatrix PCResults, ArrayList<String> chrVector,
            ArrayList<Double> posVector, ArrayList<Double> startPosVector,
            ArrayList<Double> endPosVector, GenotypeTable theGenotypesForCalculatingPCs) {

        TaxaList theTaxa = theGenotypesForCalculatingPCs.taxa();
        List<PhenotypeAttribute> myAttributes = new ArrayList<>();
        List<ATTRIBUTE_TYPE> types = new ArrayList<>();
        myAttributes.add(new TaxaAttribute(theTaxa));
        types.add(ATTRIBUTE_TYPE.taxa);
        Integer counter = 0;
        int ntaxa = theTaxa.numberOfTaxa();
        for (int i = 0; i < chrVector.size(); i++) {
            counter = counter + 1;
            if ((i > 0) && (!posVector.get(i).equals(posVector.get(i - 1)))) {
                counter = 1;
            }
            String name = "Chr_" + chrVector.get(i).toString() + "_Start_BP_" + startPosVector.get(i).toString() + "_End_BP_"
                    + endPosVector.get(i).toString() + "_PC_" + counter.toString();
            float[] data = AssociationUtils.convertDoubleArrayToFloat(PCResults.column(i).to1DArray());
            myAttributes.add(new NumericAttribute(name, data, new OpenBitSet(ntaxa)));
            types.add(ATTRIBUTE_TYPE.covariate);
        }

        String ProjectedPCsReportName = "PCs among NAM Founders";
        String ProjectedPCsReportComments = "PCs among NAM Founders";
        Phenotype ProjectedPCsAsPhenotype = new PhenotypeBuilder()
                .fromAttributeList(myAttributes, types)
                .build().get(0);
        Datum ProjectedPCsDatum = new Datum(ProjectedPCsReportName, ProjectedPCsAsPhenotype, ProjectedPCsReportComments);
        DataSet ProjectedPCsDataSet = new DataSet(ProjectedPCsDatum, this);
        fireDataSetReturned(new PluginEvent(ProjectedPCsDataSet, ProjectPcsAndRunModelSelectionPlugin.class));
        return ProjectedPCsDataSet;
    }

    public int[] identifySitesOfFlankingMarkers(int site, ArrayList<String> chrVector, ArrayList<Double> posVector,
            GenotypeTable theGenotypesForCalculatingPCs, ProjectionGenotypeCallTable pg,
            Chromosome[] chr) {
        //Look at the  "taxonDonors()" method within ProjectionGenotypeCallTable

        Chromosome testedChromosome = new Chromosome(chrVector.get(site));
        int[] chrStartAndStop = theGenotypesForCalculatingPCs.firstLastSiteOfChromosome(testedChromosome);
        GenotypeTable theGenotypesForCalculatingPCsOneChr = theGenotypesForCalculatingPCs;
        theGenotypesForCalculatingPCsOneChr = FilterGenotypeTable.getInstance(theGenotypesForCalculatingPCsOneChr, chrStartAndStop[0], chrStartAndStop[1]);

        //********************Get the flaking sites on right and left
        // Note: positive distance means the marker is to the right; negative distance means
        // the marker is to the left
        int leftFlankingMarkerSite = 0;
        int rightFlankingMarkerSite = 0;

        ArrayList distanceFromMidpointOfInterval = new ArrayList();
        ArrayList positiveDistanceFromMidpointOfInterval = new ArrayList();
        ArrayList negativeDistanceFromMidpointOfInterval = new ArrayList();
        for (int j = 0; j < theGenotypesForCalculatingPCsOneChr.numberOfSites(); j++) {
            Double testPosition = posVector.get(site);
            double distance = theGenotypesForCalculatingPCsOneChr.chromosomalPosition(j) - testPosition;
            distanceFromMidpointOfInterval.add(distance);
            if (distance > 0) {
                positiveDistanceFromMidpointOfInterval.add(distance);
            }
            if (distance < 0) {
                negativeDistanceFromMidpointOfInterval.add(distance);
            }
        }
        //Find out distance to the nearest flanking markers: NOTE THESE NEXT TWO FOR LOOPS MAY BE UNNECESSARY IF THE SNPS ARE SORTED
        // IN GENOTYPIC ORDER.
        double distanceToRightMarker = Double.MAX_VALUE;
        if (positiveDistanceFromMidpointOfInterval.size() > 0) {
            for (int j = 0; j < positiveDistanceFromMidpointOfInterval.size(); j++) {
                double positiveDistanceArrayElement = (double) positiveDistanceFromMidpointOfInterval.get(j);
                if (positiveDistanceArrayElement < distanceToRightMarker) {
                    distanceToRightMarker = positiveDistanceArrayElement;
                }
            }
        } else {
            distanceToRightMarker = 0;
        }

        double distanceToLeftMarker = Double.MAX_VALUE;
        if (negativeDistanceFromMidpointOfInterval.size() > 0) {
            for (int j = 0; j < negativeDistanceFromMidpointOfInterval.size(); j++) {
                double negativeDistanceArrayElement = (double) negativeDistanceFromMidpointOfInterval.get(j);
                negativeDistanceArrayElement = -1 * negativeDistanceArrayElement;
                if (negativeDistanceArrayElement < distanceToLeftMarker) {
                    distanceToLeftMarker = negativeDistanceArrayElement;
                }
            }
            distanceToLeftMarker = -1 * distanceToLeftMarker;
        } else {
            distanceToLeftMarker = 0;
        }

        //Obtain the sites of the flanking markers
        if (distanceToRightMarker != 0) {
            //Obtain the index of distanceFromMidpointOfInterval where the distance matches up. This
            // will be the site number
            rightFlankingMarkerSite = distanceFromMidpointOfInterval.indexOf(distanceToRightMarker);
        }
        if (distanceToLeftMarker != 0) {
            //Obtain the index of distanceFromMidpointOfInterval where the distance matches up. This
            // will be the site number
            leftFlankingMarkerSite = distanceFromMidpointOfInterval.indexOf(distanceToLeftMarker);
        }
        int[] leftAndRightFlankingMarkerSiteAndChrStartAndStop = new int[2];
        //We add chrStartAndStop[0] to these values so that pg will parse out the correct elements in pg.
        // i.e., left and rightFlankingMarkerSites are relative to one chromosome, while the
        // pg object is for all chromosomes. Thus, adding chrStartAndStop[0] to the bottom two values
        // ensures that the correct site number is being used. 
        leftAndRightFlankingMarkerSiteAndChrStartAndStop[0] = leftFlankingMarkerSite + chrStartAndStop[0];
        leftAndRightFlankingMarkerSiteAndChrStartAndStop[1] = rightFlankingMarkerSite + chrStartAndStop[0];
        return leftAndRightFlankingMarkerSiteAndChrStartAndStop;
    }      //End method here

    /*   public Datum createImputedData() {
     //int[] colsSelected = null;       // set of columns to be used to calculate distance (should be correlated columns)
     //colsSelected = tblTraits.getSelectedRows();
     //int colCount = colsSelected.length;
     //int includedCount = 0;
     //find all the rows with enough data to keep
     //int ntaxa = myCharacterAlignment.getNumberOfTaxa();
     int ntaxa = myCharacterAlignment.numberOfObservations();
     //int nsites = myCharacterAlignment.getNumberOfTraits();
     int nsites = myCharacterAlignment.numberOfAttributes();
     double[][] tempData = new double[ntaxa][nsites];
     for (int t = 0; t < ntaxa; t++) {
     for (int s = 0; s < nsites; s++) {
     //tempData[t][s] = myCharacterAlignment.getData(t, s);
     tempData[t][s] = (double) myCharacterAlignment.getValueAt(t, s);
     }
     }

     //See if there are any taxa with all missing marker data. If there are any, replace with the average
     // numeric marker value
     for (int i = 0; i < ntaxa; i++) {
     int count = 0;
     for (int j = 0; j < nsites; j++) {
     //Check to see if the tempData[i][j] is missing
     if (!Double.isNaN(tempData[i][j])) {
     break;
     }
     count++;
     }
     if (count == nsites) {
     for (int j = 0; j < nsites; j++) {
     //System.out.println("The " + i + "th taxa did not have any marker data");
     ArrayList<Double> columnValues = new ArrayList<Double>();
     for (int k = 0; k < ntaxa; k++) {
     columnValues.add(tempData[k][j]);
     }
     //Calculate the column average
     double theSum = 0;
     int theNumberOfInds = 0;
     for (int k = 0; k < ntaxa; k++) {
     if (!Double.isNaN(columnValues.get(k))) {
     theSum = theSum + columnValues.get(k);
     theNumberOfInds++;
     }
     }
     if (theNumberOfInds != 0) {
     tempData[i][j] = theSum / theNumberOfInds;
     } else {
     tempData[i][j] = 0.5;
     }

     //Set tempData[i][j] equal to the column average
     }
     }
     }

     //        int[] includedRowTemp = new int[myCharacterAlignment.getNumberOfTaxa()];
     //        for (int i = 0; i < myCharacterAlignment.getNumberOfTaxa(); i++) {
     //            double goodData = 0;
     //            for (int j = 0; j < colCount; j++) {
     //                if (!Double.isNaN(myCharacterAlignment.getData(i, colsSelected[j]))) {
     //                    goodData++;
     //                }
     //            }
     //            goodData = goodData / colCount;
     //            if (goodData >= minRequiredData) {
     //                includedRowTemp[includedCount++] = i;
     //            }
     //        }
     //        //rebuild the data set
     //        Taxon[] newIDs = new Taxon[includedCount];
     //
     //        int traitCount = colsSelected.length;
     //        java.util.List<Trait> newtraits = new ArrayList<Trait>();
     //        for (int t = 0; t < traitCount; t++) {
     //            newtraits.add(Trait.getInstance(myCharacterAlignment.getTrait(colsSelected[t])));
     //        }
     //
     //        double[][] tempData = new double[includedCount][colsSelected.length];
     //        for (int i = 0; i < includedCount; i++) {
     //            for (int j = 0; j < colCount; j++) {
     //                newIDs[i] = myCharacterAlignment.getTaxa().get(includedRowTemp[i]);
     //                tempData[i][j] = myCharacterAlignment.getData(includedRowTemp[i], colsSelected[j]);
     //            }
     //        }
     //      for(int j = 0; j < colCount; j++){
     //              newTraits[j]=aCharacterAlignment.getTraitName(colsSelected[j]);
     //              newEnvs[j]=aCharacterAlignment.getEnvironmentName(colsSelected[j]);
     //          }
     int kNeighbors = 3;
     double[][] theImputedData = KNN.impute(tempData, kNeighbors, true, true);
     DoubleMatrix values = DoubleMatrixFactory.DEFAULT.make(theImputedData);
     myPrinComp = new PrinComp(values, PrinComp.PC_TYPE.cov);

     //SimplePhenotype sca = new SimplePhenotype(new SimpleIdGroup(newIDs), theImputedData, aCharacterAlignment.getFactorNameCopy(), newtraits);
     //TaxaList tL = new TaxaListBuilder().addAll(myCharacterAlignment.getTaxa()).build();
     TaxaList tL = new TaxaListBuilder().addAll(myCharacterAlignment.taxa()).build();
     //SimplePhenotype sca = new SimplePhenotype(tL, myCharacterAlignment.getTraits(), theImputedData);
        
     Phenotype sca = new PhenotypeBuilder().fromPhenotypeList(myCharacterAlignment.).build();
        
     Phenotype(tL, myCharacterAlignment.attribute(ntaxa), theImputedData);
     StringWriter sw = new StringWriter();
     //sca.report(new PrintWriter(sw));
     String theComment = sw.toString() + "\nImputed Phenotypic Values." + "\nTaxa with insufficient data: " + (myCharacterAlignment.taxa() - sca.getNumberOfTaxa()) + "\nK = " + kNeighbors + minRequiredData + "% cutoff):\n";
     String theName = "Imputed_Data";
     Datum result = new Datum(theName, sca, theComment);
     return result;
     }*/
}

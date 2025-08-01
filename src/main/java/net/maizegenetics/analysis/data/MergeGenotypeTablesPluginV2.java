package net.maizegenetics.analysis.data;

import java.awt.Frame;
import java.net.URL;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.function.BiFunction;

import javax.swing.ImageIcon;
import javax.swing.JOptionPane;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import net.maizegenetics.dna.WHICH_ALLELE;
import net.maizegenetics.dna.map.GeneralPosition;
import net.maizegenetics.dna.map.Position;
import net.maizegenetics.dna.map.PositionList;
import net.maizegenetics.dna.map.PositionListBuilder;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.MergedGenotypeTable;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.PluginEvent;
import net.maizegenetics.plugindef.PluginParameter;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListBuilder;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.util.GeneralAnnotation;

public class MergeGenotypeTablesPluginV2 extends net.maizegenetics.plugindef.AbstractPlugin {
    
    private static final Logger myLogger = LogManager.getLogger(MergeGenotypeTablesPlugin.class);

    
    public static enum MERGE_TYPES {
      Intersect,Union,LeftJoin,RightJoin  
    };
    
    public static enum DEPTH_MERGE {
        Additive,MaxDepth
    };
    
    public static enum CALL_MERGE {
        Depth,RefereceProbability,Dosage,AlleleProbability
    };
    
    private PluginParameter<MERGE_TYPES> taxaMergeSelection = new PluginParameter.Builder<>("taxaMerge", MERGE_TYPES.Intersect, MERGE_TYPES.class)
                                                                        .description("Selection for Taxa Merge Rule")
                                                                        .guiName("Taxa Merge Rule")
                                                                        .build();
    private PluginParameter<MERGE_TYPES> positionMergeSelection = new PluginParameter.Builder<>("posMerge", MERGE_TYPES.Intersect, MERGE_TYPES.class)
                                                                        .description("Selection for Position Merge Rule")
                                                                        .guiName("Position Merge Rule")
                                                                        .build();
    private PluginParameter<CALL_MERGE> callMergeSelection = new PluginParameter.Builder<>("callMerge", CALL_MERGE.Depth, CALL_MERGE.class)
                                                                        .description("Selection for Call Merging")
                                                                        .guiName("Call Merge Rule")
                                                                        .build();
    private PluginParameter<DEPTH_MERGE> depthMergeSelection = new PluginParameter.Builder<>("depthMerge", DEPTH_MERGE.Additive, DEPTH_MERGE.class)      
                                                                        .description("Selection for How to Merge Depths")
                                                                        .guiName("Depth Merge Rule")
                                                                        .build();
    
    public MergeGenotypeTablesPluginV2(Frame parentFrame, boolean isInteractive) { 
        super(parentFrame, isInteractive);
    }
    
    public DataSet processData(DataSet input) {
        //do whatever your plugin does
        List<Datum> inputs = input.getDataOfType(GenotypeTable.class);

        if ((inputs == null) || (inputs.size() < 2)) {
            if (isInteractive()) {
                JOptionPane.showMessageDialog(getParentFrame(), "Must select at least two alignments.");
            } else {
                myLogger.warn("performFunction: Must select at least two alignments.");
            }
            return null;
        }

        try {
            GenotypeTable[] alignments = new GenotypeTable[inputs.size()];
            for (int i = 0; i < inputs.size(); i++) {
                alignments[i] = (GenotypeTable) ((Datum) inputs.get(i)).getData();
            }

            GenotypeTable merged = mergeGenotypeTables(alignments);
            DataSet result = new DataSet(new Datum("Merged Genotype Table", merged, null), this);

            fireDataSetReturned(new PluginEvent(result, MergeGenotypeTablesPlugin.class));

            return result;
        } finally {
            fireProgress(100);
        }
         
    }
    
    private GenotypeTable mergeGenotypeTables(GenotypeTable[] genotypeTables) {
        BiFunction<List,List,List> taxaMergeRule = getMergeLambda(taxaMergeSelection());
        BiFunction<List,List,List> posMergeRule = getMergeLambda(positionMergeSelection());
        return MergedGenotypeTable.getInstance(genotypeTables,taxaMergeRule,posMergeRule);
        
        //return CombineGenotypeTable.getInstance(genotypeTables,taxaMergeRule,posMergeRule);
        //return CombineGenotypeTable.getInstance(genotypeTables, true);
    }
    
    //TODO Write Unit Test for merging Taxa
    public Taxon mergeTaxon(Taxon a, Taxon b) {
        GeneralAnnotation gaA = a.getAnnotation();
        GeneralAnnotation gaB = b.getAnnotation();
        Taxon.Builder tb = new Taxon.Builder(a.getName());
        Set<String> masterKeySet = gaA.getAnnotationKeys();
        masterKeySet.addAll(gaB.getAnnotationKeys());
        //TODO check with biologists to figure out what to do where there can be only one value for the annotation biologically
        //TODO Figure out what annotations can only have one value
        //Remove any duplicate values
        for(String key: masterKeySet) {
            String[] taxaAEntries = gaA.getTextAnnotation(key);
            String[] taxaBEntries = gaB.getTextAnnotation(key);
            HashSet<String> entrySet = new HashSet<String>();
            for(String entryA : taxaAEntries) {
                entrySet.add(entryA);
            }
            for(String entryB : taxaBEntries) {
                entrySet.add(entryB);
            }
            String[] masterEntries = entrySet.toArray(new String[0]);
            for(String entry : masterEntries) {
                tb.addAnno(key, entry);
            }
        } 
        return tb.build();
    }
    
    //TODO Write Unit Test for merging positions
    public Position mergePosition(Position a, Position b) {
        System.out.println("Merging Position");
        GeneralPosition.Builder mergedBuilder = new GeneralPosition.Builder(a.getChromosome(), a.getPosition());
        mergedBuilder.snpName(a.getSNPID());
        mergedBuilder.insertionPosition(a.getInsertionPosition());
        mergedBuilder.strand(a.getStrand());
        mergedBuilder.maf(a.getGlobalMAF());
        mergedBuilder.siteCoverage(a.getGlobalSiteCoverage());
        for (WHICH_ALLELE alleleType : WHICH_ALLELE.values()) {
            mergedBuilder.allele(alleleType, a.getAllele(alleleType));
        }
        mergedBuilder.nucleotide(a.isNucleotide());
        mergedBuilder.indel(a.isIndel());
        
        GeneralAnnotation gaA = a.getAnnotation();
        GeneralAnnotation gaB = b.getAnnotation();
        Set<String> masterKeySet = gaA.getAnnotationKeys();
        masterKeySet.addAll(gaB.getAnnotationKeys());
        for(String key: masterKeySet) {
            String[] posAEntries = gaA.getTextAnnotation(key);
            String[] posBEntries = gaB.getTextAnnotation(key);
            HashSet<String> entrySet = new HashSet<String>();
            for(String entryA : posAEntries) {
                entrySet.add(entryA);
            }
            for(String entryB : posBEntries) {
                entrySet.add(entryB);
            }
            String[] masterEntries = entrySet.toArray(new String[0]);
            for(String entry : masterEntries) {
                mergedBuilder.addAnno(key, entry);
            }
        } 
        return mergedBuilder.build();
    }
    
    private BiFunction<List,List,List> getMergeLambda(MERGE_TYPES mergeType) {
        BiFunction<List,List,List> intersectFunction = (a,b) -> {
            if(a instanceof TaxaList && b instanceof TaxaList) {
                TaxaListBuilder tlb = new TaxaListBuilder();
                for(Taxon taxaA :(TaxaList)a) {
                    int bIndex = -1;
                    if((bIndex = ((TaxaList)b).indexOf(taxaA)) !=-1) {
                        tlb.add(mergeTaxon(taxaA,((TaxaList)b).get(bIndex)));
                    }
                }
                return tlb.build();
            }
            else if(a instanceof PositionList && b instanceof PositionList) {
                PositionListBuilder plb = new PositionListBuilder();
                for(Position posA:(PositionList)a) {
                    int bIndex = -1;
                    if((bIndex = ((PositionList)b).indexOf(posA)) != -1) {
                        plb.add(mergePosition(posA,((PositionList)b).get(bIndex)));
                    }
                }
                System.out.println("Size of Pos "+plb.size());
                return plb.build();
            }
            else {
                //Handle generic lists
                List<Object> listBuilder = new ArrayList<Object>();
                for(Object objA : a) {
                    if(b.contains(objA)) {
                        listBuilder.add(objA);
                    }
                }
                return listBuilder;
            }
            
        };
        
        BiFunction<List,List,List> leftJoinFunction = (a,b) -> {
            if(a instanceof TaxaList && b instanceof TaxaList) {
                TaxaListBuilder tlb = new TaxaListBuilder();
                for(Taxon taxaA :(TaxaList)a) {
                    int bIndex = -1;
                    if((bIndex = ((TaxaList)b).indexOf(taxaA)) != -1) {
                        tlb.add(mergeTaxon(taxaA,((TaxaList)b).get(bIndex)));
                    }
                    else {
                        tlb.add(taxaA);
                    }
                }
                return tlb.build();
            }
            else if(a instanceof PositionList && b instanceof PositionList) {
                PositionListBuilder plb = new PositionListBuilder();
                for(Position posA :(PositionList)a) {
                    int bIndex = -1;
                    if((bIndex = ((PositionList)b).indexOf(posA)) != -1) {
                        plb.add(mergePosition(posA,((PositionList)b).get(bIndex)));
                    }
                    else {
                        plb.add(posA);
                    }
                }
                return plb.build();
            }
            else {
                //Handle generic lists
                List<Object> listBuilder = new ArrayList<Object>();
                for(Object objA : a) {
                    if(b.contains(objA)) {
                        listBuilder.add(objA);
                    }
                    else {
                        listBuilder.add(objA);
                    }
                }
                return listBuilder;
            }
        };
       

        BiFunction<List,List,List> rightJoinFunction = (a,b) -> {
            if(a instanceof TaxaList && b instanceof TaxaList) {
                TaxaListBuilder tlb = new TaxaListBuilder();
                for(Taxon taxaB :(TaxaList)b) {
                    int aIndex = -1;
                    if((aIndex = ((TaxaList)a).indexOf(taxaB)) != -1) {
                        tlb.add(mergeTaxon(taxaB,((TaxaList)a).get(aIndex)));
                    }
                    else {
                        tlb.add(taxaB);
                    }
                }
                return tlb.build();
            }
            else if(a instanceof PositionList && b instanceof PositionList) {
                PositionListBuilder plb = new PositionListBuilder();
                for(Position posB :(PositionList)b) {
                    int aIndex = -1;
                    if((aIndex = ((PositionList)a).indexOf(posB)) != -1) {
                        plb.add(mergePosition(posB,((PositionList)a).get(aIndex)));
                    }
                    else {
                        plb.add(posB);
                    }
                }
                return plb.build();
            }
            else {
                //Handle generic lists
                List<Object> listBuilder = new ArrayList<Object>();
                for(Object objB : b) {
                    if(b.contains(objB)) {
                        listBuilder.add(objB);
                    }
                    else {
                        listBuilder.add(objB);
                    }
                }
                return listBuilder;
            }
        };
        
        BiFunction<List,List,List> unionFunction = (a,b) -> {
            if(a instanceof TaxaList && b instanceof TaxaList) {
                TaxaListBuilder tlb = new TaxaListBuilder();
                for(Taxon taxaA : (TaxaList)a) {
                    int bIndex = -1;
                    if((bIndex = ((TaxaList)b).indexOf(taxaA)) !=-1) {
                        tlb.add(mergeTaxon(taxaA,((TaxaList)b).get(bIndex)));
                    }
                    else {
                        tlb.add(taxaA);
                    }  
                }
                for(Taxon taxaB :(TaxaList)b) {
                    int aIndex = -1;
                    if((aIndex = ((TaxaList)a).indexOf(taxaB)) ==-1) {
                        tlb.add(taxaB);
                    }
                }
                return tlb.build();
            }
            else if(a instanceof PositionList && b instanceof PositionList) {
                PositionListBuilder plb = new PositionListBuilder();
                for(Position posA : (PositionList)a) {
                    int bIndex = -1;
                    if((bIndex = ((PositionList)b).indexOf(posA)) !=-1) {
                        plb.add(mergePosition(posA,((PositionList)b).get(bIndex)));
                    }
                    else {
                        plb.add(posA);
                    }  
                }
                for(Position posB :(PositionList)b) {
                    int aIndex = -1;
                    if((aIndex = ((PositionList)a).indexOf(posB)) ==-1) {
                        plb.add(posB);
                    }
                }
                return plb.build();
            }
            else {
                //Handle generic lists
                List<Object> listBuilder = new ArrayList<Object>();
                for(Object objA : a) {                   
                    listBuilder.add(objA);
                }
                for(Object objB : b) {
                    if(!listBuilder.contains(objB)) {
                        listBuilder.add(objB);
                    }
                }
                return listBuilder;
            } 
        };
        
        switch(mergeType) {
            case Union:
                return unionFunction;
            case LeftJoin:
                return leftJoinFunction;
            case RightJoin:
                return rightJoinFunction;
            case Intersect:
            default:
                return intersectFunction;
        }
    }
    
    @Override
    public ImageIcon getIcon() {
        URL imageURL = MergeGenotypeTablesPlugin.class.getResource("/net/maizegenetics/analysis/images/Merge.gif");
        if (imageURL == null) {
            return null;
        } else {
            return new ImageIcon(imageURL);
        }
    }
    @Override
    public String getButtonName() {
        return "Merge GenotypeTable V2"; 
    }
    @Override
    public String getToolTipText() {
        return "Merges Genotype Tables";
    }

    @Override
    public String pluginDescription() {
        return "This plugin takes a two or more GenotypeTables and Creates a CombinedGenotypeTable based on user-specified rules";
    }
    
    // The following getters and setters were auto-generated.
    // Please use this method to re-generate.
    //
    // public static void main(String[] args) {
    //     GeneratePluginCode.generate(MergeGenotypeTablesPluginV2.class);
    // }

    /**
     * Convenience method to run plugin with one return object.
     */
    // TODO: Replace <Type> with specific type.
    /*
    public <Type> runPlugin(DataSet input) {
        return (<Type>) performFunction(input).getData(0).getData();
    }
    */

    /**
     * Selection for Taxa Merge Rule
     *
     * @return Taxa Merge Rule
     */
    public MERGE_TYPES taxaMergeSelection() {
        return taxaMergeSelection.value();
    }

    /**
     * Set Taxa Merge Rule. Selection for Taxa Merge Rule
     *
     * @param value Taxa Merge Rule
     *
     * @return this plugin
     */
    public MergeGenotypeTablesPluginV2 taxaMergeSelection(MERGE_TYPES value) {
        taxaMergeSelection = new PluginParameter<>(taxaMergeSelection, value);
        return this;
    }

    /**
     * Selection for Position Merge Rule
     *
     * @return Position Merge Rule
     */
    public MERGE_TYPES positionMergeSelection() {
        return positionMergeSelection.value();
    }

    /**
     * Set Position Merge Rule. Selection for Position Merge
     * Rule
     *
     * @param value Position Merge Rule
     *
     * @return this plugin
     */
    public MergeGenotypeTablesPluginV2 positionMergeSelection(MERGE_TYPES value) {
        positionMergeSelection = new PluginParameter<>(positionMergeSelection, value);
        return this;
    }

    /**
     * Selection for Call Merging
     *
     * @return Call Merge Rule
     */
    public CALL_MERGE callMergeSelection() {
        return callMergeSelection.value();
    }

    /**
     * Set Call Merge Rule. Selection for Call Merging
     *
     * @param value Call Merge Rule
     *
     * @return this plugin
     */
    public MergeGenotypeTablesPluginV2 callMergeSelection(CALL_MERGE value) {
        callMergeSelection = new PluginParameter<>(callMergeSelection, value);
        return this;
    }
    
    /**
     * Selection for How to Merge Depths
     *
     * @return Depth Merge Rule
     */
    public DEPTH_MERGE depthMergeSelection() {
        return depthMergeSelection.value();
    }

    /**
     * Set Depth Merge Rule. Selection for How to Merge Depths
     *
     * @param value Depth Merge Rule
     *
     * @return this plugin
     */
    public MergeGenotypeTablesPluginV2 depthMergeSelection(DEPTH_MERGE value) {
        depthMergeSelection = new PluginParameter<>(depthMergeSelection, value);
        return this;
    }

    
}

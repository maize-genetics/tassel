/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.analysis.imputation;

import java.awt.Frame;
import java.io.File;
import java.util.ArrayList;
import javax.swing.ImageIcon;
import net.maizegenetics.dna.snp.ExportUtils;
import net.maizegenetics.dna.snp.FilterGenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableBuilder;
import net.maizegenetics.dna.snp.GenotypeTableUtils;
import net.maizegenetics.dna.snp.ImportUtils;
import net.maizegenetics.dna.snp.NucleotideAlignmentConstants;
import net.maizegenetics.dna.snp.genotypecall.GenotypeCallTable;
import net.maizegenetics.dna.snp.genotypecall.GenotypeCallTableBuilder;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.PluginParameter;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 *
 * @author kelly
 */
public class RemoveIndelsForBeaglePlugin extends net.maizegenetics.plugindef.AbstractPlugin{
    private PluginParameter<String> inFile= new PluginParameter.Builder<>("inFile",null,String.class).guiName("inFileName").inFile().required(true)
            .description("The input file to prepare for beagle input").build();
    private PluginParameter<String> outFile= new PluginParameter.Builder<>("outFile",null,String.class).guiName("outFileName").outFile().required(true)
            .description("the output file name. If not VCF, will be append .vcf to the name").build();
    private PluginParameter<Boolean> retainMono= new PluginParameter.Builder<>("retainMono",false,Boolean.class).guiName("Retain monomorphic sites")
            .description("Retain monomorphic sites in the output file").build();
    
    private static final Logger myLogger = LogManager.getLogger(RemoveIndelsForBeaglePlugin.class);

    @Override
    protected void postProcessParameters() {
    }
    
    public RemoveIndelsForBeaglePlugin() {
        super(null, false);
    }

    public RemoveIndelsForBeaglePlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }
    
    @Override
    public DataSet processData(DataSet input) {
        GenotypeTable genos= removeIndelsForBeagle(); String out= outFile.value();
        if (!outFile.value().contains(".vcf")) out= outFile.value()+".vcf.gz";
        ExportUtils.writeToVCF(genos, out, false);
        System.out.println("Wrote "+out+" with "+genos.numberOfSites()+" sites remaining");
    
        return null;
    }
       
    public GenotypeTable removeIndelsForBeagle() {
        GenotypeTable genos= ImportUtils.readGuessFormat(inFile.value());
        GenotypeCallTableBuilder newGenos= GenotypeCallTableBuilder.getInstanceCopy(genos.genotypeMatrix());
        ArrayList<String> keepSites= new ArrayList<>();
        byte gapInsert= GenotypeTableUtils.getUnphasedDiploidValue(NucleotideAlignmentConstants.GAP_ALLELE, NucleotideAlignmentConstants.INSERT_ALLELE);
        for (int site = 0; site < genos.numberOfSites(); site++) {
            byte minmaj= GenotypeTableUtils.getUnphasedDiploidValue(genos.majorAllele(site), genos.minorAllele(site));
            if (GenotypeTableUtils.isPartiallyEqual(gapInsert, minmaj)) continue;
            if (!retainMono.value() && !genos.isPolymorphic(site)) continue;
            keepSites.add(genos.siteName(site));
            if (genos.alleles(site).length>2) {
                byte badGeno= GenotypeTableUtils.getDiploidValue(NucleotideAlignmentConstants.GAP_ALLELE, NucleotideAlignmentConstants.INSERT_ALLELE);
                for (int taxon = 0; taxon < genos.numberOfTaxa(); taxon++) {
                    if (GenotypeTableUtils.isPartiallyEqual(genos.genotype(taxon, site),badGeno)) {
                        newGenos.setBase(taxon, site, GenotypeTable.UNKNOWN_DIPLOID_ALLELE);
                    }
                }
            }
        }
        GenotypeCallTable newGenoCallTable= newGenos.build();
        GenotypeTable newGT=GenotypeTableBuilder.getInstance(newGenoCallTable, genos.positions(), genos.taxa());
        GenotypeTable fa= GenotypeTableBuilder.getGenotypeCopyInstance(FilterGenotypeTable.getInstance(newGT, keepSites));
        return fa;
    }
    
     @Override
    public ImageIcon getIcon() {
        return null;
    }

    @Override
    public String getButtonName() {
        return "Remove indels for input to Beagle v.4";
    }

    @Override
    public String getToolTipText() {
        return "Removes sites with an indel as major and minor and changes indels in the third allele state to missing";
    }
}

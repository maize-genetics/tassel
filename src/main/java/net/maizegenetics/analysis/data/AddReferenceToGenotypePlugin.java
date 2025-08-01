package net.maizegenetics.analysis.data;

import net.maizegenetics.dna.WHICH_ALLELE;
import net.maizegenetics.dna.map.GeneralPosition;
import net.maizegenetics.dna.map.GenomeSequence;
import net.maizegenetics.dna.map.GenomeSequenceBuilder;
import net.maizegenetics.dna.map.Position;
import net.maizegenetics.dna.map.PositionList;
import net.maizegenetics.dna.map.PositionListBuilder;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableBuilder;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.PluginParameter;
import net.maizegenetics.util.Utils;

import javax.swing.*;
import java.awt.*;
import java.util.List;

/**
 * @author Terry Casstevens Created November 09, 2018
 */
public class AddReferenceToGenotypePlugin extends AbstractPlugin {

    private PluginParameter<String> myReference = new PluginParameter.Builder<>("reference", null, String.class)
            .description("Input Reference Fasta")
            .required(true)
            .inFile()
            .build();

    public AddReferenceToGenotypePlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    @Override
    public DataSet processData(DataSet input) {

        List<Datum> temp = input.getDataOfType(GenotypeTable.class);
        if (temp.size() != 1) {
            throw new IllegalArgumentException("AddReferenceToGenotypePlugin: processData: must select one genotype table");
        }
        GenotypeTable orig = (GenotypeTable) temp.get(0).getData();

        GenomeSequence reference = GenomeSequenceBuilder.instance(reference());

        PositionList origPositions = orig.positions();

        PositionListBuilder builder = new PositionListBuilder();
        builder.genomeVersion(Utils.getFilename(reference()));

        origPositions.forEach(position -> {
            GeneralPosition.Builder newPosBuilder = new GeneralPosition.Builder(position);
            byte ref = GenotypeTable.UNKNOWN_ALLELE;
            try {
                ref = reference.genotype(position.getChromosome(), position.getPosition());
            } catch (Exception e) {
                // do nothing, remains unknown
            }
            newPosBuilder.allele(WHICH_ALLELE.Reference, ref);
            Position newPos = newPosBuilder.build();
            builder.add(newPos);
        });

        GenotypeTable result = GenotypeTableBuilder.getInstance(orig.genotypeMatrix(), builder.build(), orig.taxa(), orig.depth(), orig.alleleProbability(), orig.referenceProbability(), orig.dosage(), orig.annotations());

        return new DataSet(new Datum(temp.get(0).getName() + "_withReference", result, null), this);
    }

    /**
     * Input Reference Fasta
     *
     * @return Reference
     */
    public String reference() {
        return myReference.value();
    }

    /**
     * Set Reference. Input Reference Fasta
     *
     * @param value Reference
     *
     * @return this plugin
     */
    public AddReferenceToGenotypePlugin reference(String value) {
        myReference = new PluginParameter<>(myReference, value);
        return this;
    }

    @Override
    public ImageIcon getIcon() {
        return null;
    }

    @Override
    public String getButtonName() {
        return "Add Reference to Genotype";
    }

    @Override
    public String getToolTipText() {
        return "Add Reference to Genotype";
    }
}

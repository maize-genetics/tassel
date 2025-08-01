/*
 * ModifyTBTHDF5Plugin
 */
package net.maizegenetics.analysis.gbs;

import cern.colt.list.IntArrayList;

import net.maizegenetics.dna.tag.TagsByTaxaByteHDF5TaxaGroups;
import net.maizegenetics.dna.tag.TagsByTaxaByteHDF5TagGroups;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;

import java.awt.Frame;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.TreeMap;

import javax.swing.ImageIcon;
import net.maizegenetics.plugindef.PluginParameter;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 * This pipeline modifies TagsByTaxa HDF5 file with data organized by taxa. It
 * can: 1. Create an empty TBT. 2. Merge two TBT 3. Combined similarly named
 * taxa 4. Pivot a taxa TBT to a tag TBT
 *
 * @author ed
 */
public class ModifyTBTHDF5Plugin extends AbstractPlugin {

    private static final Logger myLogger = LogManager.getLogger(ModifyTBTHDF5Plugin.class);

    private PluginParameter<String> myTargetTBT = new PluginParameter.Builder<>("o", null, String.class).guiName("Target TBT HDF5 File").required(true).inFile()
            .description("Target TBT HDF5 (*tbt.h5) file to be modified. (Depending on the modification that you wish to make, choose only one of the only three parameters)").build();
    private PluginParameter<String> myAdditionalTaxaFile = new PluginParameter.Builder<>("i", null, String.class).guiName("Addition Taxa File").inFile()
            .description("TBT HDF5 (*tbt.h5) file containing additional taxa to be added to the target TBT HDF5 file").build();
    private PluginParameter<String> myTransposeFile = new PluginParameter.Builder<>("p", null, String.class).guiName("Pivot TBT HDF5 File").outFile()
            .description("Pivot (transpose) the target TBT HDF5 file into a tag-optimized orientation").build();
    private PluginParameter<Boolean> myCombineTaxa = new PluginParameter.Builder<>("c", false, Boolean.class).guiName("Merge Taxa")
            .description("Merge taxa in the target TBT HDF5 file with same LibraryPrepID").build();

    IntArrayList[] taxaReads;
    int[] readsPerSample, mappedReadsPerSample;
    int goodBarcodedReads = 0, allReads = 0, goodMatched = 0;
    HashMap<String, Integer> taxaNameToIndices;
    TagsByTaxaByteHDF5TaxaGroups targetTBT = null;

    public ModifyTBTHDF5Plugin() {
        super(null, false);
    }

    public ModifyTBTHDF5Plugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    @Override
    protected void postProcessParameters() {
        int count = 0;
        if (!myAdditionalTaxaFile.isEmpty()) {
            count++;
        }
        if (!myTransposeFile.isEmpty()) {
            count++;
        }
        if (mergeTaxa()) {
            count++;
        }
        if (count != 1) {
            throw new IllegalArgumentException("ModifyTBTHDF5Plugin: must specify exactly one paramter (Addition Taxa File, Pivot TBT HDF5 File, or Merge Taxa)");
        }
    }

    @Override
    public DataSet processData(DataSet input) {
        targetTBT = new TagsByTaxaByteHDF5TaxaGroups(targetTBTHDF5File());
        if ((additionTaxaFile() != null) && !additionTaxaFile().isEmpty()) {
            addAllTaxaToNewHDF5(additionTaxaFile());
        }
        if (mergeTaxa()) {
            combineTaxaHDF5();
        }
        if ((pivotTBTHDF5File() != null) && !pivotTBTHDF5File().isEmpty()) {
            TagsByTaxaByteHDF5TagGroups tranTBT = new TagsByTaxaByteHDF5TagGroups(targetTBT, pivotTBTHDF5File());
        }
        targetTBT.getFileReadyForClosing();
        targetTBT = null;
        System.gc();
        return null;
    }

    private boolean addAllTaxaToNewHDF5(String addTBTName) {
        TagsByTaxaByteHDF5TaxaGroups addTBT = new TagsByTaxaByteHDF5TaxaGroups(addTBTName);
        for (int i = 0; i < addTBT.getTagCount(); i++) {
            if (!Arrays.equals(targetTBT.getTag(i), addTBT.getTag(i))) {
                System.err.println("Tags are not the same for the two TBT file.  They cannot be combined.");
                System.exit(0);
            }
        }
        for (int i = 0; i < addTBT.getTaxaCount(); i++) {
            String name = addTBT.getTaxaName(i);
            byte[] states = addTBT.getReadCountDistributionForTaxon(i);
            targetTBT.addTaxon(name, states);
        }
        addTBT.getFileReadyForClosing();
        return true;
    }

    private boolean combineTaxaHDF5() {
        TreeMap<String, ArrayList<String>> combineTaxa = new TreeMap<String, ArrayList<String>>();
        for (String tn : targetTBT.getTaxaNames()) {
            String ptn = parseTaxaName(tn, "#");
            ArrayList<String> taxaList = combineTaxa.get(ptn);
            if (taxaList == null) {
                combineTaxa.put(ptn, taxaList = new ArrayList<String>());
            }
            taxaList.add(tn);
        }
        for (ArrayList<String> taxaList : combineTaxa.values()) {
            if (taxaList.size() < 2) {
                continue;
            }
            byte[] di = new byte[targetTBT.getTagCount()];
            String ptn = parseTaxaName(taxaList.get(0), "" + taxaList.size());
            for (int i = 0; i < taxaList.size(); i++) {
                int j = targetTBT.getIndexOfTaxaName(taxaList.get(i));
                byte[] dj = targetTBT.getReadCountDistributionForTaxon(j);
                for (int k = 0; k < dj.length; k++) {
                    int ts = di[k] + dj[k];
                    if (ts > Byte.MAX_VALUE) {
                        di[k] = Byte.MAX_VALUE;
                        // System.out.println(di[k]+"+"+dj[k]);
                    } else {
                        di[k] = (byte) ts;
                    }
                }
            }
            targetTBT.addTaxon(ptn, di);
            for (String tn : taxaList) {
                targetTBT.deleteTaxon(tn);
            }
        }
        return true;
    }

    private String parseTaxaName(String tn, String cnt) {
        String[] s = tn.split(":");
        return s[0] + ":MRG:" + cnt + ":" + s[3];
    }

    public String targetTBTHDF5File() {
        return myTargetTBT.value();
    }

    public ModifyTBTHDF5Plugin targetTBTHDF5File(String value) {
        myTargetTBT = new PluginParameter<>(myTargetTBT, value);
        return this;
    }

    public String additionTaxaFile() {
        return myAdditionalTaxaFile.value();
    }

    public ModifyTBTHDF5Plugin additionTaxaFile(String value) {
        myAdditionalTaxaFile = new PluginParameter<>(myAdditionalTaxaFile, value);
        return this;
    }

    public Boolean mergeTaxa() {
        return myCombineTaxa.value();
    }

    public ModifyTBTHDF5Plugin mergeTaxa(Boolean value) {
        myCombineTaxa = new PluginParameter<>(myCombineTaxa, value);
        return this;
    }

    public String pivotTBTHDF5File() {
        return myTransposeFile.value();
    }

    public ModifyTBTHDF5Plugin pivotTBTHDF5File(String value) {
        myTransposeFile = new PluginParameter<>(myTransposeFile, value);
        return this;
    }

    @Override
    public ImageIcon getIcon() {
        return null;
    }

    @Override
    public String getButtonName() {
        return "Modify TBT HDF5";
    }

    @Override
    public String getToolTipText() {
        return "Modify TBT HDF5";
    }

}

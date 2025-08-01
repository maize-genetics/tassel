/*
 */
package net.maizegenetics.analysis.association;

import java.util.ArrayList;
import org.junit.Test;
import net.maizegenetics.analysis.data.FileLoadPlugin;
import net.maizegenetics.analysis.data.FileLoadPlugin.TasselFileType;
import net.maizegenetics.analysis.distance.EndelmanDistanceMatrix;
import net.maizegenetics.constants.TutorialConstants;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.phenotype.Phenotype;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.taxa.distance.DistanceMatrix;
import static org.junit.Assert.assertTrue;

/**
 *
 * @author chd45
 */
public class GenomicSelectionPluginTest {

    @Test
    public void test1() {
        FileLoadPlugin flp = new FileLoadPlugin(null, false);
        flp.setTheFileType(TasselFileType.Hapmap);
        flp.setOpenFiles(new String[]{TutorialConstants.HAPMAP_FILENAME});
        DataSet ds1 = flp.performFunction(null);
        GenotypeTable geno = (GenotypeTable) ds1.getData(0).getData();
        DistanceMatrix kin = EndelmanDistanceMatrix.getInstance(geno);
        DataSet dsKin = new DataSet(new Datum("name", kin, ""), null);

        flp.setTheFileType(TasselFileType.Phenotype);
        flp.setOpenFiles(new String[]{TutorialConstants.TRAITS_FILENAME});
        DataSet ds2 = flp.performFunction(null);
        Phenotype pheno = (Phenotype) ds2.getData(0).getData();

        ArrayList<DataSet> mergedInputs = new ArrayList<>();
        mergedInputs.add(dsKin);
        mergedInputs.add(ds2);

        DataSet mergedData = DataSet.getDataSet(mergedInputs, flp);

        GenomicSelectionPlugin pluginTest = new GenomicSelectionPlugin(null, false);
        pluginTest.kFolds(5);
        pluginTest.nIterations(2);
        DataSet outData = pluginTest.processData(mergedData);

        System.out.println(outData);
        assertTrue(true);
    }

}

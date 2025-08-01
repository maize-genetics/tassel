package net.maizegenetics.analysis.modelfitter;

import static org.junit.Assert.*;

import java.util.List;
import java.util.Optional;

import net.maizegenetics.constants.GeneralConstants;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.ImportUtils;
import net.maizegenetics.phenotype.GenotypePhenotype;
import net.maizegenetics.phenotype.GenotypePhenotypeBuilder;
import net.maizegenetics.phenotype.Phenotype;
import net.maizegenetics.phenotype.PhenotypeBuilder;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.util.TableReport;
import net.maizegenetics.util.TableReportTestUtils;
import net.maizegenetics.util.TableReportUtils;

import org.junit.Ignore;
import org.junit.Test;

public class StepwiseAdditiveModelFitterTest {
    @Ignore
    @Test
    public void test() {
        StepwiseAdditiveModelFitterPlugin samfPlugin = new StepwiseAdditiveModelFitterPlugin();
        samfPlugin.maxTerms(10);
        samfPlugin.enterLimit(1e-2);
        samfPlugin.exitLimit(2e-2);
        samfPlugin.isNested(false);
        
        //to create files for validation
//        samfPlugin.writeFiles(true);
//        samfPlugin.outputName(GeneralConstants.EXPECTED_RESULTS_DIR + "stepwise_additive_model_result");
        
        GenotypeTable geno = ImportUtils.readFromHapmap(GeneralConstants.DATA_DIR + "CandidateTests/mdp_genotype_chr9_10_nomissing.hmp.txt");
        Phenotype pheno = new PhenotypeBuilder().fromFile(GeneralConstants.DATA_DIR + "CandidateTests/mdp_traits_nomissing.txt")
                .keepAttributes(new int[]{0,2})
                .build().get(0);
        GenotypePhenotype genoPheno = new GenotypePhenotypeBuilder().genotype(geno).phenotype(pheno).intersect().build();
        
        DataSet inputData = new DataSet(new Datum("name",  genoPheno, "comment"), null);
        DataSet outData = samfPlugin.performFunction(inputData);
        
        //expected results
        TableReport expectedAnova = TableReportUtils.readDelimitedTableReport(GeneralConstants.EXPECTED_RESULTS_DIR + "stepwise_additive_model_result_anova.txt", "\t");
        TableReport expectedEffects = TableReportUtils.readDelimitedTableReport(GeneralConstants.EXPECTED_RESULTS_DIR + "stepwise_additive_model_result_effects.txt", "\t");
        TableReport expectedPrescanAnova = TableReportUtils.readDelimitedTableReport(GeneralConstants.EXPECTED_RESULTS_DIR + "stepwise_additive_model_result_prescan_anova.txt", "\t");
        TableReport expectedSteps = TableReportUtils.readDelimitedTableReport(GeneralConstants.EXPECTED_RESULTS_DIR + "stepwise_additive_model_result_steps.txt", "\t");
        
        List<Datum> reportList = outData.getDataOfType(TableReport.class);
        
        Optional<Datum> reportDatum = reportList.stream().filter(datum -> datum.getName().startsWith("Anova_")).findFirst();
        assertTrue(reportDatum.isPresent());
        TableReport obsAnova = (TableReport) reportDatum.get().getData();
        
        reportDatum = reportList.stream().filter(datum -> datum.getName().startsWith("Effects_regression_")).findFirst();
        assertTrue(reportDatum.isPresent());
        TableReport obsEffects = (TableReport) reportDatum.get().getData();
        
        reportDatum = reportList.stream().filter(datum -> datum.getName().startsWith("Prescan_anova_")).findFirst();
        assertTrue(reportDatum.isPresent());
        TableReport obsPrescanAnova = (TableReport) reportDatum.get().getData();
        
        reportDatum = reportList.stream().filter(datum -> datum.getName().startsWith("Steps_")).findFirst();
        assertTrue(reportDatum.isPresent());
        TableReport obsSteps = (TableReport) reportDatum.get().getData();
        
        TableReportTestUtils.compareTableReportValues(expectedAnova, obsAnova, 1e-4);
        TableReportTestUtils.compareTableReportValues(expectedEffects, obsEffects, 1e-4);
        TableReportTestUtils.compareTableReportValues(expectedPrescanAnova, obsPrescanAnova, 1e-4);
        TableReportTestUtils.compareTableReportValues(expectedSteps, obsSteps, 1e-4);
        
        
    }

}

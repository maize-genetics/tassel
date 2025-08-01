/*
 * BasicTassel
 */
package net.maizegenetics.testsuites;

import net.maizegenetics.analysis.association.DiscreteSitesTest;
import net.maizegenetics.analysis.association.EqtlAssociationPluginTest;
import net.maizegenetics.analysis.association.MLMTest;
import net.maizegenetics.analysis.association.PhenotypeLMTest;
import net.maizegenetics.analysis.association.ReferenceProbabilityFELMTest;
import net.maizegenetics.analysis.data.GenosToABHPluginTest;
import net.maizegenetics.analysis.data.GenotypeSummaryPluginTest;
import net.maizegenetics.analysis.data.MigrateHDF5FromT4T5Test;
import net.maizegenetics.analysis.data.PlinkLoadPluginTest;
import net.maizegenetics.analysis.distance.AMatrixPluginTest;
import net.maizegenetics.analysis.distance.AddDistanceMatrixPluginTest;
import net.maizegenetics.analysis.distance.CenteredIBSTest;
import net.maizegenetics.analysis.distance.DominanceCenteredIBSTest;
import net.maizegenetics.analysis.distance.DominanceNormalizedIBSTest;
import net.maizegenetics.analysis.distance.IBSDistanceMatrixTest;
import net.maizegenetics.analysis.distance.KinshipTest;
import net.maizegenetics.analysis.distance.NormalizedIBSTest;
import net.maizegenetics.analysis.distance.RemoveNaNFromDistanceMatrixPluginTest;
import net.maizegenetics.analysis.distance.SubtractDistanceMatrixPluginTest;
import net.maizegenetics.analysis.distance.WriteDistanceMatrixTest;
import net.maizegenetics.analysis.filter.FilterAlignmentPluginTest;
import net.maizegenetics.analysis.filter.FilterSiteBuilderPluginTest;
import net.maizegenetics.analysis.imputation.FILLINImputationPluginTest;
import net.maizegenetics.analysis.imputation.FullSibFamilyHaplotypeTest;
import net.maizegenetics.analysis.imputation.RandomGenotypeImputationTest;
import net.maizegenetics.analysis.numericaltransform.ImputationByMeanTest;
import net.maizegenetics.analysis.numericaltransform.TransformDataPluginTest;
import net.maizegenetics.analysis.numericaltransform.kNearestNeighborsTest;
import net.maizegenetics.baseplugins.DistanceMatrixPluginTest;
import net.maizegenetics.dna.map.ChromosomeTest;
import net.maizegenetics.dna.snp.AlignmentScopeTest;
import net.maizegenetics.dna.snp.ExportUtilsTest;
import net.maizegenetics.dna.snp.FilterGenotypeTableTest;
import net.maizegenetics.dna.snp.NucleotideAlignmentConstantsTest;
import net.maizegenetics.dna.snp.genotypecall.GenotypeTest;
import net.maizegenetics.dna.map.PositionHDF5ListTest;
import net.maizegenetics.dna.snp.AlignmentBuilderTest;
import net.maizegenetics.dna.snp.CombineGenotypeTableTest;
import net.maizegenetics.dna.snp.GenotypeTableBuilderTest;
import net.maizegenetics.dna.snp.GenotypeTableStreamTest;
import net.maizegenetics.dna.snp.ImportUtilsTest;
import net.maizegenetics.dna.snp.LineIndexHapmapTest;
import net.maizegenetics.dna.snp.genotypecall.AlleleFreqCacheTest;
import net.maizegenetics.dna.snp.score.AlleleDepthUtilTest;
import net.maizegenetics.dna.snp.genotypecall.BasicGenotypeMergeRuleTest;
import net.maizegenetics.dna.snp.io.BuilderFromVCFTest;
import net.maizegenetics.popgen.LinkageDisequilibriumTest;
import net.maizegenetics.prefs.TasselPrefsTest;
import net.maizegenetics.stats.linearmodels.ModelEffectTest;
import net.maizegenetics.taxa.IdentifierSynonymizerTest;
import net.maizegenetics.taxa.TaxaListBuilderTest;
import net.maizegenetics.taxa.TaxaListIOUtilsTest;
import net.maizegenetics.util.BitSetTest;
import net.maizegenetics.util.CheckSumTest;
import net.maizegenetics.util.LoggingUtils;
import net.maizegenetics.util.SuperByteMatrixTest;

import org.junit.runner.JUnitCore;
import org.junit.runner.Result;

/**
 *
 * @author Terry Casstevens
 */
public class BasicTassel {

    public static void main(String[] args) {

        LoggingUtils.setupDebugLogging();

        Result result = JUnitCore.runClasses(
                TasselPrefsTest.class,
                SetupPreferences.class,
                ExportUtilsTest.class,
                ImportUtilsTest.class,
                FilterGenotypeTableTest.class,
                CombineGenotypeTableTest.class,
                FilterAlignmentPluginTest.class,
                PlinkLoadPluginTest.class,
                DistanceMatrixPluginTest.class,
                IBSDistanceMatrixTest.class,
                LinkageDisequilibriumTest.class,
                BitSetTest.class,
                AlignmentScopeTest.class,
                DiscreteSitesTest.class,
                PhenotypeLMTest.class,
                ReferenceProbabilityFELMTest.class,
                GenotypeTest.class,
                SuperByteMatrixTest.class,
                FILLINImputationPluginTest.class,
                NucleotideAlignmentConstantsTest.class,
                CheckSumTest.class,
                AlignmentBuilderTest.class,
                GenotypeSummaryPluginTest.class,
                PositionHDF5ListTest.class,
                TaxaListBuilderTest.class,
                TaxaListIOUtilsTest.class,
                BasicGenotypeMergeRuleTest.class,
                FullSibFamilyHaplotypeTest.class,
                RandomGenotypeImputationTest.class,
                MigrateHDF5FromT4T5Test.class,
                net.maizegenetics.dna.snp.depth.AlleleDepthUtilTest.class,
                AlleleDepthUtilTest.class,
                GenotypeTableBuilderTest.class,
                KinshipTest.class,
                ImputationByMeanTest.class,
                kNearestNeighborsTest.class,
                TransformDataPluginTest.class,
                EqtlAssociationPluginTest.class,
                MLMTest.class,
                IdentifierSynonymizerTest.class,
                GenotypeTableStreamTest.class,
                GenosToABHPluginTest.class,
                ChromosomeTest.class,
                BuilderFromVCFTest.class,
                ModelEffectTest.class,
                AlleleFreqCacheTest.class,
                FilterSiteBuilderPluginTest.class,
                AMatrixPluginTest.class,
                LineIndexHapmapTest.class,
                DominanceCenteredIBSTest.class,
                DominanceNormalizedIBSTest.class,
                CenteredIBSTest.class,
                NormalizedIBSTest.class,
                SubtractDistanceMatrixPluginTest.class,
                AddDistanceMatrixPluginTest.class,
                WriteDistanceMatrixTest.class,
                RemoveNaNFromDistanceMatrixPluginTest.class);

        PrintResults.printResults(result);
    }
}

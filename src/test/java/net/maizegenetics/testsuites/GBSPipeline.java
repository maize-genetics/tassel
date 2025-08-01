/*
 * GBSPipeline
 */
package net.maizegenetics.testsuites;

import net.maizegenetics.analysis.gbs.ModifyTBTHDF5PluginTest;
import net.maizegenetics.analysis.gbs.FastqToTagCountPluginTest;
import net.maizegenetics.analysis.gbs.TagCountToFastqPluginTest;
import net.maizegenetics.analysis.gbs.SAMConverterPluginTest;
import net.maizegenetics.analysis.gbs.DiscoverySNPCallerPluginTest;
import net.maizegenetics.analysis.gbs.SeqToTBTHDF5PluginTest;
import net.maizegenetics.analysis.gbs.MergeMultipleTagCountPluginTest;
import net.maizegenetics.util.LoggingUtils;

import org.junit.runner.JUnitCore;
import org.junit.runner.Result;

/**
 *
 * @author terry
 */
public class GBSPipeline {

    public static void main(String[] args) {

        LoggingUtils.setupLogging();

        Result result = JUnitCore.runClasses(FastqToTagCountPluginTest.class, MergeMultipleTagCountPluginTest.class,
                TagCountToFastqPluginTest.class, SAMConverterPluginTest.class, SeqToTBTHDF5PluginTest.class,
                ModifyTBTHDF5PluginTest.class, DiscoverySNPCallerPluginTest.class);

        PrintResults.printResults(result);

    }
}

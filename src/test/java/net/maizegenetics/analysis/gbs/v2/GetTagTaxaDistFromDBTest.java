/**
 * 
 */
package net.maizegenetics.analysis.gbs.v2;

import org.junit.Test;

import net.maizegenetics.constants.GBSConstants;

/**
 * @author lcj34
 *
 */
public class GetTagTaxaDistFromDBTest {
    private final static String inputDB = GBSConstants.GBS_GBS2DB_FILE;
    private final static String outputFile=GBSConstants.GBS_TEMP_DIR + "tagsTaxaDepthOutput.txt";
    @Test
    public  void testTagTaxaDistFromDBPlugin() {
        new GetTagTaxaDistFromDBPlugin()
        .inputDB(inputDB)
        .outputFile(outputFile)
        .performFunction(null);
    }
}

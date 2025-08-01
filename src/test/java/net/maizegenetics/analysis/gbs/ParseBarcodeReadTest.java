/*
 * ParseBarcodeReadTest
 */
package net.maizegenetics.analysis.gbs;

import net.maizegenetics.constants.GBSConstants;
import net.maizegenetics.util.Utils;
import org.junit.Test;

/**
 *
 * @author terry
 */
public class ParseBarcodeReadTest {

    public ParseBarcodeReadTest() {
    }

    @Test
    public void testParseBarcodeRead() {

        ParseBarcodeRead pbr = getParseBarcodeRead(GBSConstants.GBS_INPUT_FASTQ_FILE, GBSConstants.GBS_TESTING_KEY_FILE, "ApeKI");

    }

    private ParseBarcodeRead getParseBarcodeRead(String fastqFile, String keyFilename, String enzyme) {

        String[] np = Utils.getFilename(fastqFile).split("_");

        String flowcell;
        String lane;
        if (np.length == 3) {
            flowcell = np[0];
            lane = np[1];
        } else if (np.length == 5) {
            flowcell = np[1];
            lane = np[3];
        } else if (np.length == 4) {
            flowcell = np[0];
            lane = np[2];
        } else if (np.length == 6) {
            flowcell = np[1];
            lane = np[3];
        } else {
            System.out.println("Error in parsing file name:");
            System.out.println("   The filename does not contain either 3 or 5 underscore-delimited values.");
            System.out.println("   Expect: flowcell_lane_fastq.txt OR code_flowcell_s_lane_fastq.txt");
            System.out.println("   Filename: " + fastqFile);
            return null;
        }

        ParseBarcodeRead result = new ParseBarcodeRead(keyFilename, enzyme, flowcell, lane);
        System.out.println("Fastq File: " + fastqFile + "  Key File: " + keyFilename + "  Enzyme: " + enzyme + "  Flowcell: " + flowcell + "  Lane: " + lane + "Total Barcodes Found: " + result.getBarCodeCount());
        return result;

    }
}
/*
 * CheckSumTest
 */
package net.maizegenetics.util;

import net.maizegenetics.constants.TutorialConstants;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 * Testing CheckSum with trivial test of hapmap file.
 * @author Dallas Kroon
 */
public class CheckSumTest {


    @Test
    public void testGetChecksum(){
        System.out.println("Testing CheckSum module");

        String protocol = "MD5";
        String knownMd5sum = "6309cda5e513b44cbbab8d3d9a77c364";

        String aChecksum = CheckSum.getChecksum(TutorialConstants.HAPMAP_FILENAME, protocol);

        assertEquals(knownMd5sum, aChecksum);
        
        System.out.println("Testing MD5 getProtocolChecksum for " + TutorialConstants.HAPMAP_FILENAME);
        
        aChecksum = CheckSum.getProtocolChecksum(TutorialConstants.HAPMAP_FILENAME, protocol);
        assertEquals(knownMd5sum, aChecksum);
        
        String knownSHA1sum = "8e6f210ecf50ec451d52d9c55a4b34a293442847";
        protocol = "SHA-1";
        System.out.println("Testing SHA-1 getProtocolChecksum for " + TutorialConstants.HAPMAP_FILENAME);
        aChecksum = CheckSum.getProtocolChecksum(TutorialConstants.HAPMAP_FILENAME, protocol);
        assertEquals(knownSHA1sum, aChecksum);
        
        System.out.println("\nFinished testGetChecksum");
    }
}
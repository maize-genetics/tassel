/*
 * TasselPrefsTest
 */
package net.maizegenetics.prefs;

import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author terry
 */
public class TasselPrefsTest {

    public TasselPrefsTest() {
    }

    /**
     * Test of setPersistPreferences method, of class TasselPrefs.
     */
    @Test
    public void testSetPersistPreferences() {

        System.out.println("Preferences...");
        TasselPrefs.setPersistPreferences(false);

        boolean retainRareAlleles = TasselPrefs.getAlignmentRetainRareAlleles();
        double maxFreq = TasselPrefs.getFilterAlignPluginMaxFreq();
        int minCount = TasselPrefs.getFilterAlignPluginMinCount();
        double minFreq = TasselPrefs.getFilterAlignPluginMinFreq();
        double maxHetFreq = TasselPrefs.getFilterTaxaPropsMaxHetFreq();
        double minHetFreq = TasselPrefs.getFilterTaxaPropsMinHetFreq();
        double minNotMissingFreq = TasselPrefs.getFilterTaxaPropsMinNotMissingFreq();
        String openDir = TasselPrefs.getOpenDir();
        String saveDir = TasselPrefs.getSaveDir();

        assertEquals("Retain Rare Alleles Default: ", TasselPrefs.ALIGNMENT_RETAIN_RARE_ALLELES_DEFAULT, retainRareAlleles);
        assertEquals("Filter Alignment Max Freq Default: ", TasselPrefs.FILTER_ALIGN_PLUGIN_MAX_FREQ_DEFAULT, maxFreq, 0.0);
        assertEquals("Filter Alignment Min Count Default: ", TasselPrefs.FILTER_ALIGN_PLUGIN_MIN_COUNT_DEFAULT, minCount);
        assertEquals("Filter Alignment Min Freq Default: ", TasselPrefs.FILTER_ALIGN_PLUGIN_MIN_FREQ_DEFAULT, minFreq, 0.0);
        assertEquals("Filter Taxa Mex Het Freq Default: ", TasselPrefs.FILTER_TAXA_PROPS_PLUGIN_MAX_HET_DEFAULT, maxHetFreq, 0.0);
        assertEquals("Filter Taxa Min Het Freq Default: ", TasselPrefs.FILTER_TAXA_PROPS_PLUGIN_MIN_HET_DEFAULT, minHetFreq, 0.0);
        assertEquals("Filter Taxa Not Missing Freq Default: ", TasselPrefs.FILTER_TAXA_PROPS_PLUGIN_MIN_NOT_MISSING_DEFAULT, minNotMissingFreq, 0.0);
        assertEquals("Open Dir Default: ", TasselPrefs.TASSEL_OPEN_DIR_DEFAULT, openDir);
        assertEquals("Save Dir Default: ", TasselPrefs.TASSEL_SAVE_DIR_DEFAULT, saveDir);

        System.out.println("Retain Rare Alleles Default: " + TasselPrefs.ALIGNMENT_RETAIN_RARE_ALLELES_DEFAULT);
        System.out.println("Filter Alignment Max Freq Default: " + TasselPrefs.FILTER_ALIGN_PLUGIN_MAX_FREQ_DEFAULT);
        System.out.println("Filter Alignment Min Count Default: " + TasselPrefs.FILTER_ALIGN_PLUGIN_MIN_COUNT_DEFAULT);
        System.out.println("Filter Alignment Min Freq Default: " + TasselPrefs.FILTER_ALIGN_PLUGIN_MIN_FREQ_DEFAULT);
        System.out.println("Filter Taxa Mex Het Freq Default: " + TasselPrefs.FILTER_TAXA_PROPS_PLUGIN_MAX_HET_DEFAULT);
        System.out.println("Filter Taxa Min Het Freq Default: " + TasselPrefs.FILTER_TAXA_PROPS_PLUGIN_MIN_HET_DEFAULT);
        System.out.println("Filter Taxa Not Missing Freq Default: " + TasselPrefs.FILTER_TAXA_PROPS_PLUGIN_MIN_NOT_MISSING_DEFAULT);
        System.out.println("Open Dir Default: " + TasselPrefs.TASSEL_OPEN_DIR_DEFAULT);
        System.out.println("Save Dir Default: " + TasselPrefs.TASSEL_SAVE_DIR_DEFAULT);
    }
}
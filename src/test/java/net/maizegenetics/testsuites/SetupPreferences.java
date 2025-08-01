/*
 * SetupPreferences
 */
package net.maizegenetics.testsuites;

import net.maizegenetics.prefs.TasselPrefs;
import org.junit.Test;

/**
 *
 * @author terry
 */
public class SetupPreferences {

    public SetupPreferences() {
    }

    /**
     * Set up Preferences for Out Tests
     */
    @Test
    public void testSetupPreferences() {
        TasselPrefs.putAlignmentRetainRareAlleles(true);
    }
}
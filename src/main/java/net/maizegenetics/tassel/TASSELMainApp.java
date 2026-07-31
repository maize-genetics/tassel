/*
 * TASSEL - Trait Analysis by a aSSociation Evolution & Linkage
 * Copyright (C) 2003 Ed Buckler
 *
 * This software evaluates linkage disequilibrium nucletide diversity and
 * associations. For more information visit http://www.maizegenetics.net
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
// Title:      TASSELMainApp
// Version:
// Copyright:  Copyright (c) 1998
// Author:     Ed Buckler
package net.maizegenetics.tassel;

import net.maizegenetics.gui.DialogUtils;
import net.maizegenetics.pipeline.TasselPipeline;
import net.maizegenetics.plugindef.ParameterCache;
import net.maizegenetics.prefs.TasselPrefs;
import net.maizegenetics.util.LoggingUtils;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class TASSELMainApp {

    private static final Logger myLogger = LogManager.getLogger(TASSELMainApp.class);

    private TASSELMainApp() {
    }

    /**
     * Installs the FlatLaf Look-and-Feel, choosing the light or dark theme based on the
     * persisted user preference. Safe to call again at runtime (e.g. from the Preferences
     * dialog) followed by {@code FlatLaf.updateUI()} to switch themes live.
     */
    public static void setupLookAndFeel() {
        try {
            if (TasselPrefs.getDarkTheme()) {
                com.formdev.flatlaf.FlatDarkLaf.setup();
            } else {
                com.formdev.flatlaf.FlatLightLaf.setup();
            }
        } catch (Exception e) {
            myLogger.warn("Could not install FlatLaf look and feel; using default.", e);
        }
    }

    public static void main(String[] args) {
        TASSELMainFrame frame = null;
        try {

            // Native macOS integration (must be set before the AWT toolkit initializes).
            System.setProperty("apple.laf.useScreenMenuBar", "true");
            System.setProperty("apple.awt.application.name", "TASSEL");

            TasselPrefs.setPersistPreferences(true);
            LoggingUtils.setupLogging();

            setupLookAndFeel();

            frame = new TASSELMainFrame();
            frame.validate();
            frame.setLocationRelativeTo(null);
            frame.setVisible(true);

            ParameterCache.load(TasselPrefs.getConfigFile());

            if (args.length > 0) {
                new TasselPipeline(args, frame);
            }

        } catch (OutOfMemoryError e) {
            e.printStackTrace();
            StringBuilder builder = new StringBuilder();
            builder.append("Out of Memory: \n");
            long heapMaxSize = Runtime.getRuntime().maxMemory() / 1048576l;
            builder.append("Current Max Heap Size: ");
            builder.append(heapMaxSize);
            builder.append(" Mb\n");
            builder.append("Use -Xmx option in start_tassel.pl or start_tassel.bat\n");
            builder.append("to increase heap size.");
            builder.append(" Included with tassel standalone zip.");
            myLogger.error(builder.toString());
        } catch (Exception e) {
            myLogger.error(e.getMessage(), e);
            DialogUtils.showError(e.getMessage() + "\n", frame);
        }
    }
}

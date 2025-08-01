/*
 *  FileBrowserUtils
 */
package net.maizegenetics.gui;

import java.io.File;

import javax.swing.JDialog;
import javax.swing.JFileChooser;

import net.maizegenetics.prefs.TasselPrefs;
import net.maizegenetics.util.Utils;

/**
 *
 * @author Terry Casstevens
 */
public class FileBrowserUtils {

    private FileBrowserUtils() {
        // utilities
    }

    public static File getOpenFile(JDialog parent) {

        final JFileChooser fileChooser = new JFileChooser(TasselPrefs.getOpenDir());

        if (fileChooser.showOpenDialog(parent) == JFileChooser.APPROVE_OPTION) {
            TasselPrefs.putOpenDir(fileChooser.getCurrentDirectory().getPath());
            return fileChooser.getSelectedFile();
        } else {
            return null;
        }

    }

    public static File getSaveFile(JDialog parent) {

        final JFileChooser fileChooser = new JFileChooser(TasselPrefs.getSaveDir());

        if (fileChooser.showSaveDialog(parent) == JFileChooser.APPROVE_OPTION) {
            TasselPrefs.putSaveDir(fileChooser.getCurrentDirectory().getPath());
            return fileChooser.getSelectedFile();
        } else {
            return null;
        }

    }

    public static File getOpenDir(JDialog parent) {

        final JFileChooser fileChooser = new JFileChooser(Utils.getDirectory(TasselPrefs.getOpenDir()));
        fileChooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);

        if (fileChooser.showOpenDialog(parent) == JFileChooser.APPROVE_OPTION) {
            File file = fileChooser.getSelectedFile();
            TasselPrefs.putOpenDir(file.getPath());
            return file;
        } else {
            return null;
        }

    }

    public static File getSaveDir(JDialog parent) {

        final JFileChooser fileChooser = new JFileChooser(Utils.getDirectory(TasselPrefs.getSaveDir()));
        fileChooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);

        if (fileChooser.showOpenDialog(parent) == JFileChooser.APPROVE_OPTION) {
            File file = fileChooser.getSelectedFile();
            TasselPrefs.putSaveDir(file.getPath());
            return file;
        } else {
            return null;
        }

    }
}

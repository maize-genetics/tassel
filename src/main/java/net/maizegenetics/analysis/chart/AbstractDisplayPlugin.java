/*
 * AbstractDisplayPlugin.java
 *
 * Created on December 22, 2006, 5:02 PM
 *
 */
package net.maizegenetics.analysis.chart;

import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.prefs.TasselPrefs;
import net.maizegenetics.util.BasicFileFilter;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.graphics2d.svg.SVGGraphics2D;
import org.jfree.graphics2d.svg.SVGUtils;

import javax.imageio.ImageIO;
import javax.swing.*;
import javax.swing.filechooser.FileFilter;
import java.awt.*;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import java.io.File;

/**
 * This plugin abstract plugin provides mechanisms for setting a savefile and for saving a variety of images.
 * This would be useful to plugin that create images and may want to automatically send them to a file.
 *
 * @author Ed Buckler
 */
public abstract class AbstractDisplayPlugin extends AbstractPlugin {

    private File mySaveFile;
    private int myImageWidth = 500;
    private int myImageHeight = 500;

    public enum Outformat {

        svg, jpg, gif, bmp, wbmp, png, printer
    }

    private Outformat myOutformat = Outformat.svg;

    public AbstractDisplayPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    protected void saveDataToFile(BufferedImage img, Outformat format, File saveFile) {
        if (saveFile == null) {
            return;
        }
        try {
            switch (format) {
                case jpg: {
                    ImageIO.write(img, "JPEG", saveFile);
                    break;
                }
                case gif: {
                    ImageIO.write(img, "GIF", saveFile);
                    break;
                }
                case bmp: {
                    ImageIO.write(img, "BMP", saveFile);
                    break;
                }
                case wbmp: {
                    ImageIO.write(img, "WBMP", saveFile);
                    break;
                }
                case png: {
                    ImageIO.write(img, "PNG", saveFile);
                    break;
                }
                case printer: {
                    break;
                }
            }
        } catch (Exception ee) {
            System.out.println("saveDataToFile:" + ee);
        }
    }

    public void saveDataToFile(Component comp) {
        saveDataToFile(comp, getSaveFileByChooser());
    }

    public void saveDataToFile(Component comp, File saveFile) {
        try {
            int d_width = comp.getWidth();
            int d_height = comp.getHeight();
            if (d_width == 0) {
                d_width = myImageWidth;
                d_height = myImageHeight;
                comp.setSize(d_width, d_height);
            }
            if (myOutformat == Outformat.svg) {
                if (comp instanceof ChartPanel) {
                    JFreeChart chart = ((ChartPanel) comp).getChart();
                    Rectangle2D plotArea = ((ChartPanel) comp).getChartRenderingInfo().getPlotInfo().getDataArea();
                    int width = (int) plotArea.getWidth();
                    int height = (int) plotArea.getHeight();
                    SVGGraphics2D g2 = new SVGGraphics2D(width, height);
                    chart.draw(g2, new Rectangle(0, 0, width, height));
                    SVGUtils.writeToSVG(saveFile, g2.getSVGElement());
                } else {
                    SVGGraphics2D g2 = new SVGGraphics2D(d_width, d_height);
                    comp.paint(g2);
                    SVGUtils.writeToSVG(saveFile, g2.getSVGElement());
                }
            } else {
                BufferedImage img = new BufferedImage(d_width, d_height, BufferedImage.TYPE_INT_RGB);
                Graphics gbox = img.getGraphics();
                comp.paint(gbox);
                saveDataToFile(img, myOutformat, saveFile);
            }
        } catch (Exception ee) {
            throw new IllegalStateException("AbstractDisplayPlugin: saveDataToFile: problem writing: " + saveFile.getName() + " as: " + myOutformat + "\n" + ee.getMessage());
        }
    }

    /**
     * Provides a save filer that remember the last location something was opened from
     */
    private File getSaveFileByChooser() {

        File tempFile = new File(TasselPrefs.getSaveDir());
        JFileChooser filerSave = new JFileChooser(tempFile);
        filerSave.resetChoosableFileFilters();
        // add suffix filters to chooser
        filerSave.removeChoosableFileFilter(null); //remove the all
        for (Outformat p : Outformat.values()) {
            filerSave.addChoosableFileFilter(new BasicFileFilter(p.toString()));
        }
        javax.swing.filechooser.FileFilter[] ff = filerSave.getChoosableFileFilters();
        filerSave.setFileFilter(ff[1]);
        filerSave.setMultiSelectionEnabled(false);
        File saveFile = null;
        int returnVal = filerSave.showSaveDialog(null);
        if (returnVal == JFileChooser.SAVE_DIALOG || returnVal == JFileChooser.APPROVE_OPTION) {
            saveFile = filerSave.getSelectedFile();
            TasselPrefs.putSaveDir(filerSave.getCurrentDirectory().getPath());
            String fileType = ((BasicFileFilter) filerSave.getFileFilter()).getExtension();
            if (!filerSave.isAcceptAllFileFilterUsed()) {
                fileType = ((BasicFileFilter) filerSave.getFileFilter()).getExtension();
            }
            myOutformat = Outformat.valueOf(fileType);
            if (!saveFile.getName().contains(fileType)) {
                saveFile = new File(saveFile.getParent(), saveFile.getName() + "." + fileType);
            }
        }
        return saveFile;

    }

    protected File getSaveFileByChooser(String[] fileExtensions) {
        return getSaveFileByChooser(fileExtensions, getParentFrame());
    }

    protected File getSaveFileByChooser(String[] fileExtensions, Component parent) {

        JFileChooser fileChooser = new JFileChooser(new File(TasselPrefs.getSaveDir()));

        for (int i = 0; i < fileExtensions.length; i++) {
            FileFilter current = new BasicFileFilter(fileExtensions[i]);
            fileChooser.addChoosableFileFilter(current);
            if (i == 0) {
                fileChooser.setFileFilter(current);
            }
        }

        fileChooser.setMultiSelectionEnabled(false);

        File saveFile = null;
        int returnVal = fileChooser.showSaveDialog(parent);
        if (returnVal == JFileChooser.SAVE_DIALOG || returnVal == JFileChooser.APPROVE_OPTION) {
            saveFile = fileChooser.getSelectedFile();
            TasselPrefs.putSaveDir(fileChooser.getCurrentDirectory().getPath());

            try {
                String fileType = ((BasicFileFilter) fileChooser.getFileFilter()).getExtension();
                if (!saveFile.getName().endsWith("." + fileType)) {
                    saveFile = new File(saveFile.getParent(), saveFile.getName() + "." + fileType);
                }
            } catch (Exception e) {
                //do nothing
            }

        }

        return saveFile;

    }

    public static String[] getPossibleGraphicOutFormats() {

        String[] s = new String[Outformat.values().length];
        int i = 0;
        for (Outformat p : Outformat.values()) {
            s[i] = "" + p;
            i++;
        }
        return s;

    }

    public int getImageWidth() {
        return myImageWidth;
    }

    public int getImageHeight() {
        return myImageHeight;
    }

    public void setImageSize(int width, int height) {
        myImageWidth = width;
        myImageHeight = height;
    }

    public Outformat getOutformat() {
        return myOutformat;
    }

    public void setOutformat(Outformat theOutformat) {
        myOutformat = theOutformat;
    }

    public File getSaveFile() {
        return mySaveFile;
    }

    public void setSaveFile(File theSaveFile) {
        mySaveFile = theSaveFile;
    }

    public void setSaveFile(String theSaveFile) {
        mySaveFile = new File(theSaveFile);
    }
}

/*
 *  WorkflowPlugin
 * 
 *  Created on March 20, 2015
 */
package net.maizegenetics.analysis.workflow;

import java.io.File;
import java.net.URL;
import java.net.URLDecoder;

import java.util.ArrayList;
import java.util.Enumeration;
import java.util.List;
import java.util.jar.JarEntry;
import java.util.jar.JarFile;

import javax.swing.ImageIcon;

import net.maizegenetics.pipeline.TasselPipeline;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.tassel.TASSELMainFrame;
import net.maizegenetics.util.Utils;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 *
 * @author Terry Casstevens
 */
public class WorkflowPlugin extends AbstractPlugin {

    private static final Logger myLogger = LogManager.getLogger(WorkflowPlugin.class);

    private final String myButtonName;
    private final String[] myArgs;
    private final TASSELMainFrame myFrame;

    private WorkflowPlugin(TASSELMainFrame parentFrame, String buttonName, String[] args) {
        super(parentFrame, true);
        myButtonName = buttonName;
        myArgs = args;
        myFrame = parentFrame;
    }

    public static List<WorkflowPlugin> getInstances(TASSELMainFrame parentFrame) {
        List<WorkflowPlugin> result = new ArrayList<>();
        List<String> filenames = getConfigFiles();
        if (!filenames.isEmpty()) {
            for (String filename : filenames) {
                String[] args = new String[]{"-configResourceFile", filename};
                WorkflowPlugin test = new WorkflowPlugin(parentFrame, buttonName(filename), args);
                result.add(test);
            }
        } else {
            filenames = getConfigFilesFromSrc();
            for (String filename : filenames) {
                String[] args = new String[]{"-configFile", filename};
                WorkflowPlugin test = new WorkflowPlugin(parentFrame, buttonName(filename), args);
                result.add(test);
            }
        }
        return result;
    }

    private static List<String> getConfigFiles() {

        List<String> result = new ArrayList<>();

        try {
            URL directory = WorkflowPlugin.class.getResource("/net/maizegenetics/analysis/workflow/");
            String directoryPath = directory.getPath();
            if (directoryPath.contains("!")) {
                String jarPath = directory.getPath().substring(5, directoryPath.indexOf("!"));
                JarFile jar = new JarFile(URLDecoder.decode(jarPath, "UTF-8"));
                Enumeration<JarEntry> entries = jar.entries();
                while (entries.hasMoreElements()) {
                    String name = entries.nextElement().getName();
                    if ((name.startsWith("net/maizegenetics/analysis/workflow/")) && (name.endsWith(".xml"))) {
                        result.add("/" + name);
                    }
                }
            }
        } catch (Exception e) {
            myLogger.debug(e.getMessage(), e);
        }

        return result;

    }

    private static List<String> getConfigFilesFromSrc() {

        List<String> result = new ArrayList<>();

        try {
            URL directory = WorkflowPlugin.class.getResource("/net/maizegenetics/analysis/workflow/");
            String directoryPath = directory.getPath();
            File temp = new File(directoryPath);
            for (File current : temp.listFiles()) {
                String name = current.getCanonicalPath();
                if (name.endsWith(".xml")) {
                    result.add("/" + name);
                }
            }
        } catch (Exception e) {
            myLogger.debug(e.getMessage(), e);
        }

        return result;

    }

    private static String buttonName(String filename) {
        String temp = Utils.getFilename(filename, ".xml");
        StringBuilder builder = new StringBuilder();
        builder.append(Character.toUpperCase(temp.charAt(0)));
        boolean newWord = false;
        for (int i = 1; i < temp.length(); i++) {
            char current = temp.charAt(i);
            if ((current == '_') || (current == '-') || (current == ' ')) {
                builder.append(' ');
                newWord = true;
            } else if (Character.isUpperCase(builder.charAt(builder.length() - 1))) {
                builder.append(current);
                newWord = false;
            } else if (Character.isUpperCase(current)) {
                builder.append(' ');
                builder.append(current);
                newWord = false;
            } else if (newWord) {
                builder.append(Character.toUpperCase(current));
                newWord = false;
            } else {
                builder.append(current);
            }
        }
        return builder.toString();
    }

    @Override
    public DataSet processData(DataSet input) {
        try {
            new TasselPipeline(myArgs, myFrame, true, myButtonName);
            return null;
        } catch (Exception e) {
            myLogger.debug(e.getMessage(), e);
            throw new IllegalStateException("WorkflowPlugin: Problem running workflow: " + myButtonName + "\n" + e.getMessage());
        } finally {
            fireProgress(100);
        }
    }

    @Override
    public ImageIcon getIcon() {
        URL imageURL = WorkflowPlugin.class.getResource("/net/maizegenetics/analysis/images/Workflow.gif");
        if (imageURL == null) {
            return null;
        } else {
            return new ImageIcon(imageURL);
        }
    }

    @Override
    public String getButtonName() {
        return myButtonName;
    }

    @Override
    public String getToolTipText() {
        return myButtonName;
    }

    @Override
    public String getCitation() {
        return "Casstevens T, Wang Y. (2015) First Annual Tassel Hackathon.";
    }

}

/*
 *  TasselLogging
 */
package net.maizegenetics.tassel;

import net.maizegenetics.gui.DialogUtils;
import net.maizegenetics.gui.FileBrowserUtils;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.prefs.TasselPrefs;
import net.maizegenetics.util.LoggingUtils;
import net.maizegenetics.util.Utils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.io.OutputStream;
import java.io.PrintStream;
import java.io.UnsupportedEncodingException;
import java.net.URL;
import java.util.*;
import java.util.List;

/**
 *
 * @author Terry Casstevens
 */
public class TasselLogging extends AbstractPlugin {

    private static TasselLogging myInstance = null;
    private static final Logger myLogger = LogManager.getLogger(TasselLogging.class);

    private final JDialog myDialog = new JDialog((Window) null, "Tassel Logging", Dialog.ModalityType.MODELESS);
    private final JTextArea myTextArea = new JTextArea();
    private final TextAreaOutputStream myTextAreaOutputStream = new TextAreaOutputStream(myTextArea);
    private final PrintStream myPrintStream = new PrintStream(myTextAreaOutputStream);

    private TasselLogging(Frame parentFrame) {
        super(parentFrame, true);
        createDialog();
        basicLoggingInfo();
        LoggingUtils.setupLogging(myPrintStream);
        myTextAreaOutputStream.clear();
        if (TasselPrefs.getLogSendToConsole()) {
            LoggingUtils.setupLogging();
        }
    }

    public static TasselLogging getInstance(Frame parentFrame) {
        if (parentFrame == null) {
            return null;
        }
        if (myInstance == null) {
            myInstance = new TasselLogging(parentFrame);
        }
        return myInstance;
    }

    public static void closeInstance() {
        if (myInstance != null) {
            myInstance.close();
        }
    }

    private void close() {
        TasselPrefs.putLogXDim(myDialog.getWidth());
        TasselPrefs.putLogYDim(myDialog.getHeight());
        myDialog.setVisible(false);
        if (TasselPrefs.getLogSendToConsole()) {
            LoggingUtils.setupLogging();
        }
    }

    public static void updateLoggingLocation() {
        if (myInstance != null) {
            myInstance.updateLogging();
        }
    }

    private void updateLogging() {
        if (TasselPrefs.getLogSendToConsole()) {
            LoggingUtils.setupLogging();
        } else {
            LoggingUtils.setupDebugLogging(myPrintStream);
        }
    }

    private void createDialog() {

        myDialog.setLayout(new BorderLayout());
        int x = TasselPrefs.getLogXDim();
        int y = TasselPrefs.getLogYDim();
        if ((x < 50) || (y < 50)) {
            myDialog.setSize(500, 400);
        } else {
            myDialog.setSize(x, y);
        }

        myTextArea.setLineWrap(true);
        myTextArea.setMargin(new Insets(10, 10, 10, 10));
        myTextArea.setEditable(false);

        final JCheckBox isDebug = new JCheckBox("Debug Level");
        isDebug.setSelected(TasselPrefs.getLogDebug());
        isDebug.setToolTipText("Set to show Debug Logging Messages");
        isDebug.addActionListener(new ActionListener() {

            @Override
            public void actionPerformed(ActionEvent e) {
                boolean debugMode = isDebug.isSelected();
                isDebug.setSelected(debugMode);
                TasselPrefs.putLogDebug(debugMode);
                LoggingUtils.setupLogging(myPrintStream);
            }
        });

        JButton closeButton = new JButton();
        closeButton.setActionCommand("Close");
        closeButton.setText("Close");
        closeButton.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                close();
            }
        });

        JButton clearButton = new JButton();
        clearButton.setActionCommand("Clear");
        clearButton.setText("Clear");
        clearButton.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                myTextAreaOutputStream.clear();
            }
        });

        JButton saveButton = new JButton("Save");
        saveButton.addActionListener(new ActionListener() {

            @Override
            public void actionPerformed(ActionEvent e) {
                try {
                    File theFile = FileBrowserUtils.getSaveFile(myDialog);
                    if (theFile != null) {
                        myTextArea.write(Utils.getBufferedWriter(theFile));
                    }
                } catch (Exception ex) {
                    DialogUtils.showError(ex.getMessage() + "\n", myDialog);
                }
            }

        });

        JPanel pnlButtons = new JPanel();
        pnlButtons.setLayout(new FlowLayout());
        pnlButtons.add(closeButton);
        pnlButtons.add(clearButton);
        pnlButtons.add(saveButton);
        pnlButtons.add(isDebug);

        myDialog.getContentPane().add(new JScrollPane(myTextArea), BorderLayout.CENTER);
        myDialog.getContentPane().add(pnlButtons, BorderLayout.SOUTH);

        myDialog.setResizable(true);

    }

    @Override
    public DataSet performFunction(DataSet input) {
        try {
            LoggingUtils.setupLogging(myPrintStream);
            myDialog.setLocationRelativeTo(getParentFrame());
            myDialog.setVisible(true);
            return null;
        } finally {
            fireProgress(100);
        }
    }

    @Override
    public ImageIcon getIcon() {
        URL imageURL = TasselLogging.class.getResource("/net/maizegenetics/analysis/images/log.gif");
        if (imageURL == null) {
            return null;
        } else {
            return new ImageIcon(imageURL);
        }
    }

    @Override
    public String getButtonName() {
        return "Logging";
    }

    @Override
    public String getToolTipText() {
        return "Logging";
    }

    public static void basicLoggingInfo() {
        myLogger.info("Tassel Version: " + TASSELMainFrame.version + "  Date: " + TASSELMainFrame.versionDate);
        myLogger.info("Max Available Memory Reported by JVM: " + Utils.getMaxHeapSizeMB() + " MB");
        myLogger.info("Java Version: " + System.getProperty("java.version"));
        myLogger.info("OS: " + System.getProperty("os.name"));
        myLogger.info("Number of Processors: " + Runtime.getRuntime().availableProcessors());
        myLogger.info("Tassel Citation: " + AbstractPlugin.DEFAULT_CITATION);
        myLogger.info("");
        Set<Map.Entry<String, TasselVersions.LibraryInfo>> infos = TasselVersions.INSTANCE.libraryInfos();
        for (Map.Entry<String, TasselVersions.LibraryInfo> info : infos) {
            myLogger.info("Tassel Using Library: " + info.getValue().getName() + ": Version: " + info.getValue().getVersion() + " Date: " + info.getValue().getDate());
            myLogger.info(info.getKey() + " Citation: " + info.getValue().getCitation());
        }
    }

    class TextAreaOutputStream extends OutputStream {

        private final byte[] myByteArray = new byte[1];
        private TextAppender myTextAppender;

        public TextAreaOutputStream(JTextArea textArea) {
            myTextAppender = new TextAppender(textArea);
        }

        public synchronized void clear() {
            if (myTextAppender != null) {
                myTextAppender.clear();
                basicLoggingInfo();
            }
        }

        @Override
        public synchronized void close() {
            myTextAppender = null;
        }

        @Override
        public synchronized void flush() {
        }

        @Override
        public synchronized void write(int val) {
            myByteArray[0] = (byte) val;
            write(myByteArray, 0, 1);
        }

        @Override
        public synchronized void write(byte[] ba) {
            write(ba, 0, ba.length);
        }

        @Override
        public synchronized void write(byte[] ba, int str, int len) {
            if (myTextAppender != null) {
                myTextAppender.append(bytesToString(ba, str, len));
            }
        }

    }

    static private String bytesToString(byte[] ba, int str, int len) {
        try {
            return new String(ba, str, len, "UTF-8");
        } catch (UnsupportedEncodingException e) {
            return new String(ba, str, len);
        }
    }

    static class TextAppender implements Runnable {

        private static final int MAX_NUM_LINES = 1000;

        private final JTextArea myTextArea;
        private final LinkedList<Integer> myLineLengths = new LinkedList<>();
        private final List<String> myBufferedText = new ArrayList<>();

        private int myCurrentLineLength = 0;
        private boolean myClearAllText = false;
        private boolean myIfNoTextQueued = true;

        TextAppender(JTextArea textArea) {
            myTextArea = textArea;
        }

        synchronized void append(String val) {
            myBufferedText.add(val);
            if (myIfNoTextQueued) {
                myIfNoTextQueued = false;
                EventQueue.invokeLater(this);
            }
        }

        synchronized void clear() {
            myClearAllText = true;
            myCurrentLineLength = 0;
            myLineLengths.clear();
            myBufferedText.clear();
            if (myIfNoTextQueued) {
                myIfNoTextQueued = false;
                EventQueue.invokeLater(this);
            }
        }

        @Override
        public synchronized void run() {
            if (myClearAllText) {
                myTextArea.setText(null);
            }
            for (String val : myBufferedText) {
                myCurrentLineLength += val.length();
                if (val.endsWith(END_OF_LINE) || val.endsWith(SYSTEM_END_OF_LINE)) {
                    if (myLineLengths.size() >= MAX_NUM_LINES) {
                        myTextArea.replaceRange("", 0, myLineLengths.removeFirst());
                    }
                    myLineLengths.addLast(myCurrentLineLength);
                    myCurrentLineLength = 0;
                }
                myTextArea.append(val);
            }
            myBufferedText.clear();
            myClearAllText = false;
            myIfNoTextQueued = true;
        }

        static private final String END_OF_LINE = "\n";
        static private final String SYSTEM_END_OF_LINE = System.getProperty("line.separator", END_OF_LINE);
    }

}

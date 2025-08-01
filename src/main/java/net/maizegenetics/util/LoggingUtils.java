/*
 *  LoggingUtils
 */
package net.maizegenetics.util;

import net.maizegenetics.prefs.TasselPrefs;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;

/**
 * @author Terry Casstevens
 */
public class LoggingUtils {

    private static final Logger myLogger = LogManager.getLogger(LoggingUtils.class);
    private static final PrintStream myOriginalOutputStream = System.out;
    private static final PrintStream myOriginalErrStream = System.err;

    private static PrintStream myPrintStream;

    private LoggingUtils() {
        // Utility Class
    }

    /**
     * Sends all logging messages (including log4j) to standard out. The logging
     * level will be INFO, unless executed by the GUI and the user preference
     * has been set level to DEBUG.
     */
    public static void setupLogging() {
        if (TasselPrefs.getLogDebug()) {
            setupDebugLogging();
        } else {
            System.setOut(myOriginalOutputStream);
            System.setErr(myOriginalErrStream);
            LoggingUtilsKt.sendLog4jToStdout();
        }
    }

    /**
     * Sends all logging messages (including log4j) to standard out. The logging
     * level will be DEBUG.
     */
    public static void setupDebugLogging() {
        System.setOut(myOriginalOutputStream);
        System.setErr(myOriginalErrStream);
        LoggingUtilsKt.sendDebugLog4jToStdout();
    }

    /**
     * Sends all logging messages (including log4j) to specified PrintStream.
     * The logging level will be INFO, unless executed by the GUI and the user
     * preference has been set level to DEBUG. This is used to direct logging
     * messages to the GUI.
     *
     * @param stream stream
     */
    public static void setupLogging(PrintStream stream) {
        if (TasselPrefs.getLogDebug()) {
            setupDebugLogging(stream);
        } else {
            System.setOut(stream);
            System.setErr(stream);
            LoggingUtilsKt.sendLog4jToStdout();
        }
    }

    /**
     * Sends all logging messages (including log4j) to specified PrintStream.
     * The logging level will be DEBUG. This is used to direct logging messages
     * to the GUI.
     *
     * @param stream
     */
    public static void setupDebugLogging(PrintStream stream) {
        System.setOut(stream);
        System.setErr(stream);
        LoggingUtilsKt.sendDebugLog4jToStdout();
    }

    /**
     * Sends all logging messages (including log4j) to specified file. The
     * logging level will be INFO, unless executed by the GUI and the user
     * preference has been set level to DEBUG.
     *
     * @param logFileName file name
     *
     * @throws FileNotFoundException
     */
    public static void setupLogfile(String logFileName) throws FileNotFoundException {
        if (TasselPrefs.getLogDebug()) {
            setupDebugLogfile(logFileName);
        } else {
            File logFile = new File(logFileName);
            myPrintStream = new PrintStream(logFile);
            System.setOut(myPrintStream);
            System.setErr(myPrintStream);
            LoggingUtilsKt.sendLog4jToStdout();
            myLogger.info("Log File: " + logFile.getAbsolutePath());
        }
    }

    /**
     * Sends all logging messages (including log4j) to specified file. The
     * logging level will be DEBUG.
     *
     * @param logFileName file name
     *
     * @throws FileNotFoundException
     */
    public static void setupDebugLogfile(String logFileName) throws FileNotFoundException {
        File logFile = new File(logFileName);
        myPrintStream = new PrintStream(logFile);
        System.setOut(myPrintStream);
        System.setErr(myPrintStream);
        LoggingUtilsKt.sendDebugLog4jToStdout();
        myLogger.info("Log File: " + logFile.getAbsolutePath());
    }

    public static void closeLogfile() {
        if (myPrintStream != null) {
            myPrintStream.close();
        }
    }

}

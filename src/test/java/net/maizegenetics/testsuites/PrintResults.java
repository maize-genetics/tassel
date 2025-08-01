/*
 * PrintResults
 */
package net.maizegenetics.testsuites;

import java.util.List;
import org.junit.runner.Result;
import org.junit.runner.notification.Failure;

/**
 *
 * @author terry
 */
public class PrintResults {

    private PrintResults() {
        // utility
    }

    public static void printResults(Result result) {
        System.out.println("\n\n");
        System.out.println("Run Count: " + result.getRunCount());
        System.out.println("Failure Count: " + result.getFailureCount());
        List<Failure> failures = result.getFailures();
        System.out.println("");
        for (int i = 0; i < failures.size(); i++) {
            Failure current = failures.get(i);
            System.out.println("Failure Description: " + current.getDescription());
            System.out.println("Failure Test Header: " + current.getTestHeader());
            System.out.println("Failure Message: " + current.getMessage());
            System.out.println("Failure Trace: " + current.getTrace());
            System.out.println("");
        }
    }
}

/*
 *  AMatrixPluginTest
 * 
 *  Created on Oct 21, 2015
 */
package net.maizegenetics.analysis.distance;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.List;
import org.junit.Test;
import static org.junit.Assert.*;
import org.junit.rules.TemporaryFolder;

/**
 *
 * @author Josh Lamos-Sweeney
 */
public class AMatrixPluginTest {

    public AMatrixPluginTest() {
    }

    @Test
    public void testInbredMatrix() {
        String[][] testingPedigree = {
            {"A", "B", "C"},
            {"C", "0", "0"},
            {"D", "B", "C"},
            {"E", "A", "D"},//Half Sibling
            {"F", "A", "D"},//Half Sibling
            {"G", "E", "F"}};

        //expected output matrix: A=0 D=1 E=2 F=3 G=4
        AMatrixPlugin test = new AMatrixPlugin(null, false);
        double[][] out = test.pedMatrix(testingPedigree);
        List<String> progenyIDs = test.myProgenyIDs;
        progenyIDs.indexOf("F");
        int A = progenyIDs.indexOf("A");
        int B = progenyIDs.indexOf("B");
        int C = progenyIDs.indexOf("C");
        int D = progenyIDs.indexOf("D");
        int E = progenyIDs.indexOf("E");
        int F = progenyIDs.indexOf("F");
        int G = progenyIDs.indexOf("G");
        assertTrue("No relationship between B and C", out[B][C] == 0);
        assertTrue("A and D are double siblings", out[A][D] == 0.5);
        assertTrue("A is the direct child of B", out[A][B] == 0.5);
        assertEquals("F's inbreeding coefficient is one and a quarter", 1.25, out[F][F], 0.0001);
        assertEquals("The relationship between F and G has to take into account F's high inbreeding coefficient", 1.0, out[F][G], 0.0001);
        assertTrue("A to F has 6 relationships, 2 .25's and 4 .0625's", out[A][F] == 0.75);
    }

    @Test
    public void testReader() throws Exception {
        String testString = "family A B C\nfamily D B  C junk data is here\n"
                + "family    E A	D junk\nfaimly	F	A	D\nFam G 	E	F";

        TemporaryFolder tf = new TemporaryFolder();
        tf.create();
        File test = tf.newFile("IOTest");
        BufferedWriter bw = new BufferedWriter(new FileWriter(test));
        bw.write(testString);
        bw.close();
        String[][] out = AMatrixPlugin.plinkToPed(test.getAbsolutePath());
        //System.out.println(Arrays.deepToString(out));
        assertTrue(out[0][0].equals("A"));
        assertTrue(out[2][0].equals("E"));
        assertTrue(out[2][2].equals("D"));
        assertTrue(out[4][0].equals("G"));
        assertTrue(out[4][1].equals("E"));
    }

}

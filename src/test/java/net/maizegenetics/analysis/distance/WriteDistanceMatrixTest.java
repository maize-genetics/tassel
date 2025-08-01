package net.maizegenetics.analysis.distance;

import static org.junit.Assert.*;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.file.Files;
import java.nio.file.Paths;

import net.maizegenetics.constants.GeneralConstants;
import net.maizegenetics.constants.TutorialConstants;
import net.maizegenetics.dna.snp.FilterGenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableUtils;
import net.maizegenetics.dna.snp.ImportUtils;
import net.maizegenetics.taxa.distance.DistanceMatrix;
import net.maizegenetics.taxa.distance.WriteDistanceMatrix;

import org.junit.Test;

public class WriteDistanceMatrixTest {

    @Test
    public void rawFileExportTest() {
        //Load up Tutorial Genotype Data
        try {
            System.out.println("Testing Raw Export");
            System.out.println("Reading in Genotype Table");
            File input = new File(TutorialConstants.HAPMAP_FILENAME);
            if (!input.exists()) {
                fail("Input File: " + TutorialConstants.HAPMAP_FILENAME + " doesn't exist.");
            }
            GenotypeTable inputAlign = ImportUtils.readFromHapmap(TutorialConstants.HAPMAP_FILENAME, null);

            GenotypeTable filteredGeno = GenotypeTableUtils.removeSitesBasedOnFreqIgnoreMissing(inputAlign, 0.01, 1, 250);
            //remove the marker with 4 alleles
            filteredGeno = FilterGenotypeTable.getInstanceRemoveSiteNames(filteredGeno, new String[]{"PZB00063.1"});

            System.out.println("Creating Kinship Matrix");
            //Create Kinship
            DistanceMatrix kin = EndelmanDistanceMatrix.getInstance(filteredGeno);
            //DistanceMatrix kin = Kinship.createKinship(filteredGeno, KINSHIP_TYPE.Endelman, GENOTYPE_TABLE_COMPONENT.Genotype);

            //Export file to tmp directory
            System.out.println("Exporting Kinship in LDAK Raw Format");
            String taxaFile = GeneralConstants.TEMP_DIR + "mdp_MultiBLUP_raw.grm.id";
            String matrixFile = GeneralConstants.TEMP_DIR + "mdp_MultiBLUP_raw.grm.raw";
            WriteDistanceMatrix.saveRawMultiBlupMatrix(kin, taxaFile, matrixFile);

            //Load both exported files and Comparison files
            //ID files
            BufferedReader exportedID = new BufferedReader(new FileReader(GeneralConstants.TEMP_DIR + "mdp_MultiBLUP_raw.grm.id"));
            BufferedReader expectedID = new BufferedReader(new FileReader(GeneralConstants.EXPECTED_RESULTS_DIR + "mdp_MultiBLUP_raw.grm.id"));

            //Raw files
            BufferedReader exportedRaw = new BufferedReader(new FileReader(GeneralConstants.TEMP_DIR + "mdp_MultiBLUP_raw.grm.raw"));
            BufferedReader expectedRaw = new BufferedReader(new FileReader(GeneralConstants.EXPECTED_RESULTS_DIR + "mdp_MultiBLUP_raw.grm.raw"));

            //Compare
            String currExportedLine = "";
            String currExpectedLine = "";
            System.out.print("Comparing IDs");
            boolean equalCheck = true;
            int counter = 0;
            String errorHolder = "";
            while ((currExportedLine = exportedID.readLine()) != null && (currExpectedLine = expectedID.readLine()) != null) {
                if (!currExportedLine.equals(currExpectedLine)) {
                    equalCheck = false;
                    errorHolder += counter + ", ";
                }
                counter++;
            }
            assertTrue("IDs do not match at line(s):\n" + errorHolder, equalCheck);
            System.out.println(",Done");
            System.out.print("Comparing Kinship Values");
            equalCheck = true;
            counter = 0;
            errorHolder = "";
            currExportedLine = "";
            currExpectedLine = "";
            while ((currExportedLine = exportedRaw.readLine()) != null && (currExpectedLine = expectedRaw.readLine()) != null) {
                if (!currExportedLine.equals(currExpectedLine)) {
                    equalCheck = false;
                    errorHolder += counter + ", ";
                }
                counter++;
            }
            assertTrue("Kinship Values do not match at line(s):\n" + errorHolder, equalCheck);
            System.out.println(",Done");

            exportedID.close();
            exportedRaw.close();
            expectedID.close();
            expectedRaw.close();
        } catch (Exception e) {
            System.out.println(e);
        }
    }

    @Test
    public void binFileExportTest() {
        //Load up Tutorial Genotype Data
        try {
            System.out.println("Testing Binary Export");
            System.out.println("Reading in Genotype Table");
            File input = new File(TutorialConstants.HAPMAP_FILENAME);
            if (!input.exists()) {
                fail("Input File: " + TutorialConstants.HAPMAP_FILENAME + " doesn't exist.");
            }
            GenotypeTable inputAlign = ImportUtils.readFromHapmap(TutorialConstants.HAPMAP_FILENAME, null);

            GenotypeTable filteredGeno = GenotypeTableUtils.removeSitesBasedOnFreqIgnoreMissing(inputAlign, 0.01, 1, 250);
            //remove the marker with 4 alleles
            filteredGeno = FilterGenotypeTable.getInstanceRemoveSiteNames(filteredGeno, new String[]{"PZB00063.1"});

            System.out.println("Creating Kinship Matrix");
            //Create Kinship
            DistanceMatrix kin = EndelmanDistanceMatrix.getInstance(filteredGeno);
            //DistanceMatrix kin = Kinship.createKinship(filteredGeno, KINSHIP_TYPE.Endelman, GENOTYPE_TABLE_COMPONENT.Genotype);

            //Export file to tmp directory
            System.out.println("Exporting Kinship in LDAK Raw Format");
            String taxaFile = GeneralConstants.TEMP_DIR + "mdp_MultiBLUP_bin.grm.id";
            String matrixFile = GeneralConstants.TEMP_DIR + "mdp_MultiBLUP_bin.grm.bin";
            String countFile = GeneralConstants.TEMP_DIR + "mdp_MultiBLUP_bin.grm.N.bin";
            WriteDistanceMatrix.saveBinMultiBlupMatrix(kin, taxaFile, matrixFile, countFile);

            //Load both exported files and Comparison files
            //ID files
            BufferedReader exportedID = new BufferedReader(new FileReader(GeneralConstants.TEMP_DIR + "mdp_MultiBLUP_bin.grm.id"));
            BufferedReader expectedID = new BufferedReader(new FileReader(GeneralConstants.EXPECTED_RESULTS_DIR + "mdp_MultiBLUP_bin.grm.id"));

            //Compare
            String currExportedLine = "";
            String currExpectedLine = "";
            System.out.print("Comparing IDs");
            boolean equalCheck = true;
            int counter = 0;
            String errorHolder = "";
            while ((currExportedLine = exportedID.readLine()) != null && (currExpectedLine = expectedID.readLine()) != null) {
                if (!currExportedLine.equals(currExpectedLine)) {
                    equalCheck = false;
                    errorHolder += counter + ", ";
                }
                counter++;
            }
            assertTrue("IDs do not match at line(s):\n" + errorHolder, equalCheck);
            System.out.println(",Done");

            System.out.print("Comparing Kinship Values");
            byte[] bytesExported = Files.readAllBytes(Paths.get(GeneralConstants.TEMP_DIR + "mdp_MultiBLUP_bin.grm.bin"));
            byte[] bytesExpected = Files.readAllBytes(Paths.get(GeneralConstants.EXPECTED_RESULTS_DIR + "mdp_MultiBLUP_bin.grm.bin"));

            equalCheck = true;
            errorHolder = "";
            assertEquals("Size of Kinship Files is not the same", bytesExported.length, bytesExpected.length);
            if (bytesExported.length == bytesExpected.length) {

                for (int i = 0; i < bytesExported.length; i += 4) {
                    //Cannot simply check byte to byte as they are floating points
                    //This allows for slight deviation between bytes
                    byte[] exportedByte = {bytesExported[i], bytesExported[i + 1], bytesExported[i + 2], bytesExported[i + 3]};
                    byte[] expectedByte = {bytesExpected[i], bytesExpected[i + 1], bytesExpected[i + 2], bytesExpected[i + 3]};

                    float exportedNumber = ByteBuffer.wrap(exportedByte).order(ByteOrder.LITTLE_ENDIAN).getFloat();
                    float expectedNumber = ByteBuffer.wrap(expectedByte).order(ByteOrder.LITTLE_ENDIAN).getFloat();

                    if (Math.abs(exportedNumber - expectedNumber) > .0000001) {
                        equalCheck = false;
                        errorHolder += i + ", ";
                        System.out.println(exportedNumber + " " + expectedNumber);
                    }
                }
            }
            assertTrue("Kinship Values do not match at byte(s):\n" + errorHolder, equalCheck);
            System.out.println(",Done");

            System.out.print("Comparing Counts");
            bytesExported = Files.readAllBytes(Paths.get(GeneralConstants.TEMP_DIR + "mdp_MultiBLUP_bin.grm.N.bin"));
            bytesExpected = Files.readAllBytes(Paths.get(GeneralConstants.EXPECTED_RESULTS_DIR + "mdp_MultiBLUP_bin.grm.N.bin"));

            equalCheck = true;
            errorHolder = "";
            assertEquals("Size of Count Files is not the same", bytesExported.length, bytesExpected.length);
            if (bytesExported.length == bytesExpected.length) {
                for (int i = 0; i < bytesExported.length; i++) {
                    if (bytesExported[i] != bytesExpected[i]) {
                        equalCheck = false;
                        errorHolder += i + ", ";
                    }
                }
            }
            assertTrue("Count Values do not match at byte(s):\n" + errorHolder, equalCheck);
            System.out.println(",Done");

            exportedID.close();
            expectedID.close();

        } catch (Exception e) {
            System.out.println(e);
        }
    }
}

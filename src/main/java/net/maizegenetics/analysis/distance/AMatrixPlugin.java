/*
 *  AMatrixPlugin
 * 
 *  Created on Oct 20, 2015
 */
package net.maizegenetics.analysis.distance;

import java.awt.Frame;
import java.io.BufferedReader;
import java.net.URL;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import javax.swing.ImageIcon;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.PluginParameter;
import net.maizegenetics.taxa.TaxaListBuilder;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.taxa.distance.DistanceMatrix;
import net.maizegenetics.util.Utils;

/**
 *
 * @author Josh Lamos-Sweeney
 * @author Yaw Nti-Addae
 * @author Kelly Robbins
 * @author Terry Casstevens
 */
public class AMatrixPlugin extends AbstractPlugin {

    private PluginParameter<String> myPedFilename = new PluginParameter.Builder<>("pedigreeFilename", null, String.class)
            .description("Pedigree Filename")
            .required(true)
            .inFile()
            .build();

    private double[][] myAMatrix;
    private HashMap<Integer, Progeny> myProgeny;
    protected List<String> myProgenyIDs;

    public AMatrixPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    @Override
    public String pluginDescription() {
        return "Create Pedigree Matrix";
    }

    @Override
    public DataSet processData(DataSet input) {

        double[][] result = pedMatrix(plinkToPed(pedFilename()));
        TaxaListBuilder builder = new TaxaListBuilder();
        for (String current : myProgenyIDs) {
            builder.add(new Taxon(current));
        }
        DistanceMatrix matrix = new DistanceMatrix(result, builder.build());

        return new DataSet(new Datum("Pedigree Matrix for " + Utils.getFilename(pedFilename()), matrix, null), this);

    }

    /**
     * Takes a filepath to a PLINK formatted file(see below). Converts to a
     * String[][] format where each row contains 3 columns: Individual
     * Identifier, Paternal Identifier, Maternal Identifier This format is
     * convenient in calculating an pedigree-based kinship (A) matrix
     *
     * Expected PLINK pedigree file format (.ped) Family ID Individual ID
     * Paternal ID Maternal ID and several more fields These fields are
     * separated by one or more whitespace characters, and each row is on its
     * own line. For our purposes, we expect every Individual ID in the file to
     * be unique
     *
     * @param ped filepath to a valid PLINK .ped file, or compatible file format
     * @return A n by 3 matrix representing a pedigree as progeny,parent,parent
     */
    public static String[][] plinkToPed(String ped) {
        try {
            List<String[]> rows = new ArrayList<>();
            BufferedReader br = Utils.getBufferedReader(ped);
            while (br.ready()) {
                String line = br.readLine();
                //Fields are Family ID, Individual ID, Paternal ID, Maternal ID, junk, junk <AFAIC>
                String[] fields = line.split("\\s+");
                if (fields.length > 3) {
                    String progenyID = fields[1];
                    String parent1ID = fields[2];
                    String parent2ID = fields[3];
                    String[] resultRow = {progenyID, parent1ID, parent2ID};
                    rows.add(resultRow);
                }
            }
            String[][] result = new String[rows.size()][3];
            for (int i = 0; i < rows.size(); i++) {
                result[i] = rows.get(i);
            }
            return result;
        } catch (Exception e) {
            throw new IllegalStateException("plinkToPed: problem reading file: " + ped);
        }
    }

    /**
     * Calculates an A matrix from a pedigree (see {@link #plinkToPed(String)})
     *
     * @param pedigree A n by 3 matrix, where each row is of the format
     * myProgeny, parent1, parent2, where each is its unique string identifier.
     * Identifiers of * blank, "0" are treated as unknown parents.
     * @return A pedigree matrix, sorted alphabetically, with one row for each
     * non-unknown parent.
     */
    public double[][] pedMatrix(String[][] pedigree) {
        //Read list of myProgeny
        myProgenyIDs = getNameList(pedigree);

        myProgeny = new HashMap<>();
        for (String[] p : pedigree) {
            int progenyID = myProgenyIDs.indexOf(p[0]);//Will always be found
            myProgeny.put(progenyID, new Progeny(progenyID, myProgenyIDs.indexOf(p[1]), myProgenyIDs.indexOf(p[2])));
        }
        for (int i = 0; i < myProgenyIDs.size(); i++) {
            if (!myProgeny.containsKey(i)) {
                myProgeny.put(i, new Progeny(i, -1, -1));
            }
        }
        int size = myProgenyIDs.size();
        myAMatrix = new double[size][size];
        for (int i = 0; i < size; i++) {
            Arrays.fill(myAMatrix[i], Double.NaN);
        }
        for (int i = 0; i < size; i++) {
            for (int j = i; j < size; j++) {
                calcMatrix(i, j);
            }
        }
        return myAMatrix;
    }

    /**
     * Given a pedigree in {@link #plinkToPed(String)} format, returns the
     * ordered list of identifiers used as indices in the matrix returned by
     * {@link #pedMatrix(String[][])}.
     *
     * @param pedigree A n by 3 matrix, where each row is of the format
     * myProgeny, parent1, parent2, where each is its unique string identifier.
     * Identifiers of blank, "0" are treated as unknown parents.
     * @return the ordered list of identifiers used as indices in the
     * pedigree-based A matrix.
     */
    public List<String> getNameList(String[][] pedigree) {
        HashSet<String> progenySet = new HashSet<>();
        for (String[] p : pedigree) {
            progenySet.add(p[0]);
            progenySet.add(p[1]);
            progenySet.add(p[2]);
        }
        //Remove unknown parents hack
        progenySet.remove("0");
        return new ArrayList<>(progenySet);
    }

    /**
     * Given an X,Y location on the matrix, calculates the location if it's not
     * calculated (also calculates any dependencies) Recursive call for
     * {@link #pedMatrix(String[][]) pedMatrix}
     *
     * @param x First relationship
     * @param y Second relationship
     * @return value of the entry (useful for recursion)
     */
    private double calcMatrix(int x, int y) {
        if (x == -1 || y == -1) {
            return 0;
        }
        double result;
        if (!Double.isNaN(myAMatrix[x][y])) {
            return myAMatrix[x][y];
        }
        Progeny X = myProgeny.get(x);
        Progeny Y = myProgeny.get(y);
        if (x == y) {
            result = 1 + (calcMatrix(X.parent1, X.parent2) / 2);
        } else {
            int F1 = X.parent1;
            int F2 = Y.parent1;
            int M1 = X.parent2;
            int M2 = Y.parent2;
            if (F1 == -1 || M1 == -1) {
                result = (calcMatrix(x, M2) + calcMatrix(F2, x)) / 2;
            } else if (F2 == -1 || M2 == -1) {
                result = (calcMatrix(F1, y) + calcMatrix(M1, y)) / 2;
            } else {
                double result1, result2;
                result1 = (calcMatrix(x, M2) + calcMatrix(F2, x)) / 2;
                result2 = (calcMatrix(F1, y) + calcMatrix(M1, y)) / 2;
                result = result1 > result2 ? result1 : result2;
            }
        }
        myAMatrix[x][y] = result;
        myAMatrix[y][x] = result;
        return result;
    }

    /**
     * Simple holder class for myProgeny-parent relationships.
     *
     * @author Josh Lamos-Sweeney
     */
    private class Progeny {

        public int progeny;
        public int parent1;
        public int parent2;

        public Progeny(int progeny, int parent1, int parent2) {
            this.progeny = progeny;
            this.parent1 = parent1;
            this.parent2 = parent2;
        }
    }

    /**
     * Create A Matrix
     *
     * @return Ped Filename
     */
    public String pedFilename() {
        return myPedFilename.value();
    }

    /**
     * Set Ped Filename. Create A Matrix
     *
     * @param value Ped Filename
     *
     * @return this plugin
     */
    public AMatrixPlugin pedFilename(String value) {
        myPedFilename = new PluginParameter<>(myPedFilename, value);
        return this;
    }

    @Override
    public ImageIcon getIcon() {
        URL imageURL = AMatrixPlugin.class.getResource("/net/maizegenetics/analysis/images/amatrix.png");
        if (imageURL == null) {
            return null;
        } else {
            return new ImageIcon(imageURL);
        }
    }

    @Override
    public String getButtonName() {
        return "Pedigree Relationship Matrix";
    }

    @Override
    public String getToolTipText() {
        return "Create Pedigree Relationship Matrix";
    }

    @Override
    public String getCitation() {
        return "Lamos-Sweeney J, Nti-Addae Y, Robbins K, Casstevens T. (Oct. 2015) Second Tassel Hackathon.";
    }

}

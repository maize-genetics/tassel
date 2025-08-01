/*
 * AbstractAlignmentMask
 */
package net.maizegenetics.gui;

import java.awt.Color;

import net.maizegenetics.dna.snp.GenotypeTable;

/**
 *
 * @author Terry Casstevens
 */
abstract public class AbstractGenotypeTableMask implements GenotypeTableMask {

    private static final long serialVersionUID = -5197800047652332969L;
    private final String myName;
    private Color myColor;
    private final GenotypeTable myAlignment;
    private static Color myLastColor = null;
    private static int myIncrementAmount = 0x7D;
    private final MaskType myMaskType;

    public AbstractGenotypeTableMask(GenotypeTable align, String name, Color color, MaskType type) {
        myAlignment = align;
        myName = name;
        myColor = color;
        myMaskType = type;
    }

    protected static Color getNextColor() {

        if (myLastColor == null) {
            myLastColor = new Color(0x7D, 0, 0);
            return myLastColor;
        }

        int newColor = myLastColor.getRGB() + myIncrementAmount;
        newColor = newColor % 0xFFFFFF;

        myIncrementAmount = (myIncrementAmount << 8) % 0xFFFFFF;
        myIncrementAmount = myIncrementAmount == 0 ? 0x7D : myIncrementAmount;

        myLastColor = new Color(newColor);
        return myLastColor;

    }

    public Color getColor() {
        return myColor;
    }

    public void setColor(Color color) {
        myColor = color;
    }

    public GenotypeTable getAlignment() {
        return myAlignment;
    }

    public MaskType getMaskType() {
        return myMaskType;
    }

    public String toString() {
        if (myColor == null) {
            return myName;
        } else {
            String color = " (Red: " + myColor.getRed() + "  Green: " + myColor.getGreen() + "  Blue: " + myColor.getBlue() + ")";
            return myName + color;
        }
    }
}

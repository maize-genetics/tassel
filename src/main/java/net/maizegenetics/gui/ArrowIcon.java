/*
 * ArrowIcon
 */
package net.maizegenetics.gui;

import java.awt.Color;
import java.awt.Component;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.geom.Path2D;

import javax.swing.Icon;
import javax.swing.SwingConstants;
import javax.swing.UIManager;

/**
 * Directional triangle for icon-only navigation buttons.
 *
 * <p>The triangle is painted as a path rather than loaded from a bitmap so that it stays sharp on
 * HiDPI displays, and it takes its color from the component being painted so that it follows the
 * active look-and-feel theme (light or dark) without needing per-theme artwork.
 *
 * @author Brandon Monier
 */
public class ArrowIcon implements Icon {

    private final int myDirection;
    private final int mySize;

    /**
     * @param direction {@link SwingConstants#WEST} or {@link SwingConstants#EAST}
     * @param size width and height of the icon in pixels
     */
    public ArrowIcon(int direction, int size) {
        if ((direction != SwingConstants.WEST) && (direction != SwingConstants.EAST)) {
            throw new IllegalArgumentException("ArrowIcon: init: unsupported direction: " + direction);
        }
        if (size < 1) {
            throw new IllegalArgumentException("ArrowIcon: init: size must be positive: " + size);
        }
        myDirection = direction;
        mySize = size;
    }

    @Override
    public void paintIcon(Component c, Graphics g, int x, int y) {

        Graphics2D g2d = (Graphics2D) g.create();
        try {
            g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
            g2d.setColor(getColor(c));
            g2d.fill(createTriangle(x, y));
        } finally {
            g2d.dispose();
        }

    }

    /**
     * Insetting the triangle keeps the point and the flat edge clear of the button border, and
     * leaves it visually balanced against text of the same nominal size.
     */
    private Path2D createTriangle(int x, int y) {

        double inset = mySize / 6.0;
        double min = inset;
        double max = mySize - inset;
        double mid = mySize / 2.0;

        Path2D result = new Path2D.Double();
        if (myDirection == SwingConstants.WEST) {
            result.moveTo(x + max, y + min);
            result.lineTo(x + max, y + max);
            result.lineTo(x + min, y + mid);
        } else {
            result.moveTo(x + min, y + min);
            result.lineTo(x + min, y + max);
            result.lineTo(x + max, y + mid);
        }
        result.closePath();
        return result;

    }

    private Color getColor(Component c) {
        if ((c != null) && !c.isEnabled()) {
            Color disabled = UIManager.getColor("Button.disabledText");
            return (disabled != null) ? disabled : Color.GRAY;
        }
        return (c != null) ? c.getForeground() : Color.BLACK;
    }

    @Override
    public int getIconWidth() {
        return mySize;
    }

    @Override
    public int getIconHeight() {
        return mySize;
    }
}

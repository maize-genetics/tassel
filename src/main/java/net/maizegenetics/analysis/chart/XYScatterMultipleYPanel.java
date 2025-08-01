/*
 * XYScatterMultipleYPanel
 */
package net.maizegenetics.analysis.chart;

import java.awt.Color;
import java.awt.Component;
import java.awt.Font;
import java.awt.event.ActionEvent;
import java.awt.font.TextAttribute;
import java.text.NumberFormat;
import java.util.Map;
import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JComponent;
import net.maizegenetics.util.TableReport;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.title.LegendTitle;

/**
 *
 * @author yz79
 * @author Terry Casstevens
 */
public class XYScatterMultipleYPanel extends BasicChartPanel {

    private final ManhattanDisplayPlugin myManhattanDisplayPlugin;
    private final ChartPanel myChartPanel;
    private final TableReportManhattanDataset dataset;

    public XYScatterMultipleYPanel(ManhattanDisplayPlugin plugin, TableReport theTable, int start, int end) {
        myManhattanDisplayPlugin = plugin;
        dataset = new TableReportManhattanDataset(theTable, start, end);
        chart = createChart(dataset);
        myChartPanel = new ChartPanel(chart);
        myChartPanel.setPreferredSize(new java.awt.Dimension(900, 500));
        jbInit();
    }

    private void jbInit() {
        this.setLayout(new BoxLayout(this, BoxLayout.PAGE_AXIS));
        this.add(myChartPanel);
        JButton saveButton = new JButton("save...");
        saveButton.setAlignmentX(Component.RIGHT_ALIGNMENT);
        saveButton.addActionListener((ActionEvent e) -> {
            myManhattanDisplayPlugin.saveDataToFile(myChartPanel);
        });
        this.add(saveButton);
    }

    public JFreeChart createChart(TableReportManhattanDataset dataset) {

        String name = "Please select numeric variables";
        String xName = "X";
        String yName = "Y";
        if (dataset != null) {
            xName = dataset.getXName();
            yName = "-Log10(P-Value)";
            name = "P-Values by Chromosome for " + dataset.getTrait();
        }

        int seriesCount = dataset.getSeriesCount();

        if (seriesCount < 200) {

            chart = ChartFactory.createScatterPlot(name,
                    xName, yName,
                    dataset,
                    PlotOrientation.VERTICAL,
                    true,
                    true,
                    false);

            if (seriesCount > 21) {
                int fontSize = 5;
                if (seriesCount < 30) {
                    fontSize = 11;
                } else if (seriesCount < 50) {
                    fontSize = 10;
                } else if (seriesCount < 70) {
                    fontSize = 9;
                } else if (seriesCount < 90) {
                    fontSize = 8;
                } else if (seriesCount < 110) {
                    fontSize = 7;
                } else if (seriesCount < 130) {
                    fontSize = 6;
                }
                LegendTitle title = chart.getLegend();
                Map attributes = title.getItemFont().getAttributes();
                attributes.put(TextAttribute.SIZE, fontSize);
                title.setItemFont(new Font(attributes));
            }

        } else {

            chart = ChartFactory.createScatterPlot(name,
                    xName, yName,
                    dataset,
                    PlotOrientation.VERTICAL,
                    false,
                    true,
                    false);

        }

        chart.setBackgroundPaint(new Color(238, 238, 238));

        NumberAxis domainAxis = new NumberAxis(xName);
        domainAxis.setUpperMargin(0.01);
        domainAxis.setLowerMargin(0.01);
        domainAxis.setRangeWithMargins(dataset.getDomainBounds(true));
        NumberFormat postionFormat = new ManhattanNumberFormat(NumberFormat.getIntegerInstance(), dataset.getActualPostions());
        domainAxis.setNumberFormatOverride(postionFormat);
        domainAxis.setVerticalTickLabels(true);
        chart.getXYPlot().setDomainAxis(domainAxis);
        chart.getXYPlot().setRangeGridlinePaint(Color.BLACK);

        NumberAxis rangeAxis = new NumberAxis(yName);
        rangeAxis.setLowerMargin(0.01);
        rangeAxis.setRangeWithMargins(dataset.getRangeBounds());
        NumberFormat decimal = NumberFormat.getInstance();
        decimal.setMinimumFractionDigits(2);
        rangeAxis.setNumberFormatOverride(decimal);
        rangeAxis.setMinorTickMarksVisible(true);
        chart.getXYPlot().setRangeAxis(rangeAxis);

        chart.getXYPlot().setBackgroundPaint(Color.WHITE);
        chart.getXYPlot().setForegroundAlpha(1.0f);
        chart.getXYPlot().getRenderer().setBaseToolTipGenerator(new XYMultipleYToolTipGenerator(postionFormat));
        return chart;

    }

    @Override
    public JComponent getMainComponent() {
        return myChartPanel;
    }
}

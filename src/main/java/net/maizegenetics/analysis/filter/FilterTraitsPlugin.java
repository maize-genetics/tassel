package net.maizegenetics.analysis.filter;

import java.awt.Frame;
import java.net.URL;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.stream.IntStream;

import javax.swing.ImageIcon;
import javax.swing.JOptionPane;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.phenotype.Phenotype;
import net.maizegenetics.phenotype.PhenotypeAttribute;
import net.maizegenetics.phenotype.Phenotype.ATTRIBUTE_TYPE;
import net.maizegenetics.phenotype.PhenotypeBuilder;

public class FilterTraitsPlugin extends AbstractPlugin {
    private static final Logger myLogger = LogManager.getLogger(FilterTraitsPlugin.class);
	private ArrayList<int[]> includeList;
	private ArrayList<Map<PhenotypeAttribute, ATTRIBUTE_TYPE>> typeChangeList;
	private boolean excludeLast = false;
	
	public FilterTraitsPlugin(Frame parentFrame, boolean isInteractive) {
		super(parentFrame, isInteractive);
		
	}
	
	@Override
	public String getButtonName() {
		return "Traits";
	}

	@Override
	public ImageIcon getIcon() {
            URL imageURL = FilterTraitsPlugin.class.getResource("/net/maizegenetics/analysis/images/FilterNew.gif");
        if (imageURL == null) {
            return null;
        } else {
            return new ImageIcon(imageURL);
        }
	}

	@Override
	public String getToolTipText() {
		return "Exclude traits or change properties";
	}

	@Override
	public DataSet performFunction(DataSet input) {
		includeList = new ArrayList<>();
		typeChangeList = new ArrayList<>();
		List<Datum> data = input.getDataOfType(Phenotype.class);
		ArrayList<Datum> outputList = new ArrayList<Datum>();

		if (isInteractive()) {
			if (data.size() == 0) {
				JOptionPane.showMessageDialog(getParentFrame(), "No Phenotype data selected.");
			}
			for (Datum datum : data) {
				FilterTraitsDialog ftd = new FilterTraitsDialog(getParentFrame(), (Phenotype) datum.getData());
				ftd.setLocationRelativeTo(getParentFrame());
				ftd.setVisible(true);
				if (ftd.getClickedOK()) {
					includeList.add(ftd.getIncludedTraits());
					typeChangeList.add(ftd.getTypeChangeMap());
				}
				ftd.dispose();
			}

			int n = includeList.size();
			for (int i = 0; i < n; i++) {
				Datum datum = data.get(i);
				Phenotype pheno = (Phenotype) datum.getData();
				int[] included = includeList.get(i);
				int numberOfOriginalTraits = pheno.numberOfAttributes();
				boolean buildnew = false;
				PhenotypeBuilder phenoBuilder = new PhenotypeBuilder().fromPhenotype(pheno);

				if (included.length > 0 && included.length < numberOfOriginalTraits) {
					buildnew = true;
					//setAttributesToKeep on the phenotype builder
					phenoBuilder.keepAttributes(included);
				}

				Map<PhenotypeAttribute, ATTRIBUTE_TYPE> typeChangeMap = typeChangeList.get(i);
				if (typeChangeMap.size() > 0) {
					buildnew = true;
					//add to the builder
					phenoBuilder.changeAttributeType(typeChangeMap);
				}

				if (buildnew) {

					String name = "Filtered_" + datum.getName();
					phenoBuilder.assignName(name);
					String comment = "";
					outputList.add(new Datum(name, phenoBuilder.build().get(0), comment));
				}
			}		
		} else {
			if (excludeLast) {
				for (Datum datum : data) {
					Phenotype pheno = (Phenotype) datum.getData();
					int numberOfOriginalTraits = pheno.numberOfAttributes();
					int[] keepAttributeIndex = IntStream.iterate(0, i -> i + 1).limit(numberOfOriginalTraits - 1).toArray();
					
					String name = "Filtered_" + datum.getName();
					Phenotype filteredPhenotype = new PhenotypeBuilder().fromPhenotype(pheno)
							.assignName(name)
							.keepAttributes(keepAttributeIndex)
							.build().get(0);
					outputList.add(new Datum(name, filteredPhenotype, "no comment"));
				}
			}
		}

		DataSet ds = new DataSet(outputList, this);
		fireDataSetReturned(ds);
		return ds;
	}
	
	public void addIncludedTraits(int[] traitsToInclude) {
		includeList.add(traitsToInclude);
	}
	
	public void setIncludeList(ArrayList<int[]> includeList) {
		this.includeList = includeList;
	}

	public void addTypeChangeMap(Map<PhenotypeAttribute, ATTRIBUTE_TYPE> typeMap) {
		typeChangeList.add(typeMap);
	}
	
	public void setTypeChangeList(ArrayList<Map<PhenotypeAttribute, ATTRIBUTE_TYPE>> typeChangeList) {
		this.typeChangeList = typeChangeList;
	}
	
	public void excludeLast(boolean exclude) {
		excludeLast = exclude;
	}
}



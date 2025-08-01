package net.maizegenetics.analysis.numericaltransform;

import java.awt.Frame;
import java.net.URL;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import javax.swing.ImageIcon;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;

import net.maizegenetics.phenotype.CategoricalAttribute;
import net.maizegenetics.phenotype.NumericAttribute;
import net.maizegenetics.phenotype.Phenotype;
import net.maizegenetics.phenotype.PhenotypeAttribute;
import net.maizegenetics.phenotype.Phenotype.ATTRIBUTE_TYPE;
import net.maizegenetics.phenotype.PhenotypeBuilder;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.util.OpenBitSet;

public class TransformDataPlugin extends AbstractPlugin {
	private static Logger myLogger = LogManager.getLogger(TransformDataPlugin.class);
	public enum BASE {natural, base_2, base_10};
	
	private List<NumericAttribute> traitsToTransform;
	private List<CategoricalAttribute> byFactor;
	private boolean logTransform = false;
	private boolean powerTransform = false;
	private boolean standardize = false;
	private BASE myBase = BASE.natural;
	private double power = 1;
	private static final double log2 = Math.log(2);
	private boolean allTraits = true;
	private String traitnames = "";
	private String factornames = "";
		
	public TransformDataPlugin(Frame parentFrame, boolean isInteractive) {
		super(parentFrame, isInteractive);
	}

	public DataSet processData(DataSet input){
            
		List<Datum> myData = input.getDataOfType(Phenotype.class);
		if (myData.size() == 1) {
			Phenotype myPhenotype = (Phenotype) myData.get(0).getData();
			if (isInteractive()) {
				allTraits = false;
				List<NumericAttribute> numericAttributes = Stream.concat(myPhenotype.attributeListOfType(ATTRIBUTE_TYPE.data).stream(), 
						myPhenotype.attributeListOfType(ATTRIBUTE_TYPE.covariate).stream())
						.map(pa -> (NumericAttribute) pa)
						.collect(Collectors.toList());
				
				List<CategoricalAttribute> catAttributes = myPhenotype.attributeListOfType(ATTRIBUTE_TYPE.factor).stream()
						.map(pa -> (CategoricalAttribute) pa)
						.collect(Collectors.toList());
				
				TransformDataDialog tdd = new TransformDataDialog(getParentFrame(), numericAttributes, catAttributes);
				tdd.setVisible(true);

				traitsToTransform = tdd.traitsToTransform();
				byFactor = tdd.factorsForStandardizing();
				logTransform = tdd.logTransformation();
				powerTransform = tdd.powerTransformation();
				standardize = tdd.standardize();
				myBase = tdd.base();
				power = tdd.exponent();
			} else {
				//add traitnames to list of attributes to be transformed
				if (traitnames.length() == 0) {
					traitsToTransform = myPhenotype.attributeListOfType(ATTRIBUTE_TYPE.data).stream()
							.map(a -> (NumericAttribute) a)
							.collect(Collectors.toList());
				} else {
					String[] attributeNames = traitnames.split(",");
					traitsToTransform = myPhenotype.attributeListOfType(ATTRIBUTE_TYPE.data).stream()
							.filter(a -> contains(a.name(), attributeNames))
							.map(a -> (NumericAttribute) a)
							.collect(Collectors.toList());
				}
				
				//add factornames ot list of factors for stanardization
				if (factornames.length() == 0) {
					byFactor = new ArrayList<>();  //do not use factors
				} else {
					String[] attributeNames = factornames.split(",");
					byFactor = myPhenotype.attributeListOfType(ATTRIBUTE_TYPE.factor).stream()
							.filter(a -> contains(a.name(), attributeNames))
							.map(a -> (CategoricalAttribute) a)
							.collect(Collectors.toList());
				}
				
			}
			
			if (logTransform || powerTransform || standardize) return transformTraits(myPhenotype, myData.get(0));
			else return null;
			
		}

		throw new IllegalArgumentException("TransformDataPlugin: Please select one Phenotype data set.");
	}

	private boolean contains(String name, String[] array) {
		return Arrays.stream(array).anyMatch(str -> str.equals(name));
	}
	
	public DataSet transformTraits(Phenotype myPhenotype, Datum myData) {
		//use a sequential stream, because the order of the attributes needs to stay the same
		List<PhenotypeAttribute> myNewAttributes = myPhenotype.attributeListCopy().stream()
				.map(a -> transformAttribute(a))
				.collect(Collectors.toList());
		
		Phenotype transformedPhenotype = new PhenotypeBuilder().fromAttributeList(myNewAttributes, myPhenotype.typeListCopy()).build().get(0);
		
		StringBuilder nameBuilder = new StringBuilder();
		nameBuilder.append("transformed_").append(myData.getName());
		
		StringBuilder commentBuilder = new StringBuilder();
		commentBuilder.append("Phenotypes transformed from ");
		commentBuilder.append(myData.getName()).append("\n");
//		commentBuilder.append(myData.getComment());
		commentBuilder.append("The following traits were transformed by ");
		if (powerTransform) commentBuilder.append("using a power ").append(power).append(" transformation:\n");
		else if (logTransform) commentBuilder.append("using a ").append(myBase.name()).append(" log transformation:\n");
		if (standardize) commentBuilder.append("standardizing.\n");
		for (NumericAttribute na : traitsToTransform) commentBuilder.append(na.name()).append("\n");
		
		return new DataSet(new Datum(nameBuilder.toString(), transformedPhenotype, commentBuilder.toString()), this);
	}
	
	public PhenotypeAttribute transformAttribute(PhenotypeAttribute myAttribute) {
		if (!(myAttribute instanceof NumericAttribute)) return myAttribute;
		NumericAttribute myNumericAttribute = (NumericAttribute) myAttribute;
		if (!traitsToTransform.contains(myNumericAttribute)) return myAttribute;
		
		if (powerTransform) myNumericAttribute = powerTransform(myNumericAttribute);
		else if (logTransform) myNumericAttribute = logTransform(myNumericAttribute);
		
		if (standardize) {
			if (byFactor.size() > 0) return standardize(myNumericAttribute, byFactor);
			return standardize(myNumericAttribute);
		}
		
		return myNumericAttribute;
	}
	
	@Override
	public ImageIcon getIcon() {
        URL imageURL = TransformDataPlugin.class.getResource("/net/maizegenetics/analysis/images/Transform.gif");
        if (imageURL == null) {
            return null;
        } else {
            return new ImageIcon(imageURL); 
        }
	}

	@Override
	public String getButtonName() {
		return "Transform Phenotype";
	}

	@Override
	public String getToolTipText() {
		return "Transform or standardize phenotypes";
	}

	public NumericAttribute powerTransform(NumericAttribute original) {
		float[] originalValues = original.floatValues();
		int n = originalValues.length;
		float[] transValues = new float[n];
		
		for (int i = 0; i < n; i++) transValues[i] = (float) Math.pow(originalValues[i], power);
		
		return new NumericAttribute(original.name(), transValues, original.missing());
	}
	
	public NumericAttribute logTransform(NumericAttribute original) {
		float[] originalValues = original.floatValues();
		int n = originalValues.length;
		float[] transValues = new float[n];
		double divisor;
		
		switch (myBase) {
		case base_10:
			divisor = Math.log(10);
			break;
		case base_2:
			divisor = Math.log(2);
			break;
		default:
			divisor = 1;
		}
		
		for (int i = 0; i < n; i++) {
			switch (myBase) {
			case natural:
				transValues[i] = (float) Math.log(originalValues[i]);
				break;
			case base_10:
			case base_2:
				transValues[i] = (float) (Math.log(originalValues[i]) / divisor);
				break;
			}
		}
		
		return new NumericAttribute(original.name(), transValues, original.missing());
	}
	
	public NumericAttribute standardize(NumericAttribute original) {
		float[] originalValues = original.floatValues();
		int n = originalValues.length;
		float[] meanSD = meanStdDev(originalValues);
		float[] transValues = new float[n];
		
		for (int i = 0; i < n; i++) transValues[i] = (originalValues[i] - meanSD[0])/meanSD[1];
		
		return new NumericAttribute(original.name(), transValues, original.missing());
	}
	
	public float[] meanStdDev(float[] data) {
		int n = data.length;
		double sum = 0;
		double sumsq = 0;
		int notMissingCount = 0;
		for (int i = 0; i < n; i++) {
			double val = data[i];
			if (!Double.isNaN(val)) {
				sum += val;
				sumsq += val * val;
				notMissingCount++;
			}
			
		}
		float mean = (float) sum / notMissingCount;
		float sdev = (float) Math.sqrt((sumsq - sum / notMissingCount * sum) / (notMissingCount - 1));
		return new float[]{mean, sdev};
	}
	
	public NumericAttribute standardize(NumericAttribute original, List<CategoricalAttribute> byFactors) {
		
		List<int[]> subsetList = subsets(byFactors);
		float[] stdData = Arrays.copyOf(original.floatValues(), original.size());
		for (int[] subset : subsetList) {
			int n = subset.length;
			float[] subsetData = new float[n];
			for (int i = 0; i < n; i++) subsetData[i] = stdData[subset[i]];
			float[] meanSD = meanStdDev(subsetData);
			for (int i = 0; i < n; i++) stdData[subset[i]] = (stdData[subset[i]] - meanSD[0]) / meanSD[1] ;
		}
		
		return new NumericAttribute(original.name(), stdData, original.missing());
	}
	
	public List<int[]> subsets(List<CategoricalAttribute> byFactors) {
		class subset {
			int[] levels;
			subset(int[] levels) { this.levels = levels; }
			public boolean equals(Object other) {
				if (other instanceof subset) return Arrays.equals(levels, ((subset) other).levels);
				return false;
			}
			
			@Override
			public int hashCode() {
				int hc = 0;
				int mult = 1;
				for (int i : levels) {
					hc += mult * Integer.hashCode(i);
					mult *= 10;
				}
				return hc;
			}
		}
		
		int nobs = byFactors.get(0).size();
		OpenBitSet missing = new OpenBitSet(nobs);
		for (PhenotypeAttribute pa : byFactors) missing.or(pa.missing());
		
		int nfactors = byFactors.size();
		Multimap<subset, Integer> subsetMap = HashMultimap.create();
		for (int obs = 0; obs < nobs; obs++) if (!missing.fastGet(obs)) {
			int[] levels = new int[nfactors];
			int count = 0;
			for (CategoricalAttribute ca : byFactors) {
				levels[count++] = ca.intValue(obs);
			}
			subsetMap.put(new subset(levels), obs);
		}
		
		ArrayList<int[]> subsetList = new ArrayList<>();
		for (subset sub : subsetMap.keySet()) {
			subsetList.add(subsetMap.get(sub)
					.stream()
					.mapToInt(Integer::intValue)
					.toArray()); 
		}
		return subsetList;
	}

	
	@Override
	public void setParameters(String[] args) {
		// TODO Auto-generated method stub
		traitsToTransform = new ArrayList<>();
		byFactor = new ArrayList<>();
		int argPtr = 0;
		while (argPtr < args.length) {
			if (args[argPtr].toLowerCase().startsWith("-trait")) {
				setTraits(args[++argPtr]);
				argPtr++;
			} else if (args[argPtr].toLowerCase().startsWith("-factor")) {
				setFactors(args[++argPtr]);
				argPtr++;
			} else if (args[argPtr].equals("-log")) {
				logTransform = true;
				powerTransform = false;
				String baseName = args[++argPtr];
				if (baseName.equals("natural")) myBase = BASE.natural;
				else if (baseName.equals("base_2")) myBase = BASE.base_2;
				else if (baseName.equals("base_10")) myBase = BASE.base_10;
				else throw new IllegalArgumentException("-log parameter value must be one of natural, base_2, base_10.");
				argPtr++;
			} else if (args[argPtr].equals("-power")) {
				powerTransform = true;
				logTransform = false;
				try {
					power = Double.parseDouble(args[++argPtr]);
				} catch(NumberFormatException nfe) {
					myLogger.error("-power parameter value must be a floating point number.",  nfe);
				}
				argPtr++;
			} else if (args[argPtr].equals("-standardize")) {
				String paramVal = args[++argPtr];
				if (paramVal.toLowerCase().startsWith("t")) standardize = true;
				else standardize = false;
				argPtr++;
			} else {
				String msg = String.format("unrecognized command line parameter for TransformDataPlugin: %s", args[argPtr]);
				myLogger.error(msg);
				throw new IllegalArgumentException(msg);
			}
		}
	}

	@Override
	public String getUsage() {
		StringBuilder usageString = new StringBuilder();
		usageString.append("The TransformDataPlugin can take the following parameters. Cannot use both log and power parameters.:\n");
		usageString.append("-traits: A comma delimited list of trait names with no embedded space. If this parameter is not specified then all traits will be transformed.");
		usageString.append("-factor: The factor name or a comma-delimited list of factor names with no embedded spaces within which values are to be standardized. ");
		usageString.append("The default is to ignore factors and use the mean and standard deviation of all observations.\n");
		usageString.append("-log: perform a log transformation. Can take one of natural, base_2, or base_10. Default = no transformation.");
		usageString.append("-power: perform a power transformation. The parameter value is the exponent to which each value should be raised. Default = no transformation.");
		usageString.append("-standardize: standardize values by subtracting the mean and dividing by the standard deviation. true or false. Default = false.");
		return usageString.toString();
	}

	public void setTraits(String namelist) {
		traitnames = namelist;
	}
	
	public void setFactors(String namelist) {
		factornames = namelist;
	}
	
	public void setTraitsToTransform(List<NumericAttribute> traitsToTransform) {
		this.traitsToTransform = traitsToTransform;
	}

	public void setByFactor(List<CategoricalAttribute> byFactor) {
		this.byFactor = byFactor;
	}

	public void setLogTransform(boolean logTransform) {
		this.logTransform = logTransform;
	}

	public void setPowerTransform(boolean powerTransform) {
		this.powerTransform = powerTransform;
	}

	public void setStandardize(boolean standardize) {
		this.standardize = standardize;
	}

	public void setMyBase(BASE myBase) {
		this.myBase = myBase;
	}

	public void setPower(double power) {
		this.power = power;
	}

}

package net.maizegenetics.analysis.numericaltransform;

import static org.junit.Assert.*;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import net.maizegenetics.analysis.numericaltransform.TransformDataPlugin;
import net.maizegenetics.phenotype.NumericAttribute;
import net.maizegenetics.phenotype.Phenotype;
import net.maizegenetics.phenotype.Phenotype.ATTRIBUTE_TYPE;
import net.maizegenetics.phenotype.CategoricalAttribute;
import net.maizegenetics.phenotype.PhenotypeAttribute;
import net.maizegenetics.phenotype.PhenotypeBuilder;
import net.maizegenetics.phenotype.TaxaAttribute;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.util.OpenBitSet;

import org.junit.After;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

public class TransformDataPluginTest {
	static Phenotype myPheno;
	
	@BeforeClass 
	public static void initializeTestData() {
		List<Taxon> taxa = new ArrayList<>();
		for (int i = 1; i <= 6; i++) taxa.add(new Taxon(String.format("t%d", i)));
		
		String[] loc = new String[] {"a","a","a","b","b","b"};
		float[] trt1 = new float[]{1,2,3,4,5,6};
		float[] trt2 = new float[] {11,12,13,14,15,16};
		List<PhenotypeAttribute> pa = new ArrayList<>();
		List<ATTRIBUTE_TYPE> attr = new ArrayList<>();
		pa.add(new TaxaAttribute(taxa));
		attr.add(ATTRIBUTE_TYPE.taxa);
		pa.add(new CategoricalAttribute("Loc", loc));
		attr.add(ATTRIBUTE_TYPE.factor);
		pa.add(new NumericAttribute("trt1", trt1, new OpenBitSet(6)));
		attr.add(ATTRIBUTE_TYPE.data);
		pa.add(new NumericAttribute("trt2", trt2, new OpenBitSet(6)));
		attr.add(ATTRIBUTE_TYPE.data);
		myPheno = new PhenotypeBuilder().fromAttributeList(pa, attr).build().get(0);
	}
	
	@Test
	public void testPowerTransformation() {
		TransformDataPlugin tdp = new TransformDataPlugin(null, false);
		tdp.setParameters(new String[] {"-traits", "trt1,trt2", "-power", "2"});
		DataSet inputData = new DataSet(new Datum("test", myPheno, "test"), null);
		DataSet resultSet = tdp.processData(inputData);
		Phenotype resultPhenotype = (Phenotype) resultSet.getData(0).getData();
		float[] resultArray = (float[]) resultPhenotype.attribute(2).allValues();
		float[] expectedArray = new float[]{1,4,9,16,25,36};
		assertArrayEquals(expectedArray, resultArray, 0.0001f);
		resultArray = (float[]) resultPhenotype.attribute(3).allValues();
		expectedArray = new float[]{121,144,169,196,225,256};
		assertArrayEquals(expectedArray, resultArray, 0.0001f);
	}
	
	@Test
	public void testLogTransformation() {
		TransformDataPlugin tdp = new TransformDataPlugin(null, false);
		tdp.setParameters(new String[] {"-traits", "trt1,trt2", "-log", "natural"});
		DataSet inputData = new DataSet(new Datum("test", myPheno, "test"), null);
		DataSet resultSet = tdp.processData(inputData);
		Phenotype resultPhenotype = (Phenotype) resultSet.getData(0).getData();
		float[] resultArray = (float[]) resultPhenotype.attribute(2).allValues();
		float[] expectedArray = new float[]{0,0.693147f,1.098612f,1.386294f,1.609437f,1.791759f};
		assertArrayEquals(expectedArray, resultArray, 0.00001f);
		resultArray = (float[]) resultPhenotype.attribute(3).allValues();
		expectedArray = new float[]{2.3978952f,2.4849066f,2.5649493f,2.6390573f,2.7080502f,2.7725887f};
		assertArrayEquals(expectedArray, resultArray, 0.00001f);
	}
	
	@Test
	public void testStandardize() {
		TransformDataPlugin tdp = new TransformDataPlugin(null, false);
		tdp.setParameters(new String[] {"-traits", "trt1,trt2", "-standardize", "true"});
		DataSet inputData = new DataSet(new Datum("test", myPheno, "test"), null);
		DataSet resultSet = tdp.processData(inputData);
		Phenotype resultPhenotype = (Phenotype) resultSet.getData(0).getData();
		float[] resultArray = (float[]) resultPhenotype.attribute(2).allValues();
		float[] expectedArray = new float[]{-1.336306f,-0.801783f,-0.267261f,0.267261f,0.801783f,1.336306f};
		assertArrayEquals(expectedArray, resultArray, 0.00001f);
		resultArray = (float[]) resultPhenotype.attribute(3).allValues();
		assertArrayEquals(expectedArray, resultArray, 0.00001f);

		tdp = new TransformDataPlugin(null, false);
		tdp.setParameters(new String[] {"-traits", "trt1", "-standardize", "true", "-factors", "Loc"});
		inputData = new DataSet(new Datum("test", myPheno, "test"), null);
		resultSet = tdp.processData(inputData);
		resultPhenotype = (Phenotype) resultSet.getData(0).getData();
		resultArray = (float[]) resultPhenotype.attribute(2).allValues();
		expectedArray = new float[]{-1,0,1,-1,0,1};
		assertArrayEquals(expectedArray, resultArray, 0.00001f);
		resultArray = (float[]) resultPhenotype.attribute(3).allValues();
		expectedArray = new float[] {11,12,13,14,15,16};
		assertArrayEquals(expectedArray, resultArray, 0.00001f);
}
	

}

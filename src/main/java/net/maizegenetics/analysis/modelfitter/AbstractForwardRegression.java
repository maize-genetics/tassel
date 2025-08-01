package net.maizegenetics.analysis.modelfitter;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.util.ArrayList;
import java.util.List;
import java.util.TreeSet;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import net.maizegenetics.analysis.association.AssociationUtils;
import net.maizegenetics.analysis.modelfitter.AdditiveSite.CRITERION;
import net.maizegenetics.dna.map.Position;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.phenotype.CategoricalAttribute;
import net.maizegenetics.phenotype.GenotypePhenotype;
import net.maizegenetics.phenotype.NumericAttribute;
import net.maizegenetics.phenotype.Phenotype;
import net.maizegenetics.phenotype.Phenotype.ATTRIBUTE_TYPE;
import net.maizegenetics.phenotype.PhenotypeAttribute;
import net.maizegenetics.phenotype.PhenotypeBuilder;
import net.maizegenetics.phenotype.TaxaAttribute;
import net.maizegenetics.stats.linearmodels.CovariateModelEffect;
import net.maizegenetics.stats.linearmodels.FactorModelEffect;
import net.maizegenetics.stats.linearmodels.ModelEffect;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListBuilder;
import net.maizegenetics.taxa.TaxaListUtils;
import net.maizegenetics.taxa.Taxon;

/**
 * @author pbradbury
 *
 */
public abstract class AbstractForwardRegression implements ForwardRegression {
    //performs forward regression
    //requires one phenotype as a double array and a genotype table
    //accepts additional fixed effects, either covariates or factors
    //returns a list of markers with p-values
    //no missing values allowed in phenotype, factors, or covariates

    private Logger myLogger = LogManager.getLogger(AbstractForwardRegression.class);
    protected double[] y;     //data for a single phenotype (no missing data allowed)
    protected final GenotypePhenotype myGenotypePhenotype;
    protected final GenotypeTable myGenotype;
    protected final Phenotype myPhenotype;
    protected double enterLimit;
    protected int maxVariants;
    protected final int numberOfSites;
    protected final int numberOfObservations;
    protected List<AdditiveSite> siteList;
    protected final List<ModelEffect> myBaseModel;
    protected List<ModelEffect> myModel;
    protected List<Object[]> myFittedVariants = new ArrayList<>();
    protected String traitname;
    protected int[] taxaIndex = null;

    /**
     * @param data              a GenotypePhenotype object
     * @param phenotypeIndex    the attribute index of the phenotype to be analyzed
     * @param baseModel         the fixed effects in the model. If null, will be set to the mean only.
     * @param enterLimit        terms will be added to the model as long as the p-value of the next term is less than or equal to enterLimit
     * @param maxVariants       at most maxVariant terms will be fit. If the enterLimit is reached first, fewer than maxVariant terms will be fit.
     */
    public AbstractForwardRegression(GenotypePhenotype data) {
        myGenotypePhenotype = data;
        myGenotype = myGenotypePhenotype.genotypeTable();
        myPhenotype = myGenotypePhenotype.phenotype();
        numberOfSites = data.genotypeTable().numberOfSites();
        numberOfObservations = data.numberOfObservations();
        myBaseModel = getBaseModel();
        myModel = new ArrayList<>(myBaseModel);

        //check for genotype vs. ref probability
       
        //Initialize the siteList
        siteList = new ArrayList<>();
        long start = System.nanoTime();
        
        if (myGenotype.hasGenotype()) {
            for (int s = 0; s < numberOfSites; s++) {
                Position pos = myGenotype.positions().get(s);
                siteList.add(new GenotypeAdditiveSite(s, pos.getChromosome().getName(), pos.getPosition(), pos.getSNPID(), CRITERION.pval, myGenotypePhenotype.genotypeAllTaxa(s), myGenotype.majorAllele(s), myGenotype.majorAlleleFrequency(s)));
            }
        	
        } else if (myGenotype.hasReferenceProbablity()) {
            for (int s = 0; s < numberOfSites; s++) {
                Position pos = myGenotype.positions().get(s);
                siteList.add(new RefProbAdditiveSite(s, pos.getChromosome().getName(), pos.getPosition(), pos.getSNPID(),CRITERION.pval, myGenotypePhenotype.referenceProb(s)));
            }
        	
        } else {
        		throw new IllegalArgumentException("Input has neither genotype nor reference probability.");
        }
        
        myLogger.debug(String.format("site list created with %d sites in %d ms.", siteList.size(), (System.nanoTime() - start) / 1000000));

    }

    /**
     * @param serialFilename    the base name of the serialization files that end in _taxa.bin and _sites.bin
     * @param pheno     the Phenotype object to be used in the analysis
     */
    public AbstractForwardRegression(String serialFilename, Phenotype pheno) {
        List<AdditiveSite> mySites = new ArrayList<>();
        TaxaListBuilder taxaBuilder = new TaxaListBuilder();
        int ntaxa = 0, nsites = 0;
        
        try (FileInputStream fos = new FileInputStream(serialFilename + "_taxa.bin")) {
            ObjectInputStream input = new ObjectInputStream(fos);
            ntaxa = ((Integer) input.readObject()).intValue();
            nsites = ((Integer) input.readObject()).intValue();
            
            for (int t = 0; t < ntaxa; t++) {
                taxaBuilder.add(new Taxon((String) input.readObject()));
            }
            input.close();
        } catch (FileNotFoundException e) {
            throw new RuntimeException("Unable to open taxa.bin input file.", e);
        } catch (IOException e) {
            throw new RuntimeException("error deserializing taxa: ", e);
        } catch (ClassNotFoundException e) {
            throw new RuntimeException("error deserializing taxa: ", e);
        }
        
        try (FileInputStream fos = new FileInputStream(serialFilename + "_sites.bin")) {
            ObjectInputStream input = new ObjectInputStream(fos);
            for (int s = 0; s < nsites; s++) {
                mySites.add((AdditiveSite) input.readObject());
            }
            input.close();
        } catch (FileNotFoundException e) {
            throw new RuntimeException("Unable to open sites.bin input file.", e);
        } catch (IOException e) {
            throw new RuntimeException("error deserializing sites: ", e);
        } catch (ClassNotFoundException e) {
            throw new RuntimeException("error deserializing sites: ", e);
        }

        myGenotypePhenotype = null;
        myGenotype = null;
        siteList = mySites;
        numberOfSites = siteList.size();
        TaxaList siteTaxaList = taxaBuilder.build();
        TaxaList jointTaxaList = TaxaListUtils.getCommonTaxa(siteTaxaList, pheno.taxa());
        myLogger.debug(String.format("The joint taxa list has %d taxa.", jointTaxaList.size()));

        //delete any phenotypes not in the site list
        myPhenotype =
                new PhenotypeBuilder().fromPhenotype(pheno).keepTaxa(jointTaxaList).build().get(0);
        numberOfObservations = myPhenotype.numberOfObservations();
        myBaseModel = getBaseModel();

        //create an index from siteList into phenotype
        TaxaAttribute myTaxa = myPhenotype.taxaAttribute();
        int[] taxaIndex =
                myTaxa.allTaxaAsList().stream().mapToInt(t -> siteTaxaList.indexOf(t)).toArray();
        TreeSet<Integer> taxaSet = new TreeSet<>();
        for (int t : taxaIndex)
            taxaSet.add(t);
        List<Integer> uniqueTaxa = new ArrayList<>(taxaSet);

        long countOfNegativeIndices = uniqueTaxa.stream().filter(I -> I < 0).count();
        myLogger.debug(String.format("siteIndices has %d negative values.", countOfNegativeIndices));

        for (AdditiveSite as : mySites)
            as.reindexTaxa(taxaIndex, uniqueTaxa);
        myLogger.debug("sites reindexed.");
    }

    @Override
    public void resetModel(int phenotypeIndex, double enterLimit, int maxVariants) {
        PhenotypeAttribute myTraitAttribute =
                myPhenotype.attributeListOfType(ATTRIBUTE_TYPE.data).get(phenotypeIndex);
        traitname = myTraitAttribute.name();
        myModel = new ArrayList<>(myBaseModel);
        y = AssociationUtils.convertFloatArrayToDouble((float[]) myTraitAttribute.allValues());

        this.enterLimit = enterLimit;
        this.maxVariants = maxVariants;
    }

    protected List<ModelEffect> getBaseModel() {
        List<ModelEffect> base = new ArrayList<>();
        int[] mean = new int[numberOfObservations];
        ModelEffect meanEffect = new FactorModelEffect(mean, false, "mean");
        base.add(meanEffect);
        for (PhenotypeAttribute factor : myPhenotype.attributeListOfType(ATTRIBUTE_TYPE.factor)) {
            ModelEffect factorEffect =
                    new FactorModelEffect(((CategoricalAttribute) factor).allIntValues(), true, factor.name());
            base.add(factorEffect);
        }
        for (PhenotypeAttribute cov : myPhenotype.attributeListOfType(ATTRIBUTE_TYPE.covariate)) {
            double[] values =
                    AssociationUtils.convertFloatArrayToDouble(((NumericAttribute) cov).floatValues());
            ModelEffect covEffect = new CovariateModelEffect(values, cov.name());
            base.add(covEffect);
        }
        return base;
    }

    @Override
    public List<Object[]> fittedModel() {
        return myFittedVariants;
    }

    @Override
    public Phenotype phenotype() {
        return myPhenotype;
    }

    protected void addVariant(AdditiveSite site, double p, int iteration, int step) {
        myFittedVariants.add(new Object[] { traitname, iteration, step, site.siteName(),
                site.chromosomeName(), site.position(), p, -Math.log10(p) });
    }

}

package net.maizegenetics.analysis.modelfitter;

/**
 * This class holds data for a variant site to provide an object to hold data used and returned by a stepwise
 * model fitter. That makes it convenient to use in a parallel stream for processing.
 */
public abstract class AbstractAdditiveSite implements AdditiveSite {
    private static final long serialVersionUID = 3032879930663240377L;
    protected final int siteIndex;
    protected final String chrName;
    protected final int position;
    protected final String name;
    protected double criterionValue;
    protected final CRITERION selectionCriterion;
    protected final int direction;

    /**
     *
     * @param site                  a zero based index of sites
     * @param chromosomeName        a chromosome name (String)
     * @param pos                   the integer position in the chromosome
     * @param id                    a name for this site
     * @param selectionCriterion    the selection criterion to be use to pick a site
     */
    public AbstractAdditiveSite(int site, String chromosomeName, int pos, String id, CRITERION selectionCriterion) {
        siteIndex = site;
        this.selectionCriterion = selectionCriterion;
        if (selectionCriterion == CRITERION.pval)
            direction = 1;
        else
            direction = -1;
        chrName = chromosomeName;
        position = pos;
        name = id;
    }

    @Override
    public int siteNumber() {
        return siteIndex;
    }

    @Override
    public String chromosomeName() {
        return chrName;
    }

    @Override
    public int position() {
        return position;
    }

    @Override
    public String siteName() {
        return name;
    }

    @Override
    public double criterionValue() {
        return criterionValue;
    }

    @Override
    public void criterionValue(double value) {
        criterionValue = value;
    }

    @Override
    public CRITERION selectionCriterion() {
        return selectionCriterion;
    }

    @Override
    public int compareTo(AdditiveSite other) {
        return direction * Double.compare(criterionValue, other.criterionValue());
    }

}

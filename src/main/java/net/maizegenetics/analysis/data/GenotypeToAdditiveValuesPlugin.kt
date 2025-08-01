package net.maizegenetics.analysis.data

import net.maizegenetics.dna.snp.GenotypeTable
import net.maizegenetics.dna.snp.GenotypeTableUtils
import net.maizegenetics.dna.snp.genotypecall.AlleleFreqCache
import net.maizegenetics.plugindef.AbstractPlugin
import net.maizegenetics.plugindef.DataSet
import net.maizegenetics.plugindef.Datum
import net.maizegenetics.util.ColumnMatrix
import java.awt.Frame
import javax.swing.ImageIcon

/**
 * @author Terry Casstevens
 * Created December 06, 2018
 *
 * This plugin takes a genotype matrix and creates a matrix of numbers.
 * The number values are the count of alleles that do not match the major allele
 * at each site / taxon.
 */

class GenotypeToAdditiveValuesPlugin(parentFrame: Frame?, isInteractive: Boolean) : AbstractPlugin(parentFrame, isInteractive) {

    override fun processData(input: DataSet?): DataSet {

        val temp = input?.getDataOfType(GenotypeTable::class.java)
        if (temp?.size != 1) {
            throw IllegalArgumentException("GenotypeToAdditiveValuesPlugin: processData: must input a genotype")
        }
        val genotype = temp[0].data as GenotypeTable

        val result = ColumnMatrix.Builder(genotype.numberOfTaxa(), genotype.numberOfSites())

        for (site in 0 until genotype.numberOfSites()) {
            val siteGenotypes = genotype.genotypeAllTaxa(site)
            val alleleCounts = AlleleFreqCache.allelesSortedByFrequencyNucleotide(siteGenotypes)
            val majorAllele = AlleleFreqCache.majorAllele(alleleCounts)

            // value assigned to site / taxon is the number of alleles
            // that do not match the major allele.
            for (taxon in 0 until genotype.numberOfTaxa()) {
                var value: Byte = 0
                val alleles = GenotypeTableUtils.getDiploidValues(siteGenotypes[taxon])
                if (alleles[0] != majorAllele) value++
                if (alleles[1] != majorAllele) value++
                result.set(taxon, site, value)
            }
        }

        return DataSet(Datum("${temp[0].name} Additive Values", result.build(), null), this)

    }

    fun runPlugin(genotype: GenotypeTable): ColumnMatrix {
        return performFunction(DataSet.getDataSet(genotype)).getDataOfType(ColumnMatrix::class.java)[0].data as ColumnMatrix
    }

    override fun getIcon(): ImageIcon? {
        return null
    }

    override fun getButtonName(): String {
        return "Convert Genotype to Additive Values"
    }

    override fun getToolTipText(): String {
        return "Convert Genotype to Additive Values"
    }

}

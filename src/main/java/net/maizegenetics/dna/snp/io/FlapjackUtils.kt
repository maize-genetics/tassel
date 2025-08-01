@file:JvmName("FlapjackUtils")

package net.maizegenetics.dna.snp.io

import net.maizegenetics.dna.snp.GenotypeTable
import net.maizegenetics.dna.snp.NucleotideAlignmentConstants.UNDEFINED_ALLELE_STR
import net.maizegenetics.util.Utils
import org.apache.logging.log4j.LogManager


/**
 * @author Terry Casstevens
 * Created July 02, 2019
 *
 * This specifics of the Flapjack format can be found here
 * http://flapjack.hutton.ac.uk/en/latest/projects_&_data_formats.html
 */

private val myLogger = LogManager.getLogger("net.maizegenetics.dna.snp.io.FlapjackUtils")

private const val FLAPJACK_MISSING = "-"

// Translation from TASSEL byte encoding to Flapjack character.
// Differences from TASSEL.  N (Unknown), - (Deletion), and + (Insertion)
// are converted to FLAPJACK_MISSING (-)
// UNDEFINED_ALLELE_STR should never end up in the file.
private val FLAPJACK_CHARS = arrayOf("A", "C", "G", "T", FLAPJACK_MISSING, FLAPJACK_MISSING,
        UNDEFINED_ALLELE_STR, UNDEFINED_ALLELE_STR, UNDEFINED_ALLELE_STR, UNDEFINED_ALLELE_STR,
        UNDEFINED_ALLELE_STR, UNDEFINED_ALLELE_STR, UNDEFINED_ALLELE_STR, UNDEFINED_ALLELE_STR,
        UNDEFINED_ALLELE_STR, FLAPJACK_MISSING)


/**
 * Writes given genotype table to Flapjack format
 */
fun writeToFlapjack(genotypes: GenotypeTable, filename: String, delimiter: Char = '\t'): String {

    if (delimiter != ' ' && delimiter != '\t') {
        throw IllegalArgumentException("FlapjackUtils: writeToFlapjack: Delimiter character must be either a blank space or a tab.")
    }

    val mapFileName = Utils.addSuffixIfNeeded(filename, ".flpjk.map")
    val genoFileName = Utils.addSuffixIfNeeded(filename, ".flpjk.geno")

    try {

        Utils.getBufferedWriter(genoFileName).use { genoWriter ->

            Utils.getBufferedWriter(mapFileName).use { mapWriter ->

                for (site in 0 until genotypes.numberOfSites()) {
                    mapWriter.write(genotypes.siteName(site))
                    mapWriter.write(delimiter.toInt())
                    mapWriter.write(genotypes.chromosomeName(site))
                    mapWriter.write(delimiter.toInt())
                    mapWriter.write(Integer.toString(genotypes.chromosomalPosition(site)))
                    mapWriter.write("\n")
                    genoWriter.write(delimiter.toInt())
                    genoWriter.write(genotypes.siteName(site))
                }

            }

            genoWriter.write("\n")

            for (taxa in 0 until genotypes.numberOfTaxa()) {

                genoWriter.write(genotypes.taxaName(taxa))

                for (site in 0 until genotypes.numberOfSites()) {

                    genoWriter.write(delimiter.toInt())

                    val alleles = genotypes.genotypeArray(taxa, site)
                    when (alleles.size) {
                        1 -> {
                            genoWriter.write(FLAPJACK_CHARS[alleles[0].toInt()])
                        }
                        2 -> {
                            genoWriter.write(FLAPJACK_CHARS[alleles[0].toInt()])
                            genoWriter.write("/")
                            genoWriter.write(FLAPJACK_CHARS[alleles[1].toInt()])
                        }
                        else -> {
                            genoWriter.write(FLAPJACK_CHARS[alleles[0].toInt()])
                            genoWriter.write("/")
                            genoWriter.write(FLAPJACK_CHARS[alleles[1].toInt()])
                            genoWriter.write("/")
                            genoWriter.write(FLAPJACK_CHARS[alleles[2].toInt()])
                        }
                    }

                }

                genoWriter.write("\n")
            }
        }

        return "$mapFileName and $genoFileName"

    } catch (e: Exception) {
        myLogger.debug(e.message, e)
        throw IllegalStateException("Problem writing Flapjack files: $mapFileName and $genoFileName\n${e.message}")
    }

}
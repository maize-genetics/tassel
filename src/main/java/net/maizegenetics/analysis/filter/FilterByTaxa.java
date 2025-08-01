/*
 *  FilterByTaxa
 *
 *  Created on Dec 21, 2016
 */
package net.maizegenetics.analysis.filter;

import net.maizegenetics.dna.snp.FilterGenotypeTable;
import net.maizegenetics.dna.snp.FilterTaxa;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.genotypecall.ListStats;
import net.maizegenetics.dna.snp.genotypecall.Stats;
import net.maizegenetics.taxa.TaxaList;

import java.util.stream.IntStream;
import java.util.stream.Stream;

/**
 * @author Terry Casstevens
 */
public class FilterByTaxa {

    private FilterByTaxa() {
    }

    public static GenotypeTable filter(GenotypeTable orig, FilterTaxa filter) {

        GenotypeTable result = orig;

        TaxaList taxaList = filter.taxaList();
        if (taxaList != null) {
            if (filter.includeTaxa()) {
                result = FilterGenotypeTable.getInstance(result, taxaList, false);
            } else {
                result = FilterGenotypeTable.getInstanceRemoveIDs(result, taxaList);
            }
        }

        TaxaList taxa = result.taxa();

        Stream<Stats> stream = null;

        if (filter.minNotMissing() != 0.0) {
            stream = stream(result, stream);
            stream = stream.filter((Stats stats) -> {
                return stats.percentNotMissing() >= filter.minNotMissing();
            });
        }

        if (filter.minHeterozygous() != 0.0 || filter.maxHeterozygous() != 1.0) {
            stream = stream(result, stream);
            stream = stream.filter((Stats stats) -> {
                double hetFreq = stats.proportionHeterozygous();
                return filter.minHeterozygous() <= hetFreq && filter.maxHeterozygous() >= hetFreq;
            });
        }

        if (stream != null) {

            TaxaList keepTaxa = stream.map((Stats stats) -> taxa.get(stats.index()))
                    .collect(TaxaList.collect());
            result = FilterGenotypeTable.getInstance(result, keepTaxa, false);

        }

        return result;

    }

    private static Stream<Stats> stream(GenotypeTable orig, Stream<Stats> stream) {
        if (stream != null) {
            return stream;
        }
        ListStats taxaStats = ListStats.getTaxaInstance(orig.genotypeMatrix());
        return IntStream.range(0, orig.numberOfTaxa()).parallel().mapToObj((int value) -> taxaStats.get(value));
    }

}

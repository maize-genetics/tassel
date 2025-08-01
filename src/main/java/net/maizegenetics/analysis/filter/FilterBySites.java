/*
 *  FilterBySites
 *
 *  Created on Dec 5, 2016
 */
package net.maizegenetics.analysis.filter;

import net.maizegenetics.dna.map.PositionList;
import net.maizegenetics.dna.snp.FilterGenotypeTable;
import net.maizegenetics.dna.snp.FilterSite;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableBuilder;
import net.maizegenetics.dna.snp.genotypecall.ListStats;
import net.maizegenetics.dna.snp.genotypecall.Stats;

import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

import static net.maizegenetics.dna.snp.GenotypeTableUtils.filterSitesByBedFile;
import static net.maizegenetics.dna.snp.GenotypeTableUtils.filterSitesByChrPos;

/**
 * @author Terry Casstevens
 */
public class FilterBySites {

    private FilterBySites() {
    }

    public static GenotypeTable filter(GenotypeTable orig, FilterSite filter) {

        GenotypeTable result = orig;

        int numSites = orig.numberOfSites();
        int numTaxa = orig.numberOfTaxa();

        if (filter.siteFilterType() == FilterSite.SITE_RANGE_FILTER_TYPES.SITES) {
            int start = filter.startSite();
            int end = filter.endSite();
            if ((start < 0) || (start > end) || (end >= orig.numberOfSites())) {
                throw new IllegalArgumentException("GenotypeTableUtils: filter: start: " + start + " or end: " + end + " site outside acceptable range.");
            }
            //result = FilterGenotypeTable.getInstance(result, start, end, filter.includeSites());
            result = FilterGenotypeTable.getInstance(result, start, end);
        } else if (filter.siteFilterType() == FilterSite.SITE_RANGE_FILTER_TYPES.POSITIONS) {

            int start = 0;
            if (filter.startPos() == -1) {
                start = orig.firstLastSiteOfChromosome(filter.startChr())[0];
            } else {
                start = orig.siteOfPhysicalPosition(filter.startPos(), filter.startChr());
                if (start < 0) {
                    start = -(start + 1);
                    if (start >= numSites) {
                        start = numSites - 1;
                    }
                }
            }

            int end = 0;
            if (filter.endPos() == -1) {
                end = orig.firstLastSiteOfChromosome(filter.endChr())[1];
            } else {
                end = orig.siteOfPhysicalPosition(filter.endPos(), filter.endChr());
                if (end < 0) {
                    end = -end - 2;
                    if (end >= numSites) {
                        end = numSites - 1;
                    }
                }
            }

            if ((start < 0) || (start > end) || (end > numSites)) {
                throw new IllegalArgumentException("GenotypeTableUtils: filter: start: " + start + " or end: " + end + " site outside acceptable range.");
            }
            result = FilterGenotypeTable.getInstance(result, start, end);
        } else if (filter.siteFilterType() == FilterSite.SITE_RANGE_FILTER_TYPES.NONE) {
            // do nothing
        } else {
            throw new IllegalStateException("GenotypeTableUtils: filter: unknown SITE_RANGE_FILTER_TYPE: " + filter.siteFilterType());
        }

        List<String> siteNames = filter.siteNames();
        if (siteNames != null) {
            if (filter.includeSites()) {
                result = FilterGenotypeTable.getInstance(result, siteNames);
            } else {
                result = FilterGenotypeTable.getInstanceRemoveSiteNames(result, siteNames);
            }
        }

        PositionList positions = filter.positionList();
        if (positions != null) {
            result = filterSitesByChrPos(result, positions, filter.includeSites());
        }

        String bedFile = filter.bedFile();
        if ((bedFile != null) && (!bedFile.isEmpty())) {
            result = filterSitesByBedFile(result, bedFile, filter.includeSites());
        }

        String chrPosFile = filter.chrPosFile();
        if ((chrPosFile != null) && (!chrPosFile.isEmpty())) {
            result = filterSitesByChrPos(result, chrPosFile, filter.includeSites());
        }

        Stream<Stats> stream = null;

        if (filter.removeSitesWithIndels()) {
            stream = stream(result, stream);
            stream = stream.filter(stats -> !stats.hasIndel());
        }

        if (filter.siteMinCount() != 0) {
            stream = stream(result, stream);
            stream = stream.filter((Stats stats) -> {
                int totalNonMissing = stats.totalGametesNonMissingForSite();
                return totalNonMissing >= (filter.siteMinCount() * 2);
            });
        }

        if (filter.siteMinAlleleFreq() != 0.0 || filter.siteMaxAlleleFreq() != 1.0) {
            stream = stream(result, stream);
            stream = stream.filter((Stats stats) -> {
                double maf = stats.minorAlleleFrequency();
                return filter.siteMinAlleleFreq() <= maf && filter.siteMaxAlleleFreq() >= maf;
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

            List<Integer> sitesToKeep = stream.map((Stats stats) -> stats.index())
                    .collect(Collectors.toList());

            int[] sites = new int[sitesToKeep.size()];
            for (int i = 0; i < sitesToKeep.size(); i++) {
                sites[i] = sitesToKeep.get(i);
            }

            result = FilterGenotypeTable.getInstance(result, sites);

        }

        if (filter.removeMinorSNPStates()) {
            result = GenotypeTableBuilder.getInstanceOnlyMajorMinor(result);
        }

        return result;

    }

    private static Stream<Stats> stream(GenotypeTable genotypes, Stream<Stats> stream) {
        if (stream != null) {
            return stream;
        }
        ListStats siteStats = ListStats.getSiteInstance(genotypes.genotypeMatrix());
        IntStream intStream = IntStream.range(0, genotypes.numberOfSites()).parallel();
        return intStream.mapToObj((int value) -> siteStats.get(value));
    }

}

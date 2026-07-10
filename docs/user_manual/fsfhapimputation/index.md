# FSFHap Imputation

FSFHap can be used to impute correct SNP calls and impute missing data in full sib families (bi-parental families). It was designed for families with inbred (or mostly inbred) parents and progeny that are at least partially inbred. It will not work at all with F1 progeny of heterozygous parents. In addition, it was designed to address issues specific to GBS data, which can have a high rate of missing data and often calls a high percentage of heterozygous sites as homozygous.

FSFHap proceeds by first attempting to identify the parent haplotypes from the progeny. It uses this strategy because when parent sequence is available, the plants that were sequenced are usually not the same ones that were used to make the full sib population and may have somewhat different genotypes. Also, the progeny taken together generally provide very high coverage for the parental sequence and are often a better source for inferring the parent sequence than the parents themselves. FSFHap has a few different algorithms for inferring the parent haplotype. The default is the most robust and should be tried first for populations derived from F1's. The "cluster" option is the method used for the analysis in the FSFHap publication, but only works well for inbred progeny that inbred to a level equivalent to four or more generations of selfing (F >= .87). The windowLD method works well for highly inbred parents, but has problems if one of the parents has residual heterozygosity. Additional options exist for backcross populations.

To infer the parent haplotypes, FSFHap analyzes one chromosome at a time. It examines windows of 50 sites by default but can be set to use a different window size. The overlap option used by the default method determines the amount of overlap between windows and is set at 25 by default. FSFHap determines the haplotypes in a window by clustering together all the sequences that have no differences. How well this works depends on the marker density, the amount of missing data, and the window size. The default method uses the overlapping sites between windows to keep track of which haplotypes come from which parent. The windowLd and cluster options use LD between windows that do not overlap to do the same thing.

Once the parent haplotypes have been inferred, each polymorphic, non-missing site in each of the progeny is coded based on the parent of origin. The Viterbi algorithm, which  uses a hidden Markov model, is used to find the most likely underlying genotype given the observed data.

The results can be output in a few different ways. After the Viterbi algorithm is run, sites are coded as A (parent 1), C (parent 2), or M (heterozygous). That can be output as parent calls with only the non-missing sites imputed. This would be a good format to use for some QTL mapping software. If the "fill gaps" option is used, any sites with missing data bordered by A calls on both sides are imputed to A. Sites border by C calls are imputed to C, and sites bordered by M calls are imputed to M.    Remaining missing sites border by different calls are left missing. To generate the "nucleotide" output, sites that are clearly monomorphic are imputed to that value in all progeny. To be called monomorphic, there must be at least two data points present in both parental haplotypes. Some sites have several non-missing data points for one of the parent haplotypes, but are entirely missing in the other parent haplotype. These are imputed in one parent and left missing in the other. Sites with too few data points to be reliably called as monomorphic are left unimputed. The non-missing data for those sites is reported unchanged.

To run FSFHap from the TASSEL GUI, select a genotype data set and choose Impute/Impute By FSFHap from the menu. The following dialog appears:

![Screen Shot 2014-12-02 at 10.56.12 AM.png](https://bitbucket.org/repo/bRpKEG/images/583650269-Screen%20Shot%202014-12-02%20at%2010.56.12%20AM.png)

The Pedigrees must be filled in by clicking Browse and navigating to a pedigree file, the format of which is described later. It is recommended that a Logfile name also be supplied for messages that will be written about the results. If the filename provided does not already exist it will be created. If it exists, new messages will be appended.

The check boxes labeled Cluster, Window LD, Bc, and Multbc control which algorithm is used to infer parent haplotypes. The default settings use the most robust algorithm for F1 populations and should be tried first. Cluster uses the algorithm described in the FSFHap manuscript (Swarts et al. 2014). It is only suitable for families of RILs inbred to the equivalent of 4 generations or more of selfing. Window LD works well on families with completely homozygous parents, but does not properly infer haplotypes in windows where one of the parents is heterozygous because it remembers exactly two haplotypes in each window. In contrast, the default algorithm keeps track of all haplotypes that it detects and tracks which parent carries each haplotype. The Bc option simply assigns the majority haplotype in each window to the recurrent parent. It is only used if the parental contribution of one of the parents is 0.75 or greater. As a result, it can be used with pedigree files that contain both F1 and BC families. Multbc was designed for multiple backcross families.

Descriptions of the other options are available by clicking the help button on the dialog or by using the command line "run_tassel.pl -FSFHapImputationPlugin".

Reference:
Novel Methods to Optimize Genotypic Imputation for Low-Coverage, Next-Generation Sequence Data in Crop Plants
Kelly Swarts, Huihui Li, J. Alberto Romero Navarro, Dong An, Maria Cinta Romay, Sarah Hearne, Charlotte Acharya, Jeffrey C. Glaubitz, Sharon Mitchell, Robert J. Elshire, Edward S. Buckler and Peter J. Bradbury
The Plant Genome  2014 7: 3: doi:10.3835/plantgenome2014.05.0023.

#Pedigree File Format

|Family|Name|Parent1|Parent 2|Contribution1|Contribution2|F|
|----|----|----|----|----|----|----|
|NAM_B97    |Z001E0001:628NHAAXX:1:250021125    |B73    |B97    |0.5    |0.5    |0.9
|NAM_B97    |Z001E0002:628NHAAXX:1:250021137    |B73    |B97    |0.5    |0.5    |0.9
|NAM_B97    |Z001E0003:628NHAAXX:1:250021149    |B73    |B97    |0.5    |0.5    |0.9

The first column contains the family name, which must be identical for each member of the family, including parents if they are included. Name is the name of individual (taxon). Parent1 and Parent2 are the parents of the population. The parent names are not currently used by the analysis, so how they are entered is unimportant. Contribution1 is the percent contribution of parent1. Examples of that value are 0.5 for a cross and 0.75 for the recurrent parent of a backcross. F is the average estimated inbreeding coefficient for the family. Contribution and inbreeding coefficient are only read from the first entry for a family and are assumed to be the same for all family members. Values should be tab separated.

#Command Line

FSFHap can be run from the command line using the FSFHapImputationPlugin. Using the plugin command alone with no options and no -endPlugin yields the following description:

FSFHapImputationPlugin Description...
The FSFHapImputation Plugin infers parental haplotypes for a full sib family then uses those haplotypes in an HMM to impute variants. It is effective at correctly imputing heterzygotes in GBS data. To use from the command line, use TASSEL's default syntax that passes data from one plugin to another (Note that this creates 2 files, one of just parental calls (A/C) and one of imputed genotypes):

	run_pipeline.pl -h input.hmp.txt -FSFHapImputationPlugin [options] -endPLugin -export output.hmp.txt

Usage:
FSFHapImputationPlugin <options>

-pedigrees <Pedigrees> : the pedigree file name (required)
-logfile <Logfile> : the name of a log file for runtime messages
-cluster <true | false> : use the cluster algorithm (Default: false)
-windowLD <true | false> : use the windowLD algorithm (Default: false)
-bc <true | false> : use the single backcross algorithm (Default: true)
-multbc <true | false> : use the multiple backcross algorithm (Default: false)
-minMaf <Min Maf> : filter out sites with less than minimumMinorAlleleFrequency [0.0‥1.0] (Default: 0.1)
-window <Window> : filter out sites with less than minimumMinorAlleleFrequency (Default: 50)
-minR <Min R> : filter out sites not correlated with neighboring sites [0.0‥1.0] (Default: 0.2)
-maxMissing <Max Missing> : filter out sites with proportion missing > maxMissing [0.0‥1.0] (Default: 0.8)
-nohets <true | false> : delete heterozygous calls before imputing (Default: false)
-maxDiff <Max Diff> : use to decide if two haplotypes are equivalent (Default: 0)
-minHap <Min Hap> : haplotype must be observed at least this often (Default: 5)
-overlap <Overlap> : overlap between adjacent windows (Default: 25)
-fillgaps <true | false> : replace missing values with flanking values if equal (Default: false)
-phet <Phet> : proportion of sites that are heterozygous [0.0‥1.0] (Default: 0.07)
-merge <true | false> : merge families and chromosomes (Default: false)
-outParents <true | false> : replace missing values with flanking values if equal (Default: true)
-outNuc <true | false> : replace missing values with flanking values if equal (Default: true)
-outIUPAC <true | false> : use IUPAC ambiguity codes for output (Default: true)

#FAQs

*Is there any way not to remove putative indels (all the data is missing for one parent) when running FSFHapImputationPlugin?*

Yes. Sort of. FSFHap can produce two types of output, parent calls and nucleotides.  Each site in the "parent calls" data set is coded A or C depending on which parent the allele came from. FSFHap only uses polymorphic sites where both parents have nucleotide calls. Where the calls from one parent are entirely missing and the other parent is G, for example, you know that G alleles came from one specific parent. But, when the value is missing the allele could have come from either parent. That makes that site less useful for making the initial parent calls. It might be possible to incorporate those sites into the algorithm but it  would not be a simple task.

The nucleotide data set includes all of the sites in the original data, not just the polymorphic ones that were included in the parent call data. If a site is monomorphic and individuals carrying haplotypes from both parents are scored at that site, then all individuals are imputed to be carrying that nucleotide at that site. On the other hand, if there are a sufficient number of alleles are scored for one of the parent haplotypes but none for the other, then individuals carrying a haplotype from the entirely missing parent are left missing and individuals carrying the haplotype from the other parent are imputed to the nucleotide at that site. Monomorphic sites with low coverage are not imputed as the completely missing parent could be missing just by chance. As a result, one can see patterns of missing data in the nucleotide data that indicate putative deletions. However, one of the alleles might be missing in one parent, not because of a deletion, but because a restriction site was missing in one parent or because the sequence in one parent did not align to reference because it was too different.

*What should I do if I get the following error: ERROR net.maizegenetics.plugindef.AbstractPlugin - Unable to find start window with only two haplotypes.*

The first thing FSFHap does to impute markers is to scan the chromosome for a window that has only two distinct haplotypes. It will only find this if both parents are homozygous in that window. If the parents are heterozygous, it will most likely generate this error message. One way to check this is to run Analysis/Geno Summary on sites then graph minor allele frequency vs. position. The points should cluster at 0.5 with scatter below that. If there is also a band at 0.25, that means one of the parents is heterozygous.

Another possible cause is that the marker density is low enough that the default window setting is too large. If that is the case, try window = 30 and overlap = 15. Another thing to try is to set the maximum missing to a lower number, say 0.7 instead of the default of 0.8.

Something else to check is whether all markers are in LD with neighboring markers by running the Tassel LD analysis. If a number of markers have been incorrectly mapped, they will show up as being in poor LD with their neighbors and, if there are enough, will interfere with the process of finding haplotypes.

*Can I include parents of a family in the pedigree file?*

That is a good idea because that allows you to check the imputation quality and determine if the genotyped lines were actually the parents of the family. However, if one of the parents is the first family member in the pedigree file and the parent contribution is entered as 1.0 and 0.0 that will generate an error. The reason is that the parent contribution for the first family member is applied to the entire family. Instead, use the same parent contribution values for the parents as for the progeny. Values of 1.0 and 0.0 may seem more correct but in this context they are not. You want to act like the parents are RILs that just happened to inherit an entire chromosome from one of the parents, which is not unusual for individual chromosomes. As the chromosomes are analyzed independently, including the parents this way creates no problem.

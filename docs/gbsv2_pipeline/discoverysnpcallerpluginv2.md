DiscoverySNPCallerPluginV2 takes a GBSv2 database file as input and identifies SNPs from the aligned tags. Tags positioned at the same physical location are aligned against one another, SNPs are called from the aligned tags, and the SNP position and allele data are written to the database.

The parameters to this plugin are:

* 	-callBiSNPsWGap <true | false> : Include sites where the third allele is a GAP (mutually exclusive with inclGaps) (Default: false)  This option has not yet been implemented.
*	-db <Input GBS Database> : Input Database file if using SQLite (REQUIRED)
*      -gapAlignRatio <Gap Alignment Threshold> : Gap alignment threshold ratio of indel contrasts to non indel contrasts: IC/(IC + NC).  Tags will be excluded from any loci that has a tag with a gap ratio that exceeds the threshold.
*	-inclGaps <true | false > : Include sites where major or minor allele is a GAP (Default: false) This option has not yet been implemented.
*	-inclRare <true | false> : Include the rare alleles at site (3 or 4th states) (Default: false)  This option has not yet been implemented.
*      -maxTagsCutSite <Max tags per cut site position> : Maximum number of tags allowed per cut site when aligning tags .  All cut site positions and their tags are stored in the database, but alignment of tags at a particular cut site position only occurs when the number of tags at this position is equal to or less than maxTagsPerCutSite.  This guards against software degradation when a position has hundreds or thousands of associated tags.  (Default: 64)
*	-mnLCov <Min Locus Coverage> : Minimum locus coverage (proportion of Taxa with a genotype) (Default: 0.1)
*	-mnMAF <Min Minor Allele Freq> : Minimum minor allele frequency (Default: 0.01)
*	-ref <Reference Genome File> : Path to reference genome in fasta format.  Ensures that a tag from the reference genome is always included when the tags at a locus are aligned against each other to call SNPs.  (Default: Don’t use reference genome)
*	-sC <Start Chromosome> : Start Chromosome:  If missing, processing starts with the first chromosome (lexicographically) in the database.
*	-eC <End Chromosome> : End Chromosome :  If missing, plugin processing ends with the last chromosome (lexicographically) in the database.
*       -deleteOldData <true | false> : Whether to delete old SNP data from the data bases.  If true, all data base tables previously populated from the DiscoverySNPCallerPluginV2 and later steps in the GBSv2 pipeline is deleted.  This allows for calling new SNPs with different pipeline parameters.  (Default: true)

DiscoverySNPCallerPluginV2 now takes string values for start and end chromosome.  The software translates the chromosome parameter value to upper case, strips off a leading "CHR" or "CHROMOSOME", and stores the remaining value up to the first space as the chromosome name.  Chromosomes that parse to an integer are compared as ints when determining order.  Otherwise chromosomes are compared lexicographically.

NOTE: When chromosomes were initially supported as string values, the software took any leading numeric portion of the chromosome, ordered based on the numeric portion, then ordered the chromosomes based on the remaining string.  In that scenario, "1A" was be less than "11C".  It was found this ordering slowed processing below an acceptable level so the software now uses strictly lexicographic ordering and "11C" is now less than "1A".  Below are examples.

A numeric chromosome list is ordered as below:

* 5
* 6
* 7
* 8
* 9
* 10

Chromosome names starting with numeric values are ordered as below:

* 10A
* 10B
* 1A
* 1B
* 1C
* 2A
* 2D
* 9A
* 9B

Chromosomes with mixed numeric/non-numeric will be ordered as below:

* 1
* 11C
* 1A
* 2D
* MYCHROM
* TRIPSACUM
* WheatChrom

Calling this plugin with numeric values for the chromosome from the command line would look like this (NOTE: not all parameters given below are required):

```
#!java

./run_pipeline.pl -fork1 -DiscoverySNPCallerPluginV2
-db /Users/lcj34/git/tassel-5-test/tempDir/GBS/Chr9_10-20000000/GBSv2.db
-sC 9 -eC 10 -mnLCov 0.1 -mnMAF 0.01 -deleteOldData true -endPlugin -runfork1
```

Calling this plugin with non-numeric values for the chromosome from the command line would look like this (NOTE: not all parameters given below are required):

```
#!java

./run_pipeline.pl -fork1 -DiscoverySNPCallerPluginV2
-db /Users/lcj34/git/tassel-5-test/tempDir/GBS/Chr9_10-20000000/GBSv2.db
-sC "1A" -eC "TRIPSACUM" -mnLCov 0.1 -mnMAF 0.01 -deleteOldData true -endPlugin -runfork1
```
To call DiscoverySNPCallerPluginV2 from program code : This example includes a reference genome, which is optional:

```
#!java

	         new DiscoverySNPCallerPluginV2()
                .inputDB(GBSConstants.GBS_GBS2DB_FILE)
                .minMinorAlleleFreq(0.01)
                .minLocusCoverage(0.1)
                .includeGaps(false)
                .includeRareAlleles(false)
                .maxTagsPerCutSite(64)
                .startChromosome(9)
                .endChromosome(10)
                .referenceGenomeFile(myRefFile)
                .performFunction(null);
```

#FAQS

**Running DiscoverySNPCallerPluginV2 hangs without an error message being printed - what is the problem?**

You might have an issue with the SQLite jar.  TASSEL is currently (as of October, 2022) compiled with an sqlite jar that supports the Apple M1 chip.  This should work for all machines (including those without an M1 chip).  However, we've found with larger dbs, the only sqlite jar that works is the 3.8.5-pre1 jar.  YOu can find that jar at https://repo1.maven.org/maven2/org/xerial/sqlite-jdbc/3.8.5-pre1

Replace the sqlite jar found at tassel-5-standalone/lib with the jar found at the maven repository above and rerun

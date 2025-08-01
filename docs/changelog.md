---
title: "Tassel 5 Change History"
---

## (V5.2.96) March 31, 2025

- Added checks to GenerateRCode.tableReportToVectors() to handle null values
- Changed TableReport.toStringTabDelim() to add string NaN instead on NULL for null values
- Display Table Report values that are null as NaN
- Fix displaying StepwiseOLSModelFitterPlugin results on GUI

## (V5.2.95) February 22, 2025

- Changing ListStats to use WeakHashMap instead of HashMap to prevent memory leak

## (V5.2.94) August 20, 2024

- New phg.jar v1.10
- Fixed calculation of numeric hybrid genotypes

## (V5.2.93) December 21, 2023

- New phg.jar v1.9
- Updated log4j v1.2.13 to v2.21.1

## (V5.2.92) November 10, 2023

- New phg.jar v1.8
- update postgresql-9.4-1201.jdbc41.jar to postgresql-42.6.0.jar

## (V5.2.91) October 24, 2023

- New phg.jar v1.7
- Ignore lines starting with # when reading bedfiles.
- Hide passwords in logging

## (V5.2.90) July 11, 2023

- New phg.jar v1.6
- Ignore lines starting with # when reading bedfiles.
- Hide passwords in logging

## (V5.2.89) May 19, 2023

- New phg.jar v1.5
- Added protobuf-java-3.23.0.jar, protobuf-java-util-3.23.0.jar, protobuf-kotlin-3.23.0.jar, and error_prone_annotations-2.19.1.jar
- removed "order by position" clause that made query too slow
- Fix GenomicSelectionPlugin progress reporting

## (V5.2.88) March 14, 2023

- New phg.jar v1.4
- Change GenerateRCode.fastassociation parameter from Double to double
- Add continuous integration script
- Adding generic way for TASSSEL to report logging about PHG or other packages

## (V5.2.87) January 10, 2023

- New phg.jar v1.3

## (V5.2.86) October 11, 2022

- New phg.jar v1.2
- Add new option to CreateHybridGenotypePlugin for numeric genotypes
- To support the new MAC M1 chip upgrading the sqlite jdbc jar from 3.8.5-pre1 to 3.39.2.1 and upgrading the snappy jar from 1.1.1.6 to 1.1.8.4

## (V5.2.85) September 22, 2022

- New phg.jar v1.1
- Added method to convert GIGWA Dataframe to GenotypeTable for rTASSEL

## (V5.2.84) August 30, 2022

- New phg.jar v1.0
- Add dependency jars commons-io-2.11.0.jar and sshj-0.32.0.jar
- Change TASSEL's PositionListBuilder.validateOrdering() to report out of order positions as INFO instead of ERROR

## (V5.2.83) August 4, 2022

- New phg.jar v0.0.40
- Moved calls to myLogger to the end of LoggingUtils methods to prevent warnings about Log4j loggers not setup
- rTASSEL - Fix SummarizedExperiment
- Updated HTSJDK from version 2.23.0 to 2.24.1

## (V5.2.82) June 17, 2022

- New phg.jar v0.0.39
- fixed chromsome type to String in DiversityAnalyses

## (V5.2.81) April 8, 2022

- New phg.jar v0.0.38
- Added MultivariateStepwisePlugin - A multi-trait multi-locus stepwise approach for conducting GWAS on correlated traits

## (V5.2.80) February 17, 2022

- New phg.jar v0.0.37
- Update to Kotlin 1.6
- For VCF imports, change -noDepth to default to false

## (V5.2.79) January 19, 2022

- Update biojava to latest version
- New phg.jar v0.0.36

## (V5.2.78) November 16, 2021

- New phg.jar v0.0.35

## (V5.2.77) October 27, 2021

- New phg.jar v0.0.34

## (V5.2.76) October 26, 2021

- New phg.jar v0.0.33

## (V5.2.75) October 4, 2021

- New phg.jar v0.0.32

## (V5.2.74) September 10, 2021

- New phg.jar v0.0.31
- Added imagej repository to pom for jhdf5 library
- modifies pom dependency to ejml-ddense
- Improved memory usage of TaxaList
- Corrected TaxaListBuilder.addAll(Taxon[]) to not rebuild taxa
- updates ejml library and DoubleMatrix classes
- Added MeanR2FromLDPlugin
- Added method GenerateRCode.asTasselDistanceMatrix() from Brandon
- Merged Add-GWASCountingPlugin into master

## (V5.2.73) June 23, 2021

- Update PHG and TASSEL 5 to Kotlin version 1.4.32
- Add method for adding Enzymes to GBS pipeline via a configuration file.
- New phg.jar v0.0.30

## (V5.2.72) April 8, 2021

- New phg.jar v0.0.29

## (V5.2.71) March 26, 2021

- New phg.jar v0.0.28

## (V5.2.70) February 10, 2021

- New phg.jar v0.0.27
- htsjdk-2.19.0.jar -> htsjdk-2.23.0.jar

## (V5.2.69) January 14, 2021

- New phg.jar v0.0.26

## (V5.2.68) January 7, 2021

- New phg.jar v0.0.25
- PHG-480 Fix HaplotypeGraphBuilderPlugin -haplotypeIds flag
- Fix bug when loading Numerical Genotypes

## (V5.2.67) November 10, 2020

- New phg.jar v0.0.24
- Added kotlin-stdlib-jdk7-1.3.50.jar and kotlin-stdlib-jdk8-1.3.50.jar to lib directory
- TAS-1338 Add method for R that filters GenotypeTable given arrays of seqName, start position, and end position

## (V5.2.66) October 29, 2020

- New phg.jar v0.0.23
- TAS-1337 Refactor GenerateRCode.genotypeTableToDosageIntArray to return byte[][] instead of int[]
- Added "tagFile" parameter to GetTagTaxaDistFromDBPlugin.

## (V5.2.65) September 10, 2020

- New phg.jar v0.0.22
- Fixed bug when loading TableReports into TASSEL

## (V5.2.64) July 9, 2020

- New phg.jar v0.0.21
- TAS-1335 changes equals method for Tuple
- TAS-1327 adds time-out tests to FastMultithreadedAssociationPlugin
- TAS-1333 Fix Manhattan Plot Tool Tips to show correct SNP ID
- Corrected GenerateRCode.tableReportToVectors() to handle âN/Aâ convert to Integer.MIN_VALUE

## (V5.2.63) June 15, 2020

- New phg.jar v0.0.20
- TAS-1332 Allowing null parameter values to work correctly.  This mirrors usageParameters()

## (V5.2.62) June 2, 2020

- New phg.jar v0.0.19
- TAS-1330 Add method to AbstractPlugin that returns Map of parameter values.
- TAS-1329 Add Plugin to return all Plugins / Parameters for given jar(s) (i.e. phg.jar)

## (V5.2.61) May 7, 2020

- This fixes issues with installing on MacOS Catalina
- New phg.jar v0.0.18
- TAS-1328 Modify TASSEL 5 Javadoc Build to Handle Kotlin Code
- Fixed problem with exporting multiple data sets.  The fix makes it possible to export different types of data sets.
- TAS-1326 Add addition parameters to GenerateRCode for fastAssociation
- TAS-1324 Add Newick Tree Support
- Modified SeqViewerPanel to automatically display NumericGenotype is available
- Added ability to CreateTreePlugin to take a DistanceMatrix as input, in addition to GenotypeTable

## (V5.2.60) February 20, 2020

- New phg.jar v0.0.17
- TAS-1321 Fix bug in FilterSiteBuilderPlugin when position list has chromosomes not in genotype data
- TAS-1318 Improve ParameterCache Handling

## (V5.2.59) December 12, 2019

- New phg.jar v0.0.16
- TAS-1317 changed method for p-value to LinearModelUtils.Ftest for fast association methods

## (V5.2.58) November 13, 2019

- New phg.jar v0.0.15
- PHG-390 Tassel updates.  Making ReadBedFile to return a 1 based closed range set which will fix the consecutive Bed region issue.
- Updating GeneratePluginCode to write out Kotlin syntax in addition to Java when auto-generating PluginParameter Methods.

## (V5.2.57) October 10, 2019

- New phg.jar v0.0.14
- TAS-1315 Add method to GenerateRCode to run GenomicSelectionPlugin from R
- implements write breakpoint file and adds to ImputeProgenyStatesPlugin

## (V5.2.56) September 17, 2019

- New phg.jar v0.0.13
- PHG-380 Updating ReadBedFile to handle non-numeric chromosomes correctly.  Also using Range.closedOpen() to follow Bed specification a bit more closely.

## (V5.2.55) September 10, 2019

- New phg.jar v0.0.12
- Updated kotlin version from 1.3.10 to 1.3.50
- TAS-1312 Fix Chromosome / GeneralPosition to prevent blocking when creating a Position
- TAS-1307 added error checks and messages for phenotype import
- Updated htsjdk library from version 2.14.0 to 2.19.0
- TAS-1307 changed error message
- TAS-1305 Add Preference, so that user can define their default locale within TASSEL

## (V5.2.54) July 25, 2019

- PHG-355 update the  positionsvcf paramet
- TAS-1306 Add TASSEL genotype import method with keepDepth and sortPositions options for execution from R
- TAS-1304 uses int[] to store history to allow more than 127 states at a position
- PHG-355 update the -positionsVCF parameter in the PathsToVCFPlugin so that it can take a position file in multiple formats
- TAS-1303 Fix Plugin Dialogs entry of decimal number for different locales (i.e. dot for grouping and comma is for decimal)
- TAS-1292 tassel 5 to r connection
- Corrections to comments in GenerateRCode
- Update comment in GenotypeToAdditiveValuesPlugin
- TAS-1292 modified createPhenotypeFromRDataFrameElements to import int as factor, data, or covariate
- Removed Boolean case from AbstractPlugin.convert() as it's not needed
- Reverting AbstractPlugin.getParameterFields to be protected.  Unneccesary to be public.
- TAS-1292 Added method GenerateRCode.exportToFlapjack() to export GenotypeTable to Flapjack format
- TAS-1301 write genotypetable export to f
- TAS-1301 Write GenotypeTable export to Flapjack

## (V5.2.53) June 27, 2019

- TAS-1299 Kmer Counting Plugin
- TAS-1291 comment main class in StepwiseAdditiveModelFitterPlugin
- adds trySplit(), uses LinearModelUtils F test instead of apache commons to get smaller p values instead of zero
- TAS-1291 reports the effects for an additive + dominance model
- Updated CreateTreePlugin button and tool tip to say âCreate Treeâ instead of Cladogram
- Replaced EqtlAssociationPlugin with FastMultithreadedAssociationPlugin in the TASSEL GUI

## (V5.2.52) March 21, 2019

- New phg.jar v0.0.8
- Update AbstractPlugin handling of exceptions to not System.exit() if no other plugins waiting on it's output
- TAS-1296 Add maxPercentNaN to RemoveNaNFromDistanceMatrixPlugin
- TAS-1295 Add utility that suppresses logging
- TAS-1294 Fixed Bug with FilterSiteBuilderPlugin when filtering by MAF -> List of Taxa -> Min Count
- TAS-1293 ConvertOldFastqToModerFormatPlugin moved from PrivateMaizeGenetics

## (V5.2.51) December 20, 2018

- New phg.jar v0.0.7
- TAS-1286 Add kotlin compile to TASSEL pom.xml
- Concept for autogenerating R code
- Updated ShowParameterCachePlugin to show plugin, parameter, and value columns
- TAS-1288 Create Plugin to add reference to GenotypeTable
- TAS-1290 modifications to fix handling of heterozygous genotypes by
- Fixed problem with converting integers and doubles in different locales.
- TAS-1290 modifications to fix handling of heterozygous genotypes by MLMplugin and WeightedMLM to treat 0/1 and 1/0 genotypes as a single class
- TAS-1289 Change MLM and WeightedMLM default compression level from Optimum to No Compression

## (V5.2.50) October 25, 2018

- New phg.jar v0.0.6
- adds correction for ln(0) in EmissionProbability
- removes recent fix to ViteriAlgorithmVariableStateNumber
- corrects distance adjustment in ViterbiAlgorithmVariableStateNumber
- Added utility to create HTML change log since last build

## (V5.2.49) October 4, 2018

- Added supporting library fastutil-8.2.2.jar
- Added supporting library kotlin-stdlib-1.2.71.jar
- Added case in AbstractPlugin to convert String to Character
- TAS-1285 Improve mechanics of -configParameters and add it's benefits to the GUI
- Added utility CheckSum.getChecksumForString()
- Added constants to NucleotideAlignmentConstants for homozygous diploid values of A, C, G, T

## (V5.2.48) August 2, 2018

- New phg.jar v0.0.3

## (V5.2.47) July 26, 2018

- Added plugin to pom.xml which creates sTASSEL.jar, so that itâs executable
- PHG-1284 adds test for bytes outside 0-3

## (V5.2.46) July 19, 2018

- GLM new error message thrown if users choose phenotypeOnly but does not supply a phenotype dataset
- PHG-226 Add PHG plugins to TASSEL menu
- TAS-1233 Adding distributionManagement, javadoc, source plugins to pom.xml

## (V5.2.45) July 12, 2018

- TAS-1233 Adding distributionManagement, javadoc, source plugins to pom.xml
- TAS-1233 updated , , and  in pom.xml
- TAS-1233 Added developers section to pom.xml
- TAS-1233 Added license definition to pom.xml
- TAS-1233 Added Source Code Management (scm) definition to pom.xml
- TAS-1282 Removed obsolete class FilterTaxaAlignmentPlugin.  Using FilterTaxaBuilderPlugin instead
- TAS-1282 Improved GetTaxaListPlugin to work with DistanceMatrix and Phenotype in addition to GenotypeTable
- TAS-1282 Improve PluginParameter to accept 1) file with simple list of taxa (.txt) 2) TASSEL formatted taxa list (.json or .json.gz) 3) comma separated list of taxa
- Updated FILLINFindHaplotypesPlugin to give better error message when input file canât be read.
- Fixed problem with HMAtrixPlugin related to weighting (update from Josh)

## (V5.2.44) May 17, 2018

- TAS-1279 Improve error handling of PLINK import
- TAS-1277 Add support to load files (i.e. VCF) without depth from the command line
- TAS-1276 Fixed problem when running plugins with DistanceMatrix as a parameter, the logging uses toString() which is a memory problem.
- TAS-1274 Fixed Genotype Mask function hanging when whole site is Unknown genotype
- TAS-1260 rewrite vcf import export to use HTSJDK
- TAS-1260 Added support for hapoid genotype values to BuilderFromVCFUsingHTSJDK, so that .g.vcf files can be imported
- TAS-1260 Added chromosome name to error message in BuilderFromVCFUsingHTSJDK
- TAS-1260 Added known variants to position list in BuilderFromVCFUsingHTSJDK
- TAS-1260 Updating VCFUtil to take care of Terry's comments.
- TAS-1260 Updating VCFUtil to clean up some commented out code.
- Changed Position attribute subPosition to insertionPosition
- TAS-1263 Add option to FilterSiteBuilderPlugin to remove sites that contain indels
- TAS-1266 Fixes Stepwise bug
- TAS-1266 changes to StepwiseOLSModelFitter
- StepwiseOLSModelFitter got genotypes directly from the GenotypeTable underlying the GenotypePhenotype, which resulted in incorrect genotypes being assigned to observations when the taxa in GenotypeTable and Phenotype used to create the GenotypePhenotype did not match. After the fix genotypes are retrieved with the appropriate GenotypePhenotype method

## (V5.2.43) February 22, 2018

- Modified ParseGVCF to accept reference blocks (END specified) with GT (genotype) equal 0 (reference) or .(missing)
- Fixed TASSEL Pipeline (CLI) for terminating early with threads wait for input from another thread.
- TAS-1262 Removed centiMorgans attribute from Position and added sub-position.  Also fixed bug with Chromosome.compareTo() method

## (V5.2.42) February 1, 2018

- TAS-1258 Improve Performance of MaskGenotypeCallTable (i.e. GenotypeTableBuilder.getInstanceMaskIndels())
- TAS-1259 Update htsjdk from 2.11.0 to 2.14.0
- Refactored TasselPipeline (run_pipeline.pl) -filterAlign flags to use FilterSiteBuilderPlugin instead of FilterAlignmentPlugin

## (V5.2.41) January 11, 2018

- Fixed TASSEL CLI (run_pipeline.pl) to work when using -sortPositions after -plink -ped  -map
- TAS-1255 Added ability to specify global config file (-configParameters config.txt) to the command line (run_pipeline.pl).
- TAS-1254 Improve Performance of GeneralPosition.compareTo()
- TAS-1256 Correct description of NumericalGenotypePlugin
- Changes to VCF Export to fix depth allele index out of bounds issue. Now will make sure the allele index is within the bound of what the GenotypeTable will return.
- Removing NumberFormat class which was not used,  This will help speed up creation of Position Objects.
- TAS-1253 Added supported for HindIII-NlaIII
- While we are not adding more enzyme support, this one is added for consistency as the code claimed it is supported.
- fix bug in AdditiveSiteStorePlugin and add static method to deserialize a data stored in a file.
- added return type of double to NumericAttribute

## (V5.2.40) October 26, 2017

- GBSSeqToTagDBPlugin: corrected error message from "maxKmerLength" to "kmerLength"
- ProductionSNPCallerPluginV2.java edited online with Bitbucket (corrected error message from "maxKmerLength" to "kmerLength")
- Adjusted DistanceMatrix to accept floating point numbers different by 0.0000001 instead of 0.00000001, to consider them equal.
- TAS-1251 Fixes bug. Checks for ReferenceProbability in GenotypeTable if genotype is not present.
- Added method Position.of(Chromosome, int)
- PHG-36 Adds class BackwarForwardVariableStateNumber
- Corrected handling of thread pools in LineIndexHapmapGenotypeCallTable
- modifies plugin parameter description for ImputationAccuracyPlugin
- adds code to index imputed Taxa
- adds code to ImputationAccuracyPlugin to calculate accuracy for an imputed data set that contains a subset of the taxa in the original dataset
- adds method to ImputationAccuracyPlugin to allow the imputed data to contain only a subset of the sites in the original data
- minor corrections to PhaseHighCoverage and ParentPhasingPlugin
- removes an pair of empty brackets and corrects a gui label. Does not affect plugin operation
- changes avgSegmentLength in TransitionProbability to protected
- fixes bugs in ViterbiAlgorithmVariableStateNumber

## (V5.2.39) September 7, 2017

- adds ViterbiAlgorithmVariableStateNumber
- Updated supporting library htsjdk-1.138.jar to htsjdk-2.11.0.jar
- Deals with a VCF edge case of invariant sites
- Not an elegant solution or particularly good.  All this code needs to be redone.
- Methods for getting bed files using Ranges
- TAS-1245 Corrected handling of thread pools in AbstractMaskMatrix and ParseGVCF
- Method to convert a string to diploid byte genotypes
- TAS-1244 Create simple methods for creating Position
- TAS-1231 Initial implementation of GVCF parser.
- Changed Reference Probability to Numeric Genotype on the TASSEL GUI to better describe data
- TAS-1202 Updated TASSEL VCF support to import haploid GT calls.  Mainly due to us needing to read in VCFs coming in from Sentieon's GVCFTyper for haplotype collapsing.
- TAS-1239 Removing xml.jar, xmlParserAPIs.jar, xercesImpl.jar.  This is now part of Java.
- TAS-1239 Removing unused library poi-3.0.1-FINAL-20070705.jar
- TAS-1239 Removing geronimo-spec-activation-1.0.2-rc4.jar.  This is now part of Java.
- TAS-1239 Updating TASSEL build from jhdf5-12.02.3 (local repo) to jhdf5-14.12.5 (maven)
- TAS-1243 Changed file write frequency to accommodate large db data.
- Adjusted type casting in DataTreePanel, so that code compiles with Java 9
- minor changes to SubtractGenotypesPlugin
- adds code to list crossovers by parent of origin to ImputationUtils
- TAS-1239 Replaced forester_1034.jar in local maven repo with forester-1.038.jar in maven repository
- TAS-1239 Replaced colt.jar in local maven repo with colt-1.2.0.jar in maven repository
- TAS-1235 Removed dependency on batik libraries.  Added jfreesvg-3.2.jar to write SVG files instead.
- BackwardForwardAlgorithm fix
- prevents the calculation of alpha and beta from becoming too small
- TAS-1240 Add additional access method that specifies protocol
- Added MergeRenameDeleteTaxaPlugin

## (V5.2.38) July 13, 2017

- Fixed bug with Manhattan Plots where numeric chromosomes werenât in order starting at 1
- TAS-1231 Initial implementation of GVCF parser.
- adds BackwardForwardAlgorithm to analysis.imputation package
- Fixing bug introduced by last commit
- Added in export genotype to fasta code for use in PHG.
- Added .gitignore with common files that shouldn't be committed
- Further cleaned up the GVCF Genome Sequence Classes.
- Cleaned up old classes.  Only one implementation of GVCFGenomeSequence is still here.
- TAS-1215 New methods in GenomeSEquence to get genotype as byte/string
- Fixed problem with Allele Depth when filtering Genotype Table
- updated guava-19.0.jar to guava-22.0.jar
- Provide Optional support for parse genotypes from String
- Other approaches would just through an error, this provides support for using the Optional construct.
- Evaluation of simple diploid nucleotide states.
- Commit of changes to GVCFGenomeSequenceBuilder which now marks hets as Ns.
- Corrected implementation of MaskGenotypeCallTable.genotypeAsStringRange()
- adds classes to subtract homozygous parent genotypes from hybrid genotypes to yield other parent haplotypes: DifferenceGenotypeCallTable and SubtractGenotypesPlugin
- Added LoggingUtils.setupStdOutLogging()
- Corrected implementation of Tuple.compareTo()
- Added IdentityRecognitionPlugin which clusters related taxa based on genetic distance
- Add plugin that masks Indels to Unknown (Data -> Change Indels to Unknown)
- Simplified DistanceMatrixPlugin
- Fixed problem with IBSDistanceMatrix Calculations when only one processor available
- Update to allow for multithreading and to compute various stats in a thread safe manner.
- Initial speed improvements on GVCF.
- Allow for importing gzipped GVCF files.
- Fixing some errors caused by real GVCF files. And adding in support to retrieve other GVCF annotations.
- Added in some additional functionality to GVCFSequence support.  Fixed a few bugs when building consecutive windows.  Also can now filter and mask GVCFSequences.
- Removed unused method GenotypeCallTableBuilder.getHomozygousInstance().  New method is GenotypeTableBuilder.getHomozygousInstance()
- Remove unused code from AbstractMaskMatrix
- Fixed problem with Kinship Calculations when only one processor available

## (V5.2.37) April 6, 2017

- Changed sqlite-jdbc-3.16.1.jar back to sqlite-jdbc-3.8.5-pre1.jar. This is to fix performance issues.
- fixed if(NaN) check, no number is ever = Double.NaN
- Corrected pom.xml
- Initial Commit of GVCF Genome Sequence support. Can Export a FASTA file from consecutive GVCF regions and handles Indels properly.
- Removed obsolete class SiteMappingInfo

## (V5.2.36) March 23, 2017

- removes default values from the description of ResamplingGWASPlugin parameters since those are added automatically by the AbstractPlugin.
- set minimum group number for compression to 4 in CompressedMLMusingDoubleMatrix and adjusted progress updating so that the progress bar would not disappear prematurely
- Changed Manhatton Plot to not have a legend if legend items greater than 200
- Changed Manhatton Plot to use smaller fonts when legend has many items.  This allows plot to be visable with many legend items
- Added -sortPositions flag to TASSEL Pipeline (CLI).  This sorts genotype positions during import.  Works for Hapmap, Plink, and VCF.

## (V5.2.35) March 2, 2017

- Added check to CombineGenotypeTable method to compare taxa between all combined GenotypeTables
- Improved TaxaListUtils getting union and intersected lists of taxa.
- TAS-1180 Added option to sort genotype taxa alphabetically
- Added support for Plink to SortGenotypeFilePlugin
- Removed Load and Export options from Data Menu.  New options to use are already on File Menu.
- Added support for Genotype filter with no positions. (Fixes problem with LDKNNiImputationPlugin)
- TAS-1193 Initial Plugin to Adjust Phasing in VCF based on Hapcut Output Files
- TAS-1067 Initial Abstract Filtering Translations for use with all Components of a GenotypeTable
- TAS-1195 Return code that checks for knownTags from SAM header line.
- TAS-1194 Fixed debug output so number of SNPs processed is correct.
- TAS-1024 Added support to include bed file range names (i.e. gene names) to the results of VCAPScanPlugin

## (V5.2.34) February 23, 2017

- Fixed problem with errors popping up in GUI when focus lost on empty text fields
- TAS-1191 Fixed Intersect Joining multiple Genotype Tables where Taxa Names became associated with incorrect genotypes.
- Corrected Error when creating Manhattan Plot from a previously exported GLM/MLM results
- TAS-1191 Fixed CombineGenotypeTable.genotypeRange method
- TAS-1189 : methods to store just barcode without cut site to trie
- TAS-1188 UPdate ColumnsToBinary*TablePlugins for monetdb chars
- Added support for bed files to VCAPScanPlugin
- Updated GenotypeSummaryPlugin to display reference / alternate allele in site summary
- Added method to GenotypeTable to retrieve alternate allele
- Updated BuilderFromVCF to populate alternate allele if defined.
- TAS-1024 Added stepSize parameter to VCAPScanPlugin
- TAS-1024 Added blockingWindowSize parameter to VCAPScanPlugin, Multi-threaded execution of ldak to 3 threads, added code to skip completed steps in case of a restart
- Fixed ListPlugins, as previous commit made it output nothing.
- Updated sqlite-jdbc-3.8.5-pre1.jar to sqlite-jdbc-3.16.1.jar
- TAS-1183 Fix exception thrown due to no data in snpPosition table.
- fixes bug in PhenotypeLM. The model method was attempting to cast a float[] to a double[].
- Fix the depth and taxa number on sequence processing
- ToDo to look at load factor
- TAS-1175 Updated RepGenLoadSeqToDBPlugin to handle combined reads.
- The changes should work for any reads, paried end ot not.  User should now set both minimunKmerLength and kmerLength (which is preferred maximum)
- TAS-1175 Upload changed schema
- TAS-1175 Now using 3 tables for tag alignments
- TAS-1175 Update RegGenSQLite method to pull tag correlations.
- 1175 Update RepGenLDAnalysis code
- Remove duplicate correlation calculations, write to DB more frequently to speed things up.  Adjusted indexes on database tables.
- TAS-398 Changed TASSEL Logging Dialog to collect all logging by default even when itâs not open.  File -> Preferences must be changed to send to console.
- TAS-466 Initial Implementation of FilterTaxaBuilderPlugin
- Modified Utils.getDirectory() to work correctly with forward and backward slashes.
- TAS-1181 New rAmpSeq Alignment plugin that uses BLAST output for ref tags

## (V5.2.33) January 12, 2017

- Improved performance of FilterGenotypeCallTable.genotypeForAllTaxa() which also improved performance of LD on filtered genotypes.
- TAS-1142 recent changes to imputation classes ImputeCrossProgeny, PhaseHighCoverage, and SelfedHaplotypeFinder
- TAS-1067 Initial Abstract Masking Translations for use with all Components of a GenotypeTable
- TAS-1067 Initial Abstract Filtering Translations for use with all Components of a GenotypeTable
- TAS-1178 Add variation of RepGenAlignerPlugin
- New algorithm for creating ref tags
- Added User Manual URL to CreateTreePlugin
- TAS-466 Improved performance of FilterTaxaBuilderPlugin
- TAS-466 Improved performance of FilterSiteBuilderPlugin
- Improved FileLoadPlugin to only show sort check box active when file types guess, vcf, hapmap, and plink chosen.
- TAS-466 Changed "Site Names List" in FilterSiteBuilderPlugin to use new PluginParameter type
- TAS-1142 additions to analysis.imputation package classes ImputationUtils, ParentPhasingPlugin, and SelfedHaplotypeFinder
- Added support to PluginParameter to be type .  It converts comma separated list of string (no spaces) to List
- TAS-944 Updated ConvertAlignmentCoordinatesPlugin to reorder depth and other site scores in addition to genotypes and position list
- Improved Genetic Distance heat map to sort taxa based genetic distances
- TAS-466 With FilterSiteBuilderPlugin, made it easier to filter whole chromosomes.  Now simply don't specify start and end positions, and it will keep whole chromosomes
- Improved ExportPlugin to eliminate genotype formats when GenotypeTable has no genotypes
- Deprecated FilterSiteNamePlugin
- TAS-1051 Add ability to mask by depth to MaskGenotypePlugin
- Added BitSet Collector
- Corrected check for Float.NaN in LinkageDisequilibrium when accumulating r-squared values
- Added convenience methods to ListSiteStats
- TAS-1175 ADded method to retrieve correlations for tags.

## (V5.2.32) December 15, 2016

- Added support to FileLoadPlugin to sort positions of Hapmap file during import.
- TAS-1175 Changed tagCorrelations table entry names.
- TAS-1175 Made change so code doesn't run analysis of tag against itself
- TAS-1177 Fixed SortGenotypeFilePlugin to sort depth in addition to genotypes and positions for VCF Files
- TAS-1175 Added r-squared calculation.
- TAS-1175 Add missing class TagCorrelationInfo.java
- TAS-1175 Initial RepGenLDAnalysisPlugin and related code.
- TAS-1176 Add additional error message when SNPCutPosTagVerification fails
- TAS-1152 Allowed for pulling a subset of the RefTag alignments.
- TAS-1091 Modified error message duplicate taxa test in AbstractFixedEffectLM
- TAS-1152 Added method to retrieve alignments for non-ref tags from the db
- TAS-466 Initial Implementation of FilterTaxaBuilderPlugin
- TAS-466 Faster implementation of FilterSiteBuilderPlugin
- Improved Manhattan Plot Position Labels and Plot Tool Tips.
- TAS-1152 Added aligning to ref-complement of refTag
- TAS-875 Improved performance of LineIndexHapmapGenotypeCallTable.genotype()
- TAS-1152 Added separate parameter to specify cluster hitCount Previously the minCount parameter was used both to determine which DB tags to use for creating kmer alignment seeds, and to determine number of hits required in a cluster to for a ref tag to be created.  For greater flexibility this has been split into 2 parameters:  minTagCount and minHitCount
- Removed obsolete class src/net/maizegenetics/analysis/imputation/QualityChecksPlugin.java
- TAS-1152 Turned off chrom9-only processing - left over from testrun
- TAS-1152 FIxed null pointer exceptions with tagAlignment loading.
- TAS-1152 more RepGenAlignerPlugin updates. Not quite ready for prime-time - getting a null-pointer exception in putTagAlignments
- TAS-1152 Updates to the repgen schema, table access methods and aligner code
- catches error when calculation of pvalue throws an exception in StepwiseAdditiveModelFitter
- Updated DialogUtils to print logging message instead of opening dialog when running from command line
- TAS-1152 ADded more modularity to aligner.
- TAS-1152 Updated clustering and maxima algorithm, added neobio calls
- TAS-1169 Adding in support for Covariate terms in Phenotype preprocessing Plugins.  Also added in some additional descriptions for plugin parameters.
- TAS-1173 Adding neobio methods from sourceforge for use in RepGen
- TAS-830 Add Export TaxaList as table
- TAS-814 Add Average Minor Allele Freq to Genotype Summary
- Update TASSEL About dialog
- Moved Load and Export to File Menu and renamed to Open, Open As..., and Save As...
- TAS-1167 Fixed H matrix plugin failing symmetry test due to rounding error
- TAS-1166 New plugin to create GOBII IFL files from Terry's gwas data.
- TAS-1018 Added initial implementation of GOBII Avro Format
- TAS-1004 Added User Manual URL to CreateHybridGenotypesPlugin
- TAS-1152 Cleanup of finding clustered hit counts
- TAS-1152 Updated RepGenAlignerPlugin, added moving-Sum and finding maxima
- Refactored FileLoadPlugin to follow PluginParameter design. And added Load and Load As... options to the Data menu. Load immediately browses to files to open with Guess option.
- TAS-1155 Updating AveragePhenotypeByTaxa and SubtractPhenotypeByTaxa to allow for more flexibility
- TAS-1150 Updating tag table in schema to include isReference
- TAS-1162 Fixed problem when creating IUPAC from VCF files
- Added convenience method runPlugin() to FileLoadPlugin
- TAS-1159 and TAS-1160 Initial implementations of the AverageByTaxa and SubtractByTaxa Phenotype Transformation plugins.
- TAS-1152 FIxed SQLite schema issues.
- TAS-1152 first pass at RepGenAlignerPlugin, some db changes
- TAS-1136 Updating B4RPlugin to handle empty strings for WHERE clause building. Allows for use in CLI without specifying all parameters.
- TAS-1151 Updated RegGenLoadSeqToDBPlugin using RepGen schema
- TAS-1150 Adding initial RepGen schema and SQLite access methods
- TAS-1136 Initial implementation of B4R Phenotype Plugin.
- TAS-1018 Add Avro 1.8.1 supporting library
- Improved reading of Tab-delimited files represent floating point columns as doubles.  This allows sorting of those columns to work correctly
- Added export option to write DARwin (.dis and .don) from a Distance Matrix.
- Changed FileLoadPlugin guess function to default to Tab-delimited Table if no other format can be determined.
- Add flag to FileLoadPlugin to optionally sort positions when importing.
- Refactored FileLoadPlugin to follow PluginParameter design.  And added Load and Load As... options to the Data menu.  Load immediately browses to files to open with Guess option.
- TAS-1018 Initial Implementation of GOBII Plugin
- TAS-1142 fix bugs in ImputeProgenyStatesPlugin and ImputeCrossProgeny
- TAS-1141 fixes bugs in PhaseHighCoverage and SelfedHaplotypeFinder
- adds ParentPhasingPlugin, ImputeProgenyStatesPlugin and a few supporting classes

## (V5.2.31) October 20, 2016

- Added TASSEL Pipeline parameter -exportIncludeDepth
- TAS-1085 Add Support for Maven to TASSEL
- TAS-913 Added export option for LDAK (PLINK) formatted phenotype file
- Refactored ExportPlugin to follow new Self-Describing Plugin Design.
- Fixed GenotypeTable filtering of site range where range includes whole table.  No new FilterGenotypeTable created.
- Added User Manual Link to RemoveNaNFromDistanceMatrixPlugin
- TAS-1139:  new plugin that pulls tag/taxadist from the GBSV2 DB
- TAS-1138 Tweak to storing reference tags in db when created on the fly
- Only add ref tag data if the ref tag was already in the DB.
- Update User Manual URL for VCAPScanPlugin
- Add User Manual URL to GenotypeSummaryPlugin
- TAS-336 Improved layout of Help tab for Plugins
- TAS-1024 Added icon for VCAPScanPlugin
- Reorganized the Analysis Menu.
- Implemented FilterGenotypeCallTable.genotypeForAllSites() for more efficient handling when no sites have been filtered.
- TAS-1130 Fixed Export fails for filtered HDF5 with depth
- Add option to MultiDimensionalScalingPlugin that removes NaNs from matrix before MDS performed.
- TAS-336 Added button to Plugins dialog that open URL to User Manual.  If a Plugin overrides pluginUserManualURL(), the button will open specific web page for that plugin.
- TAS-1127 Toss tags when adjusted start or end position is negative.
- This was coded correctly for forward strands.  It needed adjustment for sequences that were reverse complemented by the aligners.
- TAS-1120 Clear tags when clearing tagTaxaDist table
- This method is only called from GBSSeqToTagDBPlugin when a user wants a repeat run without deleting data.  Tags and the taxa distribution are both added back  when the plugin finishes.
- TAS-1120 Fixed creation of tag/taxadist map when pulling data from the db.
- TAS-1098 Updates to HapBreakpoints_IFLFIlePlugin
- THis plugin creates intermediate files for loading into GOBII
- TAS-1098 initial plugin to load hap breakpoints
- Also updates to other GOBII IFL plugins
- TAS-1066 GenomeSequenceBuilder with start value of 0
- GenomeSequenceBuilder.chromosomeSequence() expects 1-based input.  If the "start" parameter is less than 1, it will return null.  Otherwise, the user/calling-method must get this right.
- TAS-558 Move removeSitesâ¦() methods from GenotypeTableUtils to GenotypeTableBuilder
- TAS-612 Migrate net.maizegenetics.dna.snp.depth to net.maizegenetics.dna.snp.score

## (V5.2.30) September 1, 2016

- Updated supporting library guava-14.0.1.jar to guava-19.0.jar
- TAS-1105 More updates to GOBII plugins to handle GOBII IFL changes.
- Added order by clauses to MonetDB_IFLFilePLugin when creating marker and dnarun id lists.
- TAS-1105 updated BL IFL plugins to match curenct GOBII default values.
- Also fixed dnarun mapping to NOT use "num" as comparing with this to get dnasample name returns wrong results from postgres when num col is null This required changes to dnarun.dupmap as well.
- TAS-1018 Initial Implementation of GOBII Plugin
- Fixed bug when commas used instead of decimal points (i.e. some locale) during export of Table Reports (i.e. LD Results)
- TAS_1105 added "order by" to sql statement that pulls dnarun and marker ids.  THis ensures proper ordering.
- TAS-1018 Initial Implementation of GOBII Connection Utility to Postgres database.
- TAS-1105: final updates to class that creates marker_idx and dnarun_idx
- This class creates intermediate files for loading when the marker_idx and dnarun_idx values were not previously populated in the dataset_marker and dataset_dnarun tables.  It should become obsolete as these fields are now populated with the rest of the marker/dnarun values.
- TAS=1105 Updated scripts to handle dnarun_idx properly.
- TAS-1105 Added script to update marker_idx and dnarun_idx per Kevin's instructions
- TAS-1074 updated Marker and germplasm plugins.
- Added ".trim() when processing mapping file columns.  Changed markerdnarun plugin to use "num" for dnasample rather than "platename". The dnarun.nmap file has been altered to use neither as the GOBII IFL scripts do not handle comparing to postgres IS NULL.
- TAS-1107 removes .inDir() from donorFile PluginParameter so that isOutputProjection option will work.
- TAS-1018 Initial Implementation of GOBII Connection Utility to BMS MySQL database.

## (V5.2.29) August 4, 2016

- TAS-1085 Add postgresql-9.4-1201.jdbc41 to Maven Pom
- TAS-1105 Plugins moved from PMG to TASSEL for creating GOBII intermdediate files for IFL scripts
- TAS-1018 Initial Implementation of GOBII Connection Utility to Postgres database.
- TAS-1104 Added code to remove kmers not meeting minimum length requirement
- TAS-1103 Fixed FilterGenotypeCallTable incorrectly fills in missing genotypes for unknown taxa
- TAS-1101 Add Ability to Export Depth from GenotypeTable
- TAS-432 Initial display of Depth in Tassel GUI.
- TAS-1018 Change logging of TASSEL Pipeline arguments to print ????? after any flag containing "password" instead of the real value.
- TAS-1018 Initial Implementation of GOBIIGenotypeCallTable
- TAS-1023 Fixed creation of Position object to include strand information.
- TAS-1092 Updated comments
- TAS-1092 updated comments for ColumnsToBinary plugins.
- TAS-1092 added support for allele column, cleaned up comments and debug
- TAS-1092 Updated byte processing, made missToZero and negToZero apply to all numeric types
- TAS-1094 Added additional user messages
- TAS-1092 Added plugins to create monetdb binaries.
- adds ModelEffectUtils change needed for FixedEffectLM changes
- adds code to FixedEffectLM, AbstractFixedEffectLM, DiscreteSitesFELM, and FixedEffectLMPlugin to display additive and dominance effect estimates in stats output

## (V5.2.28) July 7, 2016

- TAS-1075 Further adjustments to likelihoodRatioThreshAlleleCnt.
- Changed maxCountAtGeno to 100 as the faulty algorithm needed over 20% of the depth to be minor allele when depth was 200.  This algorithm still needs syncing with VCF once that algorithm is straightened out.
- TAS-1087 Drop alignments that go off the chromosome end.
- Bowtie-2 in local mode uses soft clipping.  Our adjustment algorithm called from SAMToGBSdbPlugin adds back the clipping, changing the sequence start position relative to the reference.  If this start is off the chromosome, drop it.
- Added AlleleDepth.depthForAllele(taxon, site, allele) convenience method.
- modified AbstractFixedEffectLM to fix column head error
- deletes AbstractOLSModel and DiscreteSiteOlsSolver because they are no longer used by anything else
- adds new methods to Haplotype, HaplotypeCluster, and BiparentalHaplotypeFinder
- TAS-1018 Initial Implementation of GOBII Connection Utility to Postgres database.

## (V5.2.27) June 23, 2016

- Fixed bug when loading VCF files. (i.e. number of values: 4096 doesn't equal number of sites: 647941)
- TAS-1075 More adjustments to handle the binomial dist calculation error in BasicGenotypeMergeRule
- TAS-806 Changed outputFormat to NOT required, default to ABH format.
- Fixed bug when converting from HDF5 to VCF.
- TAS-806 allow additional output formats for GenosToABHPlugin
- TAS-1085 Add Support for Maven to TASSEL
- TAS-1018 Initial Implementation of GOBII Connection Utility to Postgres database.

## (V5.2.26) June 16, 2016

- Implemented GenotypeCallTable.genotypeAsStringArray(int, byte).
- TAS-1078 Fix null pointer exception when ref tag returns null due to long string of T's.
- Fix -export command line flag when specifying multiple outputs to a directory.
- TAS-612 Migrate net.maizegenetics.dna.snp.depth to net.maizegenetics.dna.snp.score
- TAS-1076 VCF Updates to fix a number of bugs.
- TAS-1075 Fixed math error when comparing int to float.
- Implemented TaxaListBuilder.contains() method.
- Implemented FilterGenotypeCallTable.allelesSortedByFrequency() method
- Implemented CombineGenotypeCallTable.diploidAsString() method.
- Minor code corrections.
- TAS-794 Optimization to Endelman Kinship Matrix calculation for improved progress reporting.
- TAS-1067 Initial Abstract Masking for use with all Components of a GenotypeTable
- Update comments in TasselPrefs
- Added method GenotypeCallTable.isSiteOptimized()
- TAS-944 Updates to Genome version conversion plugin
- Added new method BitSet.getAndClear()
- TAS-1067 Initial Abstract Filtering Translations for use with all Components of a GenotypeTable
- TAS-336 Added PluginParameter instance that's just a label.
- TAS-1004 New plugin to create hybrid genotypes.
- TAS-217 Corrected Filter Genotype Table Sites to exclude Site Names when "Include Sites" not checked.
- TAS-1063 Initial Implementation of Detect Inversion Plugin
- TAS-1066 String support for GenomeSequence
- adds method to AlleleProbabilityFELM to conform to FixedEffectLM interface
- TAS-1064 changes to FixedEffectLM and its sub classes to implement site filtering and to add a minor observation column count to the site report
- TAS-1060 Fixing bug with new GenomeSequence builder instance
- TAS-1060 new GenomeSequenceBuilder instance useful for testing
- Concept for rapid block based imputation of KNNi
- TAS-884 Fixed Diversity Plugin
- TAS-1024 Initial Generalize VCAP Scan Plugin
- TAS-466 Added support to FilterSiteBuilderPlugin to filter genotypes of a GenotypePhenotype data set.

## (V5.2.25) April 28, 2016

- Adjustments to AbstractPlugin to output error logging message at bottom.
- Better error handling for GetPositionListPlugin and GetTaxaListPlugin
- TAS-988 Function to Add one or more Normalized_IBS matrices
- TAS-1041 Changed button name for GBSv2 GBSSeqToTagDBPlugin for the GUI
- TAS-1041 GBSv2 Pipeline on GUI and deprecate GBSv1
- Changed menu choice Help -> Help Manual to open Tassel Wiki in browser.
- TAS-466 Updated FilterSiteBuilderPlugin to show user messages if resulting dataset is same as original or empty.
- Improved FSFHapImputation error messages.  Suggestions from Jason
- Add better error message when loading files and the guess function can't determine file format.
- TAS-579 Added logging to show number of positions or taxa imported from JSON formatted position lists or taxa lists.
- Updated MergeAlignmentSameSitesPlugin to handle annotations in the hapmap files.
- Added -printMemoryUsage flag to command line.  This can be used at any point in the pipeline to print current amount of memory used.
- TAS-1024 Initial Generalize VCAP Scan Plugin
- TAS-817 Corrections to Comments in GCTADistanceMatrix.
- TAS-1029 Toss tags in Discovery when the tag created from the reference genome contains ambiguous alleles.
- Changed ignoreDepth to keepDepth as we were really keeping depth when true
- VCF fixes for PL and GQ Field and Progress Bar
- TAS-779 Added logging to ImputationAccuracyPlugin to show which genotype was used for original, masked, and imputed.
- TAS-1019 Fixed FILLIN when NULLPointerException when output file doesn't specify absolute or relative path.

## (V5.2.24) April 7, 2016

- TAS-1024 Initial Generalize VCAP Scan Plugin
- TAS-779 Initial implementation of MaskGenotypePlugin
- moved TeosinteHaplotypeFinder to privatemaizegenetics
- TAS-1021 Fixed Filter by Position ends up with 1 position on next chromosome
- adds code that tests for replicated taxa. If taxa are replicated GLM throws a runtime error.
- TAS-391 :  Turned off debugging
- TAS-391 Change readQualityScoreFile to accept additional columns beyond the first 3 required columns.
- TAS-391 : Delete old data defaults to true; Unimplemented params in DiscoverySNPCallerV2 commented out.
- TAS-1019 Add "Browse" Button to "Donor" field of FILLIN imputation.
- TAS-1012 Reduce the default LD search for KNNi to 10Mb
- TAS-1019 Fixed FILLIN from GUI throws null error message
- TAS-371 Minor change to adjustCoordinates() to handle hard clipping as well as soft.
- TAS-391 : Using cutPosition instead of SAMPosition for reverse strand tags
- TAS-391 Updated SAMToGBSDBPlugin to store correct cut position for reverse strands.  Jeff's DiscoverySNPCallerPluginV2 changes to follow.
- TAS-391 Added parameter to SAMToGBSDBPlugin that specifies minimum mapping quality for a tag to be included in the DB.
- TAS-1010 Reorganized Filter Menu to make clear preferred filters to use and future direction.
- TAS-1017 Initial concepts for repeat GBS processing
- TAS-520 Fixing of the indel problem and code to debug the SNP caller
- TAS-1011 Default Bins on Histogram to 25 (previously 5)
- TAS-1010 Fixed Filter Genotype Table Sites - "S N P" has weird spaces

## (V5.2.23) March 31, 2016

- adds TeosinteHaplotypeFinder used for imputing Teosinte population
- TAS-1006 Add reference allele to GBSv2 data base
- Fixed VCF Loading Bug.  Due to incorrect export.
- adds AbstractOLSModel and DiscreteSiteOlsSolver classes which will be used to implement a multithreaded version of GLM
- Remove unnecessary logging during VCF export

## (V5.2.22) March 17, 2016

- TAS-1001 Fixed code that creates a reference tag when a reference genome is included with DiscoverySNPCallerV2Plugin. Forward and reverse tags for a cut position are now aligned separately.
- corrects taxa labels on residuals in CompressedMLMusingDoubleMatrix
- modifies MLMPlugin to calculate residuals from a kinship and phenotype (without a Genotype).
- TAS-466 Now allowing optional "chr" and "pos" column headers in file used with FilterSiteBuilderPlugin -chrPosFile. Also improved error messages.
- TAS-794 Optimization to Endelman Kinship Matrix calculation for improved progress reporting.
- fixes memory leak by moving ForkJoinPool to class variable in AdditiveModelForwardRegression
- TAS-1001 Initial changes to SamToGBSDBPlugin to correctly add strand informtation to DB.
- adds code to handle potential error conditions in StepwiseAdditiveModelFitter

## (V5.2.21) March 3, 2016

- Added restriction enzymes PstI-BfaI and SphI-EcoRI
- TAS-931 Added Dominance Normalized relationship matrix (Dominance_Normalized_IBS)
- TAS-993 Improve Memory Efficiency of DistanceMatrix by storing only half the matrix, and storing values as floats instead of doubles.
- fixes bug in StepwiseAdditiveModelFitter correctly matching genotypes and phenotypes
- Speed improvements to taxa summary option of GenotypeSummaryPlugin.
- TAS-966 Added option to use Proportion of Heterzyogotes for Kinship method Centered_IBS_Dominance
- TAS-891 Create Function to remove NaN from Distance Matrix
- changes getLogger to use the correct class in GenomicSelectionPlugin
- Fixes ClassicMds bug causing it to fail when using the EJML matrix library.
- Deprecated to whole class net.maizegenetics.analysis.distance.Kinship.
- TAS-944 : Now retains depth and GenotypeTable annotations
- TAS-925 Added Distance Matrix Annotations to Tassel GUI Information Window.
- Added DistanceMatrixUtils.getIBSDistance() that was previously there.
- TAS-988 Function to Add one or more Normalized_IBS matrices
- TAS-946 Added convenience method SubtractDistanceMatrixPlugin.wholeMatrix(DistanceMatrix) for setting whole matrix programmatically.
- TAS-987 Function to Add one or more Centered_IBS matrices
- Changes to FILLINDonorGenotypeUtils that affect changes in TAS-985 and TAS-986
- Changes to FILLINImputationPlugin that adds support for current directory "short" pathnames for inputs and creates the necessary directories for the specified outpath (TAS-985) and fixes a bug for projection where the high-density file donating haplotypes could only be hmp or h5 (TAS-986)
- TAS-378 Corrected implementation of CombineGenotypeTable.genotypeAllTaxa() which was causing problems when exporting union joined genotype tables.
- TAS-946 Added better error message when whole matrix used by SubtractDistanceMatrixPlugin doesn't have necessary annotations.
- Increases size of siteQueue and attempts to capture a couple of error types in FastMultithreadedAssociationPlugin
- fixes bug in GenotypePhenotype that affected TableReport display. Changes the method indexOfGenotype.

## (V5.2.20) February 4, 2016

- TAS-800 Implemented FilterGenotypeCallTable.diploidAsString() method
- TAS-800 Increase speed of Hapmap Export
- TAS-946 Function to Subtract one or more Normalized_IBS sub-matrices from larger matrix
- TAS-960 Added Support to Import Kinship Matrix from MultiBLUP format
- TAS-463 Updated javadoc for GenotypeTable.genotypeAsString()
- TAS-797 Improved Support to Export Kinship Matrix to MultiBLUP format
- TAS-946 Added subclass DistanceMatrixWithCounts for storage of counts when available for the distance matrix.

## (V5.2.19) January 14, 2016

- TAS-955 Fix for ArrayOutOfBoundsIssue caused by not checking the alleles in the GenotypeTable when Variants are around.
- TAS-954 HDF5AlleleDepth constructor needs to create TaxaArrayList using buildFromHDF5Genotypes() to ensure only taxa with depths are added to the list.
- TAS-924 Added error message to EndelmanDistanceMatrix.subtractEndelmanDistance() when sub-matrix created with older build of Tassel.
- Updated getCitation()
- Fixed missing author in getCitation()
- Fixed bug in runPlugin() and updated getCitation()
- Adds plugin that merges breakpoint files
- TAS-938 Added method to convert seq to int with padding at beginning for monetdb.
- Initial commit of VCF fixes to better handle indels.
- corrects taxa count for breakpoint file in WritePopulationAlignmentPlugin
- Adds method to FSFHap imputation that uses parent sequence as haplotypes rather than inferring them from the progeny. The method was added to handle EU-NAM 50K array data. Adds the calss UseParentHaplotypes and modifies CallParentAllelesPlugin, PopulationData. Also adds a method to write breakpoint files to WritePopulationAlignmentPlugin.
- updates to SolveByOrthogonalizing to get data and svd results. updates to FastMultithreadedAssociationPlugin to improve performance.
- finishes implementation of FastMultithreadedAssociationPlugin. Implements additive only models. Dominance may be added later.
- TAS-898 TAS-912 Adjustments for 0-based bed files.
- TAS-925 Removed IBSDistanceMatrix as a datatype.  Same data is now stored as a DistanceMatrix with annotations.
- TAS-925 Added code to guess file import type of newly annotated Distance Matrix files.
- Additions to FastMultithreadedAssociationPlugin
- Adds FastMultithreadedAssociationPlugin
- Update EqtlAssociationPlugin and SolveByOrthogonalizingPlugin make some efficiency improvements
- Add new constructors for NumericAttribute that do not require a BitSet for missing.

## (V5.2.18) December 10, 2015

- TAS-936 Make likelyReadendStrings have "package" visibility for junit tests.
- TAS-936 Fixed identification of read-end cut site
- The index of the likely read end cut site was not identified when it was not the first option in the array of likely read ends.  This resulted in both the cut site remmants and all other data beyond remaining as part of the tag, thus preventing some tags from aligning properly with the exernal aligners.
- Removed DistanceMatrix.setDistance() method.  DistanceMatrixBuilder should be used now.
- TAS-935 Adding nucleotideBytetoString() method TOo many files have a copy of nucleotideBytetoString() method - it should live in the NucleotideAlignmentConstants.java file.
- TAS-932  Added error message on fastQ file name when culling files not represented in the key file.
- TAS-924 Function to Subtract one or more Centered_IBS sub-matrices from larger matrix
- TAS-925 Add GeneralAnnotation to DistanceMatrix (including import / export to file)
- TAS-858 Changed KinshipPlugin method names.  Scaled_IBS -> Centered_IBS; GCTA -> Normalized_IBS; Dominance -> Dominance_Centered_IBS; Dominance_Normalized_IBS (Not yet implemented)
- Updated PreferencesDialog to be a Self-Describing Plugin.
- Changed the TableReport viewer to show decimal numbers with full precision.
- TAS-898 Changed -bedFile option of FilterSiteBuilderPlugin to make end position exclusive.  Start position already inclusive per BED file specification.
- Improved use of Tassel Pipeline.  If defined fork doesn't receive input from another plugin, then it will be started regardless if -runfork was specified.  Therefore, -runfork flags are no longer needed.  But they are still valid flags as before.
- fixes progress bar bug caused by small data sets in RandomGenotypeImputationPlugin

## (V5.2.17) November 12, 2015

- TAS-921 modifies RandomGenotypeImputationPlugin to optionally write output directly to a hapmap file rather than an in-memory GenotypeTable
- TAS-858 Added Dominance relationship matrix calculation to KinshipPlugin
- TAS-918 Create Global flag that can limit number of threads a function will use.
- Minor improvements to exporting DistanceMatrix as text square matrix.
- TAS-876 Adds scroll bars (JScrollPanes) to selection lists in the dialog for QQPlot in QQDisplayPlugin
- TAS-908 Adds method to DistanceMatrix to calculate the Hadamard product of two DistanceMatrix's. Implements the method in KinshipPlugin. Adds equals() method to TaxaArrayList because it is needed to check input matrices prior to computing the Hadamard product.
- TAS-888 Added in compression for Square matrix export
- TAS-915 Changed the PositionHDF5List constructor to handle a .h5 file with the last loci having just 1 mapped position.  Current logic skipped this scenario, resulting in excpetions thrown on .h5 file load into tassel.
- TAS-916 SNPQualityProfilerPlugin output presented incorrect data - it was off by 1 position.  Now fixed.
- TAS-912 Create plugin that creates a BED file given list of features
- changes EqtlAssociationPlugin to use report headers from AssociationConstants
- Improved information window for Genome viewer to list chromosomes, start, and end positions.
- Fixed bug with Tassel Pipeline when machine has only one core.  The pipeline is using half the cores, but with one that results in zero.
- Changed the TableReport viewer to show decimal numbers to 5 decimal places to make it easier to see value.
- TAS-900 Submitting for Karl: FindMatchByWordHash - fixed bug in recording matches with less than or equal to maxWordCopies. RNADeMultiplexProductionPlugin: Fixed counters of good reads with matches in conti/gene DB.  Also cleaned up references to ProdctionSNPCallerPluginV2.
- Fixed bug with AbstractGenotypeCallTable.minorAlleles() implementation.  Problem occurred when site had only unknown values.
- TAS-898 Add support for subsetting PositionList with BED file
- TAS-897 Add ability to filter given list of chr/positions
- TAS-902 - minor changes to RNADeMultipLexSeqtoDBPlugin when deleting exisint db.
- TAS-903 Add HMatrixPlugin that creates Hybrid Genomic Matrix from Pedigree Relationship Matrix and Genotype Relationship (Kinship) Matrix
- TAS-891 Create Function to remove NaN from Distance Matrix
- TAS-466 Initial implementation of flexible filter plugin
- Minor updates to CombineDataSetsPlugin
- TAS-900 moved debug statement to correct place.
- TAS-900 Added check for null tag created in LoadRNAContigsToGBSDBPlugin:processData.  Karl was seeing a nullPointerException from TagDataSQLite:putAllNamesTag() due to null tags in the map.  These were created when "N" was in the sequence.  Code in in LoadRNAContigs now checks for non-null tag before adding to the contigNameMap.  We should consider adding a count of null tags to be output at the end of processing.
- TAS-900 Bug in isPresent of Match fixed
- changes the older style format import in PhenotypeBuilder to set only -99 and -999 to Nan instead of any number starting with -99.
- TAS-900 Added code to handle tissue/taxon/tag not found in DB for the new methods getAllCountsForTagTissue()/getAllCountsForTaxonTissue()/getAllCountsForTissue() in TagDataSQLite
- Initial implementation of MaskGenotypeCallTable
- adds ClusterGenotypesPlugin
- changes haplotype clustering in BiparentalHaplotypeFinder
- TAS-900 Added method to return set of Tissue from DB.
- updates to ClusterGenotypesPlugin
- Updated run_pipeline.bat and start_tassel.bat to accept -Xmx and -Xms flags from the command line.  Thank you Josh for this contribution.
- TAS-875 Initial Implementation of LineIndexHapmapGenotypeCallTable
- TAS-900 Extraction of tissue based data from the DB
- TAS-900 Support for name field in Tags.
- TAS-900 Providing function call for tuple so it can be used more easily in lambdas
- TAS-903 Add HMatrixPlugin and AMatrixPlugin to Tassel GUI Menu
- TAS-903 Added PluginParamater.distanceMatrix() type that allows parameter on GUI to select between multiple Distance Matrices if selected.  And import from file when not selected.
- TAS-903 Added DoubleMatrixFactory.make(DistanceMatrix) method for more efficient conversion from DistanceMatrix
- TAS-903 Add AMatrixPlugin that creates AMatrix from PLINK Pedigree
- TAS-902 Changed Discovery so find taxon from keyfile instead of requiring fastQ file names to be encoded with underscores.
- TAS-904 : Added JavaDoc, added a convenience runPlugin method that takes a GenotypeTable object as input and outputs one as well, added a log4j logger, reported on % total and % non-missing SNPs set to missing,
- TAS-900 Refactoring to make builders and clarify naming
- TAS-902 Updated RNA discovery pipeline to calculate taxa distribution as before.  Through we don't store it in DB, we use it to determine whether tags are represented a minimum amount of time.  If not, they are deleted from our Set and not stored in the DB.
- TAS-906 fixed an import error
- TAS-906 commit Discovery for RNASeq
- TAS-900 Renaming plugin
- TAS-900 Moving RNA packages to TASSEL source
- TAS-900 Support for RNA processing pipeline
- TAS-904 : Reimplemented using a MaskGenotypeTableBuilder to save RAM
- adds code to BiparentalHaplotypeFinder, CallParentAllelesPlugin, and FSFHapImputation to add a cluster report for BiparentalHaplotypeFinder
- TAS-904 : Changed the code to use a GenotypeCallTableBuilder instead of a GenotypeTableBuilder, and to put the taxa loop inside the site loop.
- TAS-906 Stub for writing taxa/tissue data to SQL tables.  Includes TaxaTissueDist class
- Added convenience method DataSet.getDataSet(Object) for when only one data item.
- TAS-906 Fix readTissueAnnotationFile - earlier change caused an exception
- TAS-906 Added missing file from last submission.  GBSUtils should have been delivered with initializeBarcodeTrie changes.
- TAS-906 Initial stuff for laoding tissue/taxon/tag values for RNA seqs
- TAS-904 : Make the SetLowDepthGenosToMissingPlugin available on the GUI
- TAS-904 : Function plugin (command line or GUI)
- TAS-906 Added "tissueName" to barcode.  Changed GBSUtils.initializeBarcodeTrie to call BarCode constructor based on existing of tissue annotation attached to taxon.
- TAS-906 Fixed readTissueAnnotationFile() to make more efficient.
- TAS-906 Change keyfile constants to global for general use.
- TAS-906 Define tissueNameField in GBSUtils, and add method to TaxaIOUtils to return a set of tissue values from a keyfile.
- TAS-900 Clarification on tables names
- TAS-904 : Day 2 Hackathon progress
- TAS-902 created readDeMultiPlexFastQBlock and moved to GBSUtils as both RNA discovery and production need this functionality.
- TAS-900 Addition of RNA support for GBS DB schema
- TAS-904 : Initial empty class for plugin
- TAS-900 Implementation of getTagsNameMap
- TAS-900 GBS database changes to support named tags
- Allowing SQL queries to be set from the DataSet
- Make two plugin more fluent in their interfaces.  More similar to modern plugin spec.
- TAS-895 - made parallel change to GBSSeqToTagDBPlugin to speed up removal of 2nd cut site (same change as in Production)
- TAS-466 Fixed FilterSiteBuilderPlugin to accept more than 3 columns in Bed file.

## (V5.2.16) October 15, 2015

- TAS-466 Initial implementation of flexible filter plugin
- Corrected bug were distance matrix algorithms (GCTA, Endelman, IBS) hangs when number of sites is lower than the number of available processors.
- Added -printGenoSummary option to the Pipeline that prints number of taxa, number of sites, and chromosome summary of a Genotype Table.
- TAS-855 Filter GenotypeTable using PositionList
- TAS-895 Changed code to use and indexOf method to find the second cut site instead of Aho-Corasick Trie.  This runs much faster.
- adds maxThreads to ResamplingGWASPlugin
- increase FastqChunk limit to 1M
- updates to AdditiveModelForwardRegression and ResamplingGWASPlugin to allow alternate method to control max threads
- TAS-891 Create Function to remove NaN from Distance Matrix
- minor changes, mostly adding debug statements, to AbstractForwardRegression, AdditiveModelForwardRegression and ResamplingGWASPlugin in the analysis.modelfitter package
- adds static method to serialize additive sites to GenotypeAdditiveSite
- Adds check to fireProgress in EndelmanDistanceMatrix to make sure that percent

## (V5.2.15) September 24, 2015

- TAS-790 Optimization to IBS Distance Matrix calculation to only split workload into number of available virtual processors.
- TAS-794 Optimization to Endelman Kinship Matrix calculation to only split workload into number of available virtual processors.
- TAS-817 Optimization to GCTA Kinship Matrix calculation to only split workload into number of available virtual processors.
- TAS-875 Added LIXPlugin to create index for block zipped (bgzip) files
- TAS-886 ALtered ProductionSNPCallerPluginV2 to create VCF file as default output.  IF the output file suffix is ".h5", HDF5 file will be created instead.
- TAS-868 : Now can connect to a remote DB by specifying host in the connection config file (e.g., host=cbsuss06.tc.cornell.edu)
- TAS-868 Added DefaultPluginListener implementation of PluginListener that can be attached to any Plugin to get progress logging messages when writing custom code.
- TAS-466 Initial implementation of flexible filter plugin
- TAS-875 Initial Implementation of LineIndexHapmapGenotypeCallTable
- TAS-875 Implemented FilterGenotypeCallTable.genotypeForAllTaxa() for efficient calling of base GenotypeCallTable (i.e. LineIndexHapmapGenotypeCallTable) when only sites removed.
- TAS-865 Changes to GBSSeqToTagDBPlugin and ProductionSNPCallerPluginV2 to make the v2 pipeline tag calling results closer to the v1 pipeline. New library ahocorasick-0.2.4.jar added for quicker trie processing.
- Adjusted forester.jar in sTASSEL.jar manifest, so that running Tassel with this jar sets the Classpath correctly for forester.
- Added library htsjdk-1.138.jar (A Java API for high-throughput sequencing data (HTS) formats.)
- Added -printGenoSummary option to the Pipeline that prints number of taxa and number of sites of a Genotype Table.
- Added logging message for "Number of Processors" to basic environment logging.
- Added logging messages for start and finish of loading files.
- TAS-579 Changed importing of PositionList (json) to streaming design.  Before it created DOM in memory.  The streaming design is much faster and required much less memory.
- TAS-863 VCF Import Export fixing for Ref-alt allele switching.
- changes command line parameter for MultiDimensionalScaling Plugin from "Number of axes" to just axes.

## (V5.2.14) August 27, 2015

- TAS-849 Increasing cache size for HDF5 files to number of taxa times number available processors (but no more than 1/3 available memory) for better performance of parallel site based algorithms.
- fixes permutation testing in ReferenceProbabilityFELM
- changes report headers to use AssociationConstants where appropriate to avoid QQ-plot problems in the future
- fixes bug that prevented QQ plot from displaying results from ReferenceProbabilityFELM by changing the column headers to p and Pos.
- TAS-742 implements a random access of a hapmap file
- TAS-872 Changed Taxon class to not be serializable.  That class is not serializable.
- Removed unused plugin src/net/maizegenetics/analysis/imputation/GenotypeImputationPlugin.java
- TAS-348 Removed "Save Data Tree As..." functionality.  This doesn't work and has had many backward-compatibility problems in the past.
- TAS-832 Changes to several classes in analysis.modelfitter package to add the ability to use serialized sites to speed up the analysis.
- TAS-874 Alter barcode trie to store both barcode and initial cut site. This allows for correct identity of barcode when 2 barcodes differ by 1 added BP at the end (e.g. TCAC and TCACC).
- TAS-336 Added condition to AbstractPlugin.printParameterValues() to print number of positions in PositionList instead of the whole PositionList information.
- Add quality score function in Fastq reads
- TAS-832 adds code to set enter and exit limits in StepwiseAdditiveModelFitter from StepwiseAdditiveModelFitterPlugin
- TAS-466 Initial implementation of flexible filter plugin
- TAS-870 Usability Improvements to Cladogram Function
- TAS-790 Improved performance of IBS Distance Matrix Calculation using a parallel stream to process blocks of sites.
- convert Fastq to FastA for Blast alignment
- updates StepwiseModelFitterPlugin getters and setters, modifies StepwiseAdditiveModelFitter reports to replace "" with "--" in cells with no data.
- Corrected name of logger in SNPQualityProfilerPlugin - it erroneously said "SAMToGBSdbPlugin".
- TAS-832 Adds AdditiveSiteStorePlugin to provide a fast method to store and retrieve additive site lists that are re-used. Modifies AdditiveSite and the model fitter classes that use it to implement storage and retrieval of site lists.
- makes AdditiveSite extend Serializable, adds a comment to StepwiseAdditiveModelFitter, fixes debug timing output in AbstractForwardRegression

## (V5.2.13) July 30, 2015

- TAS-853:  New GBSv2 plugin to print tag sequences from specified database.
- TAS-867 Fixed Manhattan Plot problems where position exceeds max value of integer and where first position of chromosome more than last position of previous chromosome. These makes position relative simply numbering them 0 through N.
- Added 5 new restriction enzyme pairs
- TAS-834 Added missing concurrency to method's returned hashmap.
- TAS-790 Improved performance of IBS Distance Matrix Calculation using a parallel stream to process blocks of sites.
- TAS-867 Fixed Manhattan Plot problems where position exceeds max value of integer and where first position of chromosome more than last position of previous chromosome.  These makes position relative simply numbering them 0 through N.
- TAS-866 Changed gui-names for parameters in GBSSeqToTagDBPlugin and ProductionSNPCallerPluginV2 that have switched from "tag" to "kmer".
- TAS-866 Changed parameter "mxTagL" to "kmerLength" and "mnTagL" to "minKmerL" to better reflect the use/purpose of the parameter in both GBSSeqToTagDBPlugin and ProductionSNPCallerPluginV2.  Others references to "tag" were changed to "kmer".
- TAS-824 Make start/end chromosome non-required parameters.  If these parameters are not specified they are inferred from the sorted chromosome list obtained from the database.
- Fixed bug in GenomeSequenceBuilder - added early return when chromosome is not found when attempting to get sequence from chromosomeSequence(). Also added minor tweak to fullRefCoordinateToChromCoordinate() to use a parallel stream.
- TAS-832 adds AdditiveResidualForwardRegression class
- TAS-834 Changed fullRefCoordinateToChromCoordinate to return the overall coordinate(full reference sequence coordinate) as well as a Tuple containing the Chromosome and relative position within that chromosome.
- TAS-834 Added method to GenomeSequenceBuilder to convert whole genome index to relative chromosome/position
- TAS-832 makes changes to ResamplingGWASPlugin and implementing classes to debug and reorganize. Makes efficiency improvements to CovariatePermutationTestSpliterator and NestedCovariatePermutationTestSpliterator. Deletes ResidualForwardRegression.
- TAS-861 Move check of position list to before we create genotype table to prevent null pointer exception.
- TAS-832 Add code that attempts to parse SNP name to get chromosome and position when reading a numeric genotype
- TAS-832 adds StepwiseAdditiveModelFitterPlugin to modelfitter package. Modifies StepwiseAddititiveModelFitter to work with the plugin. Modifies UnaryOperator shuffleDouble() in BasicShuffler.
- TAS-832 moves ResamplingGWASPlugin from net.maizegenetics.analysis.association package to net.maizegenetics.analysis.modelfitter
- TAS-832 Adds StepwiseAdditiveModelFitter which provides a multi-threaded version of StepwiseOLSModelFitter. It does not include VIF test (yet). Also adds CovariatePermutationTestSpliterator and NestedCovariatePermutationTestSpliterator, which help parallelize the permutation test.
- TAS-851 - Added hashcode check to Chromosome "equals" method.
- Removed warning message when AS (alignment score) field is missing. THis is not a required field and there were too many messages spewing.
- Added in error checks for FORMAT issues with VCF.
- fixes bug in StepwiseOLSModelFitter that converted nucleotide genotypes to covariates incorrectly. It was not correctly coding heterozygotes.
- TAS-832 changes to AbstractForwardRegression, AdditiveModelForwardRegression, and ResidualForwardRegression to implement parallel stepwise regression
- TAS-832 add ForwardStepAdditiveSpliterator and three sub-classes to help parallelize stepwise regression
- TAS-832 cleanup and reformat AbstractAdditiveSite
- TAS-832 Cleanup and reformat AdditiveSite and derivatitive classes
- TAS-794 Implementation of Endelman Kinship Matrix Calculation
- TAS-832 Added new AdditiveSite interface and implementing classes. Modified AbstractForwardRegression, ResidualForwardRegression, and AdditiveModelForwardRegression to use them.
- TAS-832 remove AdditiveSite before replacing with a new version
- TAS-832 Adds a UnaryOperator that shuffles a double[] to BasicShuffler. Changes copy function of CovariateModelEffect so that makes a copy of the covariate instead of just using the original without copying. Part of effort to parallelize StepwiseModelFitter
- Added date/time to logging messages related function progress.
- TAS-794 Initial implementation of Endelman Kinship Matrix Calculation
- TAS-817 Optimization to GCTA Kinship Matrix calculation to completely skip sites where major allele frequency is 100%.
- Fixed bug where projection didn't output correctly (TAS-851). Also automatically makes sure major minor matches when running projection
- Fixes int overflow bug in AbstractFixedEffectLM
- TAS-851: Revert Chromosome comparison processing.  If the chromosome name does not parse to an int, use java's string compare.  The Chromsome.compareNonNumericChrom() method is too slow.
- TAS-850 Adding checks to verify a tag's sequence length minus the barcode length is equal to or greater than maxtag length.  If not, processing stops with an error message.  THis matches processing behavior in GBSSeqToTagDBPlugin.  Without the check, tag building throws a StringOutOfBoundsException when attempt to substring the read to create the sequence.

## (V5.2.12) June 25, 2015

- TAS-817 Implementation of GCTA Kinship Matrix Calculation
- TAS-832 Commented out the chromosomeResidual parameters in StepwiseOLSModelFitterPlugin since they are not implemented yet
- TAS-832 parameter post-processing added to StepwiseOLSModelFitterPlugin
- TAS-832 fixes problems caused by renaming StepwiseOLSModelFitterPlugin2 to StepwiseOLSModelFitterPlugin.
- TAS-832 Move ForwardRegression and implementing classes to Tassel 5. Also, ResamplingGWASPlugin.java added to Tassel.
- TAS-832 renamed StepwiseOLSModelFitterPlugin2 to StepwiseOLSModelFitterPlugin
- TAS-832 removes StepwiseOLSModelFitterDialog
- TAS-832 removes old style plugin for Stepwise
- TAS-832 changes to StepwiseOLSModelFitterPlugin2 needed to implement the pipeline
- TAS-728 Updates to "Taxa Quality Control" Workflow.
- Removing unworking / obsolete class src/net/maizegenetics/analysis/gbs/MergeIdenticalTaxaPlugin.java
- Adds method to ModelEffect for subsampling and adds implementations to FactorModelEffect, CovariateModelEffect, and NestedCovariateModelEffect. Adds a new integerLevel method to ModelEffectUtils.
- Adds method for generating chromosome residuals from a model in StepwiseOLSModelFitter
- Fixed VCF Multithreading Bug.
- TAS-599 Initial Commit of Merging GenotypeTables
- Setting to control the maximum number of cores for LDKKNiImputation
- TAS-817 Add Comments to Implementation of GCTA Kinship Matrix Calculation
- KNNimputationV2 with a full range of options
- KNN Imputation deals with HM3.1 duplications
- TAS-839 Alphabetical sorting of the keys to VCF, plus lambda fun
- TAS-839 Updating VCF Export for INFO exports
- LD methods that permits control on minimum minor allele used for the calculation
- BitUtil to nice looking string
- pre-computes the genotype index for GenotypePheno, which makes the class much more efficient
- Fix LinkageDisequilibrium. getLDForSitePair method instantiation
- Implement runPlugin for FilterTaxaPropertiesPlugin
- fixes bug introduced recently in PartitionedLinearModel which caused the Stepwise analysis to crash

## (V5.2.11) June 4, 2015

- adds constructor to PartitionedLinearModel that takes List as an argument and reformats the class.
- Adds new constructor to SweepFastLinearModel that takes List as an argument.
- TAS-817 Initial implementation of GCTA Kinship Matrix Calculation
- TAS-680 Depth storage has UNKNOWN encoding (0x80) now.  That translates to an integer value of -1.  This reduces the max depth value to 9763.
- TAS-702:  Delivered forgotten file - made TaxaDistBuilder:getDepthMatrixForEncodedDepths public so it is accessible to create a non fixed taxa distribution when getting all tagtaxadist for append mode in GBSSeqToTagDBPlugin.
- TAS=702:  Added code allowing for delete old db, or append to db to GBSSeqToTagDBPlugin step of GBSV2 pipeline.
- TAS-797 Export Kinship to MultiBLUP and LDAK.  Can export in raw and binary form.
- TAS-741: Updated Chromsome comparisons so a 9 is less than 10 and a 9a is less than a 10D.
- TAS-741  Fixed an error in TagDataDQLite:getCutPositionTagTaxaMap(). With non-numeric chromosomes, the value set to "chromosome" in the SQL query must be surrounded by single quotes.
- TAS-794 Refactoring Kinship and DistanceMatrix Calculations.
- Removes getCitation() over-ride from ImputationPlugin that return null. This causes an error in Tassel v5.2.10, though it did not in earlier versions.
- TAS-820 Fixes issue with changed results for DiscreteSitesTest by increasing TOL to 1e-13 in SweepFast.
- TAS-811 Updated GenotypeTableBuilder.getGenotypeCopyInstance() to accept any GenotypeTable (instead of just FilterGenotypeTable and CombineGenotypeTable).
- TAS-702 Initial code to allow user to clear DB tables at various stages of the GBSv2 pipeline.  SAMToGBSdbPlugin, DiscoverySNPCallerPluginV2 and SNPQualityProfilerPlugin have been handled.  GBSSeqToTagDBPlugin has complications that need further attention.
- TAS-762 : Fixed the regex for bad query words
- TAS-762 : Plugin works except for the regex-based check for bad query words (such as "DROP").
- TAS-779 Tripled the speed of KNN imputation, add some poor quality reporting
- Change TOL in SweepFast to 1e-14, because the previous value (1e-12) was too high in one case.
- TAS-789 Add simple guess loading of genotype files.

## (V5.2.10) May 21, 2015

- TAS-668 Adds MultiDimensionalScalingPlugin and ClassicMds to implement MDS. Modifies TASSELMainFrame to add the plugin to the Tassel analysis menu.
- Added code to automatically select data sets after creation.
- TAS-741:  Changes to accomodate non-numeric chromosome processing.
- TAS-769 : Upgrade to PostgreSQL 9.4
- TAS-722 Fixed problem of SAMToGBSdbPlugin not identifying tags when an aligner (e.g. BWA-MEM) performs hard clipping.  TagExportToFastQ now uses the original tag sequence as the "description".  When processing the .sam file, if the aligner presented sequence is not present as a tag in the db, we check if the original sequence is present.  If yes, create a tag from the original sequence, store this tag and position.
- TAS-788 Fixed use of Chromosome entered incorrectly in last update.
- TAS-808 Modified CompressedMLMusingDoubleMatrix to handle duplicate taxa in the phenotype data. In particular it now runs when using mdp_phenotype.txt in the tutorial data set.
- TAS-788 Changed ThinSitesByPosition to take a genotype table rather than file as input.  Added this plugin to the GUI "data" menu tab.
- TAS-786 Changed GenosToABHPlugin to take a GenotypeTable as input instead of a file.  Also added this plugin to the menu.
- TAS-741: added initial methods needed for converting chromosome processing from int to string.
- TAS-795 Created SelectFromAvailableSitesDialog and modified FilterSiteNamePlugin to complete this task. Changed some private fields to protected in SelectFromAvailableDialog to make it easier to extend.
- TAS-770 Converted GBSUtils methods to public so that can be used by lab efforts
- TAS-779 Cleaning up LDKNNiImputationPlugin copies and reinstate tests on distance
- TAS-771 Improve error wording for GenomeSequenceBuilder
- TAS-804	Fixed whitespace splitter bug for VCF headers.  Now follows VCF.
- TAS-794 Updates to Kinship method calculateKinshipFromMarkers to improve efficiency a bit. For example, replaced DoubleMatrix with double[].
- TAS-770: Created GBSUtils.java file to hold methods/constants used by multiple classes in the GBSv2 pipeline.  Methods/constants were moved from GBSSeqToTagDBPlugin to this new utility class.
- TAS-728 Fixed Progress Bars for Workflows.
- TAS-728 Minor corrections to MLM workflows.
- TAS-789 Add simple guess loading of genotype files.
- TAS-632 Corrected IBSDistanceMatrix.computeHetBitDistances() to mask all bit sets
- TAS-779 Initial implementation of MaskGenotypePlugin
- TAS-728 Added ability to find Workflow configuration files (.xml) in source code (when running with no sTASSEL.jar).
- TAS-728 Added "Taxa Quality Control" Workflow.
- TAS-696 Add code to fail immediately if SNPQuality table already contains data for the taxa subset.
- TAS-790 Correct new IBSDistanceMatrix algorithm with better estimateSize() method.
- TAS-788 Commit plugin that thins positions from genotype table based on user specified minimum distance.
- TAS-779 Improved documentation of KNN imputation
- TAS-781 Added in command line functionality
- TAS-786 Added citation for GenosToABHPlugin
- commented out debug code that prints Phenotypes in GenomicSelectionPlugin
- TAS-779 GUI integration of imputation
- replace RidgeRegressionEmmaPlugin with GenomicSelectionPlugin
- Added prediction with PEV output to GenomicSelectionPlugin
- TAS-779 Imputation accuracy improvements by using a different distance matrix calculation
- TAS-781 Added in Full Weighting Functionality and set the MLM and WeightedMLM to use the new plugin style
- TAS-786 Changed input parameters "inputFile" and "outputFile" to be -i and -o respectively for consistency with other TASSEL plugins.
- TAS-786 Fixed imports
- Added pev calculation to EMMAforDoubleMatrix. Added genomic prediction function to GenomicSelectionPlugin.
- TAS-784 Refactored TableReportBuilder.addElements().
- TAS-786:  COmmit code for GenosToABHPlugin.  Code written during hacakthon by Stefan Reuscher, Jeff Glaubitz and Lynn Johnson
- TAS-783 Remove 2D plot from GUI and codebase
- debug changes for GenomicSelectionPlugin and EMMAforDoubleMatrix
- TAS-728 Add citation to Workflow Dialog
- TAS-728 Add Workflow name to top of dialog.
- TAS-785 Add Citation to Plugin Dialog
- TAS-781 Added in Weighting for Training Phase of MLM.
- Modified GenomicSelectionPlugin to delete missing data for each phenotype individually prior to calling EMMA
- TAS-779 Imputation Accuracy evaluation updates to TableReport
- Method for alphabetically sorting genotypes
- Addition of variable length row to TableReportBuilder
- Modified EMMAforDoubleMatrix and added GenomicSelectionPlugin for genomic selection cross validation
- Updata java doc comments for DistanceMatrixUtils
- TAS-779 Initial imputation accuracy plugin
- TAS-728 Added -filterTaxaProperties flag to Tassel Pipeline
- TAS-728 Added -archaeopteryx flag to Tassel Pipeline
- TAS-728 Add -filterTaxaNames flag to Tassel Pipeline
- Added java doc to keepTaxa in DistanceMatrixUtils
- Modified keepTaxa methods in DistanceMatrixUtils
- Added method to DistanceMatrixUtils for selecting a subset of taxa from a distance matrix
- TAS-781 Adding WeightedPlugin to Analysis Menu
- TAS-781 Initial Commit of Weighted MLM Plugin
- TAS-779 LD-KNNi imputation algorithm working with TASSEL input system
- TAS-779 Functional version of LD-KNNi imputation algorithm
- TAS-728 Implement WorkflowPlugin.getCitation().
- TAS-742 Initial implementation of distance matrix calculation that parallelizes on sites.
- TAS-779 Added comparable to Tuple
- TAS-728 Add xml workflow configuration files for MLM with and without (population structure and kinship)
- TAS-728 Refactored ManhattanDisplayPlugin to select TableReport that has p values when given multiple TableReports.
- TAS-728 Removed excludeLastTrait from GLM (PCA) workflow.
- TAS-728 Add Workflow Menu to Tassel GUI
- TAS-728 Add xml workflow configuration files for GLM with and without population structure
- TAS-779 Initial LD KNN plugin stub
- TAS-728 Added code to ManhattanDisplayPlugin to select table report with p-value.
- TAS-728 New -mhd (Manhattan Display) flag to Tassel Pipeline
- TAS-728 Added new -pca flag to Tassel Pipeline.

## (V5.2.9) April 23, 2015

- TAS-758 Added column TagAsForwardStrand, and added cut position offset. This makes it easier for user to verify SNPs as they are relative to the forward strand.
- TAS-772 Support for lower case letters.  Conversion methods for byte arrays.
- TAS-771 Genome index system for genome sequence, filter for characters
- BaseEncoder improvement to convert byte[]
- TAS-758 Changed snp position out file to show "ForwardStrand" based on "forward" annotation rather than "strand", which seems to always be 1. Printing "forward" from the annotation gives us a true/false indicating if the tag was reverse complimented.
- TAS-758 Added cut position and strand to the SNP output file.
- TAS-489 Fixed bug with Flowcell_Lane annotation in GenotypeTableBuilder.
- TAS-758 Adjustments to tab-deliminited output file columns as per Ed's suggestions for SNPCutPosTagVerificationPlugin
- TAS-758 Added new plugin for debugging GBSv2 db data. This plugin allows used to specify a snp or cut position.  GIven this position, a tab-delimited output file is created showing the tag/taxon that includes this position (for Cut Position) and the Alleles with their tag/taxa distribution (for SNP position).  Output is written to file specified by user.
- TAS-708 Move callGenotype() processing into the batch loop.  This allows for clearing the tagCntMap after each batch is processed - we no longer need to keep the entire map of tags in memory.
- TAS-728 Added Workflow Menu
- TAS-735 Updated sqlite jdbc driver to sqlite-jdbc-3.8.5-pre1.jar
- TAS-764 Added in handling of no minor allele for SNPQualityProfiler. Just skipped making stats in this case.  Also improved speed of db writing.
- TAS-748 Fixed error when no minimum quality score was stored for positions.  If requested score is 0, we skip to method that grabs all positions. In all cases, ProductionSNPCallerPluginV2 now checks size of positionList after the database call.  If the list is empty a message is printed indiated no positions were found with a quality score at or above the specified value.
- TAS-708: Added minimum quality score for tag reads.  Changed existing "minimum quality score" to be minimum position quality score.  Made the parameter name for tag read quality score "mnQS" to match GBSSeqToTagDBPlugin.
- TAS-742 Speed improvement for the IBSDistance - Removed LongAdder counter.

## (V5.2.8) April 9, 2015

- TAS-735 Fixes the hard drive full issue from SQLite.
- TAS-735 Updated sqlite jdbc driver to sqlite-jdbc-3.8.7.jar
- TAS-708:  Remove extraneous debug code and duplicate batch execute command from TagDataSQLite:putSNPPositionQS
- TAS-754 Stepwise model following new plugin model
- Concept for making a TableReport implement toString()
- TAS-428 Added warning in Tassel pipeline to use NumericalGenotypePlugin instead of -numericalGenoTransform
- Refactored ExportPlugin to not create Swing component when running from command line.
- TAS-690 Fixed exporting Reference Probabilities where Float.NaN exported as ? (now NA).
- TAS-761  To DiscoverySNPCallerPluginV2: Added code to prevent storing reference tags in the db that were created from the reference genome. Tags flagged as reference that were pre-existing remain in the db.
- TAS-744  Implemented the "depthOutput" option, changing the parameter from "ndo" to "do" and clarifying the description. It defaults to "true".
- TAS-742 Speed improvement for the IBSDistance - increase cache and parallelize
- TAS-336 Added method to Plugins that resets all parameters to default values.
- TAS-708 Added batching to ProductionSNPCallerPluginV2.  Also changed the culledFiles() method of GBSSeqToTagDBPlugin to use stream without a need to sort.
- TAS-378 Implement CombineGenotypeTable.positions()
- TAS-759 Added summary output of tags processed for specified quality score to end of ProductionSNPCallerPluginV2.
- TAS-757:  Added checking for file size > 1 before sorting GBSSeqToTagDBPlugin fastq files in method culledFiles().
- TAS-757  Added alphabetical sorting to the file list prior to its return in "culledFiles()" method.
- Simplification of adding annotations to Positions
- Refactored FileLoadPlugin constructor for simpler use.
- TAS-689 Added in delimiter support to Synonymizer. This allows for handling of long taxa names which should be split into shorter names.
- TAS-415 Improved loading of PLINK files.
- Add object version of adding annotations to GeneralPosition
- Implement string() for GenomeFeature
- Implement indexOf for PositionHDF5List
- TAS-748 Added minimumQualityScore as a parameter to ProductionSNPCallerPluginV2. Default is 0 (ie, grab all snps from snpposition table).  Added checking to the tagDataReader.getAlleleMap() stream to ensure allele position exists in the PositionList obtained from the database.
- TAS-748 Removed commented out code.
- TAS-748 Changed  method TagDataSQLite:putSNPPositionQS() to take a positionList that includes a quality score annotation.  This annotation is stored in the snpposition table.  Moved readQualityScoreFile() out of UpdateSNPPositionQUalityPlugin and into PositionListIOUtils.  IT now returns a PositionList that stores the quality score in an annotation.
- Also added new signature to getSNPPOsitions() that takes a qualityscore value and returns only those positions from the snpposition table whose quality score meets or exceeds the specified value.
- TAS-751 Make GenotypeSummaryPlugin friendly for code level use
- Improved appearance of menu icons.
- Added GenotypeTableUtils.calcBitPresenceOfDiploidValueFromGenotype() to support project at Iowa State.
- TAS-748 New plugin UpdateSNPPositionQualityPlugin.java added that takes a csv file and GBSv2 database as input.  Quality score from the csv file is stored in the snpposition table in the database for the specified chromosome/position.
- TAS-748 Added Chromosome string to the csv output for SNPQualityProfiler
- TAS-728 Adding simple workflow configuration.
- TAS-600 Conversion of Synonymizer to use TaxaList structures/Fixing Bugs
- TAS-747 Create a collector for TaxaList.
- TAS-747 Create a collector for PositionList.
- TAS-715: Created synchronized maps, printed Read Counts values to output log for ProductionSNPCallerPluginV2
- TAS-728 Adding Dialog for Workflows in Tassel GUI.

## (V5.2.7) March 26, 2015

- Fixed non in-memory VCF build bug.
- TAS-709 Added message to DiscoverySNPCallerPluginV2 indicating which chromosome is being processed.
- TAS-738 Several changes to GBSv2SeqToTagDBPluginV2 to handle errors when processing fastQ files in a directory.  Code now checks for 0 tags added to map to prevent a "divide by 0" error.  It also "culls" the files in the directory, creating a list for processing that only contains files represented in the key file.
- TAS-740 Added command-line option '-ldTestSiteName' and functions to do site-by-all LD based on a site's SNP ID
- Fixed bug where information does not write to DB.
- TAS-728 Added method Plugin.wasCancelled() to indicate whether plugin was cancelled by user.
- Added new enzymes to GBS Pipeline: AseI, AvaII, and KpnI-MspI
- TAS-730 VCF import, removed timeout as the Futures will block and wait till the thread finishes.
- TAS-728 Initial Workflow Plugin
- TAS-713: Added calls to "sortTaxaAlphabetically()" when building the masterTaxaList.  Also alphabetically sorted the ArrayList collection from readTaxaAnnotationFileAL().  This provides consistency with legacy GBS 5, and provides alphabetized taxa when displayed through the HDFView program.  It does not change the order in which the TASSEL GUI displays the taxa when loading an HDF5 file.
- TAS-733 fixes output from RandomGenotypeImputationPlugin
- TAS-717: Added support for a reference genome to DiscoverySNPCallerPluginV2.
- TAS-720 Added in better reporting of statistics.  Prints out a message every 10k and allows for user to output information to csv.
- TAS-736:  Renamed all database parameters for GBSv2 to "db".  They were previously inconsistently accessed as -i, -o, -db.
- TAS-733 Creates RandomGenotypeImputationPlugin that imputes missing genotypes by randomly drawing a genotype from the genotype distribution for the genotype site.
- Improved Datum.equals() implementation
- TAS-578 Initial implementation of export / import TaxaList
- Refactored FileLoadPlugin to prevent creation of GUI components when executing in non-interactive mode.
- TAS-723:  Created new method TaxzListIOUtils "readTaxaAnnotationFileAL" that returns an arrayList rather than a TaxaList.  This allows duplicate taxa in the same lane.  The existing method returning a TaxaList is still available and still called to store unique taxa where appropriate.
- TAS-684 Removed empty VARIANT annotation from GeneralPostion and fixed toString() method when no variants present.
- TAS-579 Added ability to export PositionList to JSON formatted file.
- TAS-731 Fixed RidgeRegressionEmmaPlugin. It was incorrectly testing for the presence of ReferenceProbability.
- TAS-712: Added missing depth information when building taxon from HDF5 file from  ProductionSNPCallerPluginV2.  Code tested by Ramu.
- TAS-728 Initial Infrastructure to TasselPipeline for adding Workflows to Tassel GUI
- TAS-579 Added ability to import PositionList from JSON formatted file.
- TAS-667 modified PrincipalComponentsPlugin to use imputation by mean when the input data is a Genotype without a ReferenceProbability.
- Refactored Sizeof Utility
- TAS-725 Modified CompressedMLMUsingDoubleMatrix and EMMAforDoubleMatrix to calculate the additive effect as the difference between the homozygous classes divided by two.
- TAS-727 Added in Apache 2.0 license to the code which was missing it. Also cited the source of the code and the author.
- TAS-720 Added in position counter which outputs every 10k positions and display of initial parameters
- TAS-510 Added Genotype Stream for given Taxon.
- TAS-726 Fix Chart display of PCA results.
- Removing unused classes src/net/maizegenetics/analysis/association/DiscreteSitesFullModelOnlyFELM.java and src/net/maizegenetics/analysis/association/DiscreteSitesMatrixFELM.java
- TAS-699: Added maxTagLength parameter to ProductionSNPCallerPLuginV2. This value needs to match the maxTagLength used when creating the database or incorrect allele values will be stored in the hdf5 file.
- TAS-721 Added Strand constants to Position Interface.
- TAS-656 Initial Commit of Database Stream.  Tested using the GBSv2.db. Might need some more validation checks.

## (V5.2.6) March 12, 2015

- TAS-711 Fix Manhattan and QQ Plots broken because of MLM column name change.
- TAS-710: Added maxTagsPerCutSite user option to DiscoverySNPCallerPluginV2.  This option allows the user to specify how many tags may be associated with a position for alignment.  Too many tags and biojava getMultipleSequenceAligment() grinds to a halt. Default is equal to the maximum tag length (which if not defined by user, is 64)
- TAS-701: Corrected last putTagAlignment metric query to use ">" instead of "="
- TAS-701:  Added metrics to TagDataSQLite:putTagAlignments indicating number of cut sites, number of cut sites with  1, 2, 3 and > 3 tags.
- TAS-489 Added import / export of TaxaList to / from JSON.
- TAS-40 Changed Site Scores of type float stored as byte to round float values to 3 decimal places.  This reduces file size when exporting.
- TAS-690 Add Ability to Export ReferenceProbability from GenotypeTable on Command Line (-exportType ReferenceProbability)
- TAS-590 Added ability to get allele presence for Diploid Unknown (WHICH_ALLELE.Unknown).
- TAS-526 Fixed problems filtering GenotypeTable that doesn't have genotypes.
- TAS-657 Show UNEAK Plugins in Tassel 5 as Deprecated to using Reference Pipeline

## (V5.2.5) March 5, 2015

- change the sign of the dominance effects in MLM
- TAS-701 Removed hard-coded index value for AS aligner field.  This value is optional and bowtie2/bwa put it in different locations.  New method is written to find and return AS value if it exists.  Default to 0 when non-existent.
- Enable save to file option in AbstractFixedEffectLM. The code was not writing the results to a file
- TAS-690 Add Ability to Export ReferenceProbability from GenotypeTable on Command Line (-exportType ReferenceProbability)
- TAS-602 Updating Identifier Synonymizer to display similarity scores correctly.
- TAS-698 Fixed exitLimit parsing for StepwiseOLSModelFitterPlugin
- TAS-690 Add Ability to Export ReferenceProbability from GenotypeTable
- [TAS-602] New methods for Synonymizer scoring validated by unit tests
- TAS-691 SAMToGBSdbPlugin Sam file doesn't incorporate into tagDB in GBSV2 pipeline The issue here is that when the sequence is backward aligned and it is shorter than minAlignLength, the returned Tuple(tag,Optional.of(position)) doesn't match any tag in the DB, since it is returned before getting reverse complementary sequence.
- TAS-642 Added command line flag -exportIncludeAnno to indicate whether to include annotations in exported file.
- TAS-688 CompressedMLMusingDoubleMatrix and EMMAforDoubleMatrix modified to decompose marker effects into additive and dominance components. ModelEffectUtils method added to create a genotype effect with genotype classes sorted with homozygous classes first.
- TAS-642 Added option to export Hapmap files without the Taxa annotations.
- TAS-669 Corrected action when canceling file chooser during export.  Also, added message whether to overwrite files if they already exist.
- TAS-488 Continue improving Hapmap Importing Error Messages.
- TAS-686 Added helpful message if use selects a factor index that isn't actually a factor
- TAS-679 Fixed Poor allele cache performance in core genotype table with randomly order sites
- TAS-510 Initial implementation of Byte Stream for Genotypes.
- TAS-680 Depth storage has UNKNOWN encoding (0x80) now.  That translates to an integer value of -1.  This reduces the max depth value to 9763.
- TAS-682:  Change mappingapproach annotation for bowtie2 to be "Bowtie2" for consistency with expectations.
- TAS-681 fixes bug caused by error in factor merging code in PhenotypeBuilder and changes handling of missing in CategoricalAttribute
- TAS-659 Completed handling of CompoundNotFoundException for files assigned to Lynn, Jeff and Fei.  After discussing with Jeff, we decided any occurrence of this exception should result in the halting of processing as it implies there is a prior code error.  In other words, this really should NOT happen and indicates the code stored something other than the expected DNA sequence in the "tag".
- Code was tested via the junits for DiscoverySNPCallerPluginV2 and EvaluateSNPCallQualityOfPipelineTest.
- Removed the maxDivergence option since has not been properly implemented and is not recommended.
- Improved the usage statement regarding the ref option
- TAS-659 Added code to DiscoverySNPCallerPluginV2 to handle CompoundNotFound error thrown from biojava 4
- TAS-659 Adding slf4j-api-1.7.10.jar and slf4j-simple-1.7.10.jar that the new new biojava 4.0.0 depend.
- TAS-659 Update biojava libraries to 4.0.0
- TAS-658 Changed "Remove Minor Allele States" function in filter sites to make both alleles UNKNOWN instead of possibly only one being UNKNOWN.
- GBSSeqToTagDBPlugin modify parameters of maxTagNumber and batchSize
- GBSSeqToTagDBPlugin remove sortFileBySize
- GBSSeqToTagDBPlugin remove reduceTagByCount

## (V5.2.4) February 19, 2015

- TAS-675 Fixed command line where Tassel 4 to 5 HDF5 conversion causes hang up without GUI
- GBSSeqToTagDBPlugin Stop using Trove library. The pipeline works fine now. Looking for other libraries to implement ShortIntMap
- TAS-652 Update Commons Math Library to commons-math3-3.4.1.jar
- TAS-670 Corrected problem loading NaN (unknown) values in matrix files.  Also added user friendly error messages.
- TAS-511 Added EqtAssociationPlugin (Fast Association) to the Analysis menu. Output added to Association folder of data tree.
- TAS-489 Initial redesign of General Annotations
- TAS-572 Adding a filter for aligned tags in DiscoverySNPCallerPluginV2. A new plugin parameter myGapAlignmentThreshold (a double - between 0-1) is added to allow a user to specify a threshold value for tag alignments.  A ratio of indel contrasts and non-indel contrasts is used to calculate the number of gaps that differentiate the aligned tag from the reference tag.  A gap ratio that exceeds the defined threshold results in all tags from that loci being excluded from SNP analysis.
- TAS-637 Updated GUI frames to have a full end-to-end workflow
- TAS-511 icon for EqtlAssociationPlugin
- TAS-511 Adds functions to GenotypePhenotype to handle some additional genotype return types. Makes changes to EqtlAssociationPlugin to handle union joins properly.
- GBSSEqToTagDBPlugin bug fix at removing tags without duplicates
- GBSSEqToTagDBPlugin Set max tag number at beginning, map entry memory size = 100 + 2* taxaNumber. Please set -Xms > 4 times of map memory, -Xmx should be slightly larger than -Xms. Using large memory because 1. prevent rescaling in tagMap (slowing down at least 10 times if that happens). 2. Memory using keep increasing somehow, I think this is due to inefficiency of Java memory management.
- 8 (or 16) fastq files as a batch. At the end of batch processing, reduce tagsWithoutReplicate if map size > 50% max tag number, reduce map size by count if map size = max tag number. System.gc().
- after all the files are processed, reduce map size by minCount first. Then cut the second cut site. Write map to memory.
- Note: make sure the max tag number is large enough to capture rare alleles and prevent map rescaling. More fastqs in a batch will need larger tag number. Considering fix the batch number as 8.
- TAS-511 changes to EqtlAssociationPlugin and SolveByOrthogonalizing so that eqtl methods pass junit tests
- added test for p
- TAS-664 Made changes to GenotypePhenotype fix bug that caused problems when a taxon has a phenotype but no genotype
- TAS-571: Previously the mappingAnnotation defaulted to bowtie2. Code is added to detect "bowtie2" in fastq header.  "mappingAnnotation" is set to bowtie2 or bwa based on presence of "bowtie2" in input file's header lines.
- fixes problem with Kinship calculation from ReferenceProbability. The Kinship class previously did not recognize ReferenceProbability and use the proper method for calculating kinship from it.
- TAS-663 Fixed bug on initializing TagDistributionMap size
- improve Kinship - Endelman method efficiency by a factor of 2
- TAS-662 Fixed a bug where FILLINIMputationPlugin wouldn't accept files as donor inputs (like is done for projection). I had streamlined the inputs/outputs earlier and generated this bug. Now, donor input is a string that should be either a file or folder.
- README.md edited online with Bitbucket
- TAS-660 Add commons-codec-1.10.jar to supporting libraries
- TAS-489 Consolidating reading / writing annotations to HDF5.
- TAS-578 Added support to TaxaListTableReport to show multiple values an annotation separated by commas.
- TAS-649 Need to move when locking Taxa module in HDF5, as building the HDF5 file reads the taxa module during the building.
- TAS-649 Added error message when trying to open HDF5 produced by ProductionSNPCallerPluginV2 -ko option
- TAS-503 Increasing the flexibility of the SNPProfiler
- TAS-511 updates to address junit test issues in EqtlAssociationPlugin and SolveByOrthogonalizing
- TAS-503:  snpQualityInsertPS prepareStatement updated to write to "minorAlleleFreqGE2" (vs "minorAlleleFreq").  This problem found when running the GBSv2 pipeline test from EvaluateSNPCallQualityOfPipelineTest.java.
- TAS-489 Removing unused class src/net/maizegenetics/util/AbstractAnnotation.java
- TAS-503 Quality scoring of SNPs and reporting to the database - statistic still need to be checked
- Message displayed using JOptionPane in NumericalGenotypePlugin now only displayed when isInteractive = true. Otherwise, logging is used.
- TAS-635
- TAS-635, TAS-643, TAS-644
- TAS-640 fixes bug in CompressedMLMusingDoubleMatrix that failed to line up genotype and phenotype properly when the order of taxa in the data sets was different.
- TagsOnGeneticMap, output FastA for alignment
- TAS-391 Added code to ant build script that includes .sql files in sTASSEL.jar.
- TAS-511 made methods parallel, added minR2 test, and moved p-value calculation to updateOutputWithPvalues in EqtlAssociationPlugin, SolveByOrthogonalizing classes
- TAS-578 Changed TaxaListTableReport to show only "taxa" column.  "taxa" and "name" were duplicated before.
- TAS-415 Improved error message when importing PLINK files.
- TAS-598 SynonymizerPlugin with GUI frames added.
- TAS-511 several fixes to EqtlAssociationPlugin and SolveByOrthogonalizing classes
- fixes bug. average = sum of column / number of non-NaN values. Was average = sum of column / number of rows.
- TAS-628 Added IBSDistanceMatrix.computeHetBitDistancesThirdState() written by Kelly.
- TAS-40 Corrected bug with counters in AlleleDepth class.
- TAS-511 add setters and getters to EqtlAssociationPlugin and degree of freedom method to SolveByOrthogonalizing

## (V5.2.3) January 8, 2015

- Added RemoveIndelsForBeaglePlugin to Impute Menu
- TAS-630 Fixed NullPointerException in ProjectionGenotypeCallTable when getting Allele Frequency
- Fixes multiple phenotype problem in StepwiseOLSModelFitter and adds model Rsq to the output
- Fixes bug in NumericalGenotypePlugin introduced because file import now sets Datum comment to null instead of "".
- added code to set minMaf from CallParentAllelesPlugin in BiparentalHaplotypeFinder
- replaced System.outs with Logging in BiparentalHaplotypeFinder. added code to set minMaf in BiparentalHaplotypeFinder from CallParentAllelesPlugin
- TAS-628 Added GenotypeCallTable.thirdAlleleForAllSites() to get second minor allele.
- TAS-629 Fixed MergeGenotypeTables duplicates positions
- TAS-622 Updated snappy jar from 1.1.0.1 to 1.1.1.6
- TAS-625 Show IUPAC Nucleotide Codes in Report Panel when viewing Nucleotide Data
- fixes bug in AbstractFixedEffectLM that resulted in division by zero error when the number of tests was less than 100
- TAS-623 Fixed QQ and Manhattan Plots Fail due to "Position" and "Chr" column name change
- When displaying text on Tassel Main Panel, made scrollbar be at the top.
- TAS-422 Remove net.maizegenetics.util.Report
- TAS-509 Fixing VCF error caused by .:.:.:.:. values
- TAS-615 Rename DiscoveryTBTPlugin to GBSSeqToTagDBPlugin.  TagsByTaxa (TBT) is no longer relevant to the new pipeline - we use a database instead.
- TAS-422 Preparing to remove net.maizegenetics.util.Report
- TAS-584 Remove obsolete IdGroupUtils
- Added Distance Matrix Plugin to Analysis Menu.
- TAS-40 Continued implementation of Site Scores package with added getInstance() methods for use with an existing HDF5 file and more methods for AlleleDepth.

## (V5.2.2) December 18, 2014

- Fixed bug with GenotypeTableBuilder.getInstanceOnlyMajorMinor(), that caused "Remove minor SNP states" from working correctly in Sites Filter.
- Added -e ignore option for unspecified restriction enzyme. In this case barcodes and common adapter start sequences will be recognized, but chimeric DNA fragments (or partial digests) will not be trimmed.
- TAS-381 Cleaned up FisherExact.java, removed commented out lines, altered resizeArrays to make more efficient.
- TAS-378 Continued implementation of CombineGenotypeTable
- write function in TagsOnGeneticMap
- TAS-503 SNPQualityProfilerPlugin now actually does some profiling of the SNP quality.
- TAS-381 FischerExact:resizeArray()  Added missing code to redeclare array with new size before calculating array values.
- TAS-381 Make FisherExact a singleton, add resizing factorial arrays.
- TAS-399 Changed GeneratePluginCode to use variable name for method signatures instead of GUI Name.
- TAS-611 Corrects values used for degrees of freedom in calculation of p-values for additive and dominance terms in DiscreteSitesTest and comments out 3 lines of debug code in AbstractFixedEffectLM
- Added MslI.  Fixed MspI common Y-adapters to match Poland et al. 2012.  Added PstI-MspI-GDFcustom.
- TAS-511 Added attributeStream() to Phenotype and subclasses and added method to remove observations with any missing values to PhenotypeBuilder. The changes were made to make it easier to test EqtlAssociationPlugin but should be useful for other applications as well.
- TAS-511 deleted function using ReferenceProbabilitySpliterator
- TAS-511 further implements EqtlAssociationPlugin and removes ReferenceProbabilitySpliterator
- sort function in TagsOnGeneticMap
- TAS-578 Initial implementation of export / import TaxaList
- fixes bugs with save to file option in CompressedMLMusingDoubleMatrix
- TAS-511 Additional work on EqtlAssociationPlugin implementation.
- TAS-503 Concepts for SNPQualityProfilerPlugin.
- TAS-591 Clarification of methods.  Corrected bug with invariant cutPositions.
- TAS-571 Correcting command line parameters with same prefix.  Changed proportion to go from 0 to 1.
- TAS-591 Requesting comment/clarification for findAlleleByAlignment.  If the syntax of "//TODO TAS-591" is not correct, please let me know.
- TAS-511 start of matrix eqtl plugin implementation
- TAS-511 Fixed errors in SolveByOrthogonalizing so that it passes unit tests. Apply to permutation testing in AbstractFixedEffectLM.
- TAS-488 Continue improving Hapmap Importing Error Messages.
- TAS-563 PluginParameter for choosing from a runtime list.
- TAS-371 Added helpful message to ProductionPipeline for if user points to a directory without subdirectories
- TAS-511 Fix bug and add ability to solve add plus dom models to SolveByOrthogonalizing
- In DiscreteSitesFELM, if F = Infinity (usually because df = 0) set F = NaN and p = NaN.
- fixes bug in PhenotypeBuilder that was failing to set "NaN" as missing data in version 4 style trait import
- Added Updated ProjectPcsAndRunModelSelectionPlugin from Alex
- TAS-573 Improve numerical imputation efficiency
- add test for parentPlugin != null in progress bar update in AbstractFixedEffectPlugin
- implements progress bar in FixedEffectLinearModelPlugin
- TAS-563 Initial implementation: PluginParameter for choosing from a runtime list.
- TAS-586 Fixes bugs in permutation code in AbstractFixedEffectLM
- TAS-511 SolveByOrthogonalizing is the initial implementation of methods for linear model solving outlined in the matrix Eqtl paper.
- TAS-573 adds imputeFromNonMissingNeighbors method which only uses neigbors with non-missing data for the value being imputed.
- TAS-571 Added new optional filter variables to SAMToGBSdbPlugin allowing user to specify minimum alignment length and minimum alignment proportion.  Methods to calculate are called from parseRow().
- TAS-582 Added taxa name to phenotype import error message in PhenotypeBuilder.ImportPhenotypeFile()
- TAS-577 Implements spliterator that iterates through ReferenceProbability by site as ReferenceProbabilitySpliterator. Adds a method to GenotypeTableUtils that uses that to create a stream.
- TAS-582 Make Numerical Import functions consistently accept only NA, NaN, and . as Unknown values.
- TAS-456 Changed Data Menu option name from "Transform" to "Transform Phenotype".  And improved error message when incorrect input selected.
- TAS-456 removed obsolete class src/net/maizegenetics/analysis/numericaltransform/Conversion.java
- TAS-456 Removed obsolete class KNN

## (V5.2.1) November 20, 2014

- add isPAV in TagsOnGeneticMap
- FSFHap improvements. The value for maximum missing is now used. Previously it was not being passed to BiparentalHaplotypeFinder. The variable minNotMissing was changed to minNotMissingProportion so that the amount of not missing data required is a function of window size.
- TAS-521 GBSv2 bug fixed on sample coverage
- TAS-569 fixes bug in MLM that not correctly matching the kinship taxa to the phenotype taxa when some of the phenotype taxa were not in the kinship matrix
- fixes problem bug in PhenotypeBuilder file import that was not recording -999 as missing data.
- fixes a bug that was attempting to cast a float[] to double[] in AbstractFixedEffectLM
- fixed FilterTraitsPlugin bug that causes a failure if traits are filtered twice.

## (V5.2.0) November 13, 2014 (Improved Phenotype Support)

- TAS-533 adds code to FilterTraitsPlugin and TasselPipeline to handle -excludeLastTrait pipeline flag.
- TAS-428 Validated PCA. Made a change to PrinComp to get output to exactly match R method prcomp. Finished implementing
- TAS-565 Override contains for PositionArrayList
- TAS-544 Support different Illumina quality score strategies
- TAS-544 Eliminate tags not found multiple times in flowcell
- TAS-564 Fixed problems with ModifyTBTHDF5Plugin and SeqToTBTHDF5Plugin where checking for optional parameters failed because the value was empty string instead of null.
- TAS-559 Changed new getSNPPositionsForChromsomes() method to return PositionList.  Added parameter validation.
- TAS-559  Added new class PositionListIOUtils with initial method to read the SNP Conserve File.
- TAS-428 fixed bug in PhenotypeBuilder that caused intersect join to give incorrect results depending on the order in which the phenotypes were joined.
- Fixed a typo in the usage statement
- TAS-456 Moved numerical imputation from the Data menu to the Impute menu.
- TAS-555 fix MLM handling of ReferenceProbability
- TAS-559 Code to pull SNPs for a specified chromosome from the database.  This should be under a separate issue, but I thought I was creating an issue and instead got a sub-task.
- TAS-428 correct bug in CompressedMLMusingDoubleMatrix to print correct residuals
- TAS-334 implemented alleleDefinitions() and alleleDefinitions() in TOPMGenotypeTable.
- TAS-40 Added more user friendly error messages to ReferenceProbabilityBuilder and SiteScoreUtils.
- TAS-428 Added code to importing phenotypes to recognize "NA" and more user friendly error message.
- TAS-456 Changed tool tip for TransformDataPlugin.
- TAS-456 Deprecate KNN class.
- TAS-428 added getters and setters to PrincipalComponentsPlugin
- TAS-428 finished implementation of PCA in the PrincipalComponentsPlugin.  Corrected handling of missing neighbor data in kNearestNeighbors.
- TAS-554 Corrected implementation of depth(), dosage(), alleleProbability() and referenceProbability() in FilterGenotypeTable.
- TAS-428 Fixed a type conversion error and add meaningful name and comments to output
- TAS-428 Implement plugins for data transformation: NumericalGenotypePlugin and TransformDataPlugin. Adds them to the Data menu.
- TAS-456 Added Numerical Transform classes from Janu
- TAS-552 Make Float.NaN is ReferenceProbability view has gray color
- TAS-428 made conversion methods public in NumericalGenotypePlugin
- TAS-552 Change color range for heat maps to be easier to read text when values are 0.
- TAS-428 added stream methods to Phenotype to create Streams of attributes of a given ATTRIBUTE_TYPE. Implemented the new methods in CorePhenotype and FilterPhenotype. Added a method to convert a FilterPhenotype to a CorePhenotype.
- TAS-428 Modified TaxaAttribute so that the name is always "Taxa".
- TAS-428 Finished Phenotype export
- TAS-428 Added Transform to data menu, completed TransformDataPlugin and TransformDataDialog, and added code to PhenotypeBuilder to recognize -999 as missing data in phenotype files.
- TAS-428 added code for setting output filename for compatibility with TASSEL 4 pipeline
- TAS-428 did initial testing of MLM code. This fixes bugs in MLMPlugin and CompressedMLMusingDoubleMatrix so that limited tests using the tutorial data run correctly.
- TAS-428 fix a couple of bugs in CompressedMLMusingDoubleMatrix
- TAS-428 removed permformFunction method from KinshipPlugin so that it stopped overriding the AbstractPlugin method.
- TAS-523 Added -log [] to Tassel Pipeline.  Same as -debug except no debug messages.
- TAS-428 Commenting out use of ProjectPcsAndRunModelSelectionPlugin until migrated to new Phenotype package.
- TAS-428 Removing obsolete net.maizegenetics.trait package
- TAS-428 made changes to NumericalTransformPlugin to maintain compatibility with pipeline options.
- TAS-521 Fix DiscoveryTBTPlugin:processFastQ to handle Tag build error consistently as per change in commit 321fd6f
- TAS-428 debug setAlternateMinorAllelesToMinor in NumericGenotypePlugin
- TAS-428 Initial implementation of NumericalGenotypePlugin that creates a ReferenceProbability from nucleotides. Alternative coding schemes including creation of AlleleProbability will be added later. This version should be sufficient for an initial release.
- TAS-428 finished implementing TransformDataPlugin and TransformDataDialog. Still need to write junit tests.
- TAS-428 move TransformDataDialog and TransformDataPlugin from net.maizegenetics.analysis.data package to net.maizegenetics.analysis.data.numericaltransform
- TAS-521 Change SAMToGBSdbPlugin:processData() and parseRow() to ignore Tags that null.  BaseEncoder:getLongFromSequence translates a sequence of 32 T's to value -1.  This sequence was seen in Zea_may.AGPv3 chromosome files.   For now, ignore them.
- TAS-521 Change TagBuilder instance(..) code to create a 0 length seq2Bit array instead of returning null when getLongArrayFromSeq() returns null.
- THis allows calls to TabBuilder.instance(..).build() to succeed without a null pointer exception. (though the result from build() may be null)
- TAS-521 Fix another null pointer exception and clean up unused code.
- TAS-521 Continued fixed serious bugs in the TaxaDistribution not incrementing and failed units tests
- This works for the smaller taxa distributions
- TAS-521:  Fix null pointer exception in TaxaDistBuilder.
- TAS-428 Initial implementation of Phenotype export
- TAS-521 Fixed serious bugs in the TaxaDistribution not incrementing and failed units tests
- TAS-521 Approaches testing map size less frequently
- TAS-521 Recommit of fixes bug in compression approach of TaxaDistributions
- TAS-541 : populated "Effect" from SnpEff output
- TAS-521: Changes to TaxaDistExpandable to improve fastQ processing
- performance.  Trove Library added.
- TAS-456 Removed obsolete numerical transform Conversion methods
- TAS-456 Removed obsolete numerical transform classes
- TAS-456 Corrected class name when getting icon resource in TransformDataPlugin
- TAS-177 Altered FILLIN help notification to mention including all sites
- TAS-521 Fixes bug in compression approach of TaxaDistributions
- If the first taxa was not within the first 256 taxa, it caused a bug.
- TAS-521 Method to handle various fast quality score systems
- TAS-428 Work on conversion to phenotype package. Fixes error related to use of non-existent NumericAttribute.stream() method.
- TAS-428 Add toString() method to PhenotypeAttribute classes that returns name()
- TAS-428 work on implementing data transformation
- Output change in PanA
- Added some new restriction enzymes
- TAS-539 Added harmless comment to test commit
- TAS-539 Altered LD plugin so command-line default is to treat hets as missing
- TAS-531 : minor changes to schema
- TAS-518 Functional but not tested version of ProductionSNPCaller for V2
- TAS-521 Fix bug in DiscoveryTBTPlugin (v2) processFastQ() that halted reading of the file on a null pointer exception.
- TAS-521 one more optimization of DiscoveryTBTPlugin.
- TAS-521 Fixed method header comment from last submission.
- TAS-521:  Reworked BaseEncoder:getSequenceFromLong() to use "append" vs "insert".  This improves execution of DiscoveryTBTPlugin by about 9%.  Tested for accuracy of converted string via junit testing of the new method, and tested performance via the DiscoveryTBTPluginV2 junit and JProfiler.
- TAS-531 : Modified gene table slightly
- TAS-428 Added logging message to Tassel Pipeline to warn users that -glm flags are deprecated if those flags are used.
- TAS-182 Initial implementation for viewing Projection Genotype Table in GUI.
- TAS-518 Initial ProductionSNPCaller for GBSv2
- TAS-507 GenomeSequenceBuilder:  Added range check to chromosomeSequenceIntInt(), added extraneous text from fasta file Chromosome line to the Chromosome object's annotation list.
- TAS-534 : jar for postgresQL jdbc
- TAS-531 : added FOREIGN KEY chromosome_id to qtl_peak table
- TAS-507 Simplification of compression and chromosome naming
- TAS-428 partial implementation of PhenotypeTransformPlugin and TransformDataPlugin
- TAS-507 Changes to return empty set if no chromosomes.
- TAS-507 chromosome sequence with recording of length
- Recording of length is key to simplify the code.
- TGAS-507 Removed unused method from GenomeSequenceBuilder
- TAS-428 adding NumericalGenotypePlugin in place of GenotypeTransformPlugin
- TAS-428 fix bug in PhenotypeBuilder join code that caused joins to fail for some data sets
- TAS-428 added code to PhenotypeBuilder to handle -excludeLastTrait pipeline command
- TAS-507 Additional changes to GenomeSequenceBuilder, renaming ReferenceGenomeSequence class to HalfByteGenomeSequence and moving functionality out of that class, into the builder.
- TAS-532 Updated maxDonorHypotheses help message to indicate that heterozygous taxa need this value to be set very high
- TAS-532 Added pluginparameter -impAllHets for FILLIN to impute as hets even when original call is homozygous. Also renamed donorFile to donorDir (since it's now a directory)
- TAS-428 TransformDataPlugin created to deal with transformations. It will call either the GenotypeTransformPlugin or the PhenotypeTransferPlugin. Most of the implementation has yet to be done.
- TAS-507 Fix for handling unknown allele (N).
- remove src/net/maizegenetics/dna/map/ReferenceGenomeSequence.java
- TAS-507 Adding GenomeSequenceBuilder, moving ReferenceGenomeSequence class inside the builder.
- TAS-428 Debug issues with StepwiseOLSModelFitter conversion
- TAS-468 Corrected ProjectionLoadPlugin to not fire data set return event, as abstract class does that.
- TAS-398 Added TasselLogging.basicLoggingInfo() that prints basic logging message about system.  Used by GUI, CLI, and Logging Dialog.
- TAS-428 modified methods to convert genotypes to reference probability to code major allele as 1 to conform to Tassel 4 coding.
- TAS-428 Deleted two PCA classes no longer being used. Start implementation of PCA plugin.
- TAS-507 More changes to GenomeSequence code.  - Added builder class - Javadoc style comments for GenomeSequence interface - Altered ReferenceGenomeSequence read method to return allele/byte.
- TAS-531 : initial PostgreSQL schema for storing genome annotations
- TAS-428 add methods for converting Genotypes to numerical values
- TAS-507 More changes to ReferenceGenomeSequence to add clarity - Changed chromPositionMap to have key type "Chromosome" instead of String.  - Expanded GenomeSequence interface comments and changed to javadoc style - Removed "getAlleleEncoding" method, replaced with method from NucleotideAlignmentConstants
- TAS-428 added code to handle Tassel 4 style pipeline commands
- StepwiseOLSModelFitter, fixed bug that was crashing analysis
- TAS-428 Modified GenotypePhenotype TableReport to display short genotype Strings when the genotypes are SiteScores.
- TAS-428 Modified exclude all in FilterTraitsDialog so that Taxa are not excluded.
- TAS-428 Added method to set random seed for testing permutation
- TAS-428 Added getters and setters to FixedEffectLMPlugin
- TAS-507 Adding changes I forgot to make on initial submission.  For ReferenceGenomeSequence:chromosomeSequence(): - return should be nothing more than the chromEntry value.
- TAS-428 Fixes bugs in PhenotypeLM and in PhenotypeBuilder uncovered by unit test. PhenotypeBuilder fixes affect fromFile and fromLists methods.
- Initial code for GenomeSequence (interface) and ReferenceGenomeSequence(class).  Some of this will change when the GenomeSequenceBuilder is implemented.  - GenomeSequence.java:  new file defining the interface - ReferenceGenomeSequence.java:  new file implementing the functionality - DiscoverySNPCallerPluginV2.java:  added keyFile functionality, but have NOT added the new plugin value.  Currently only used in junits.  - NucleotideAlignmentConstants: added method to grab allele char from number.  Used in junits when decoding a byte string.
- Still to do:  Implement the GenomeSequenceBuilder class.
- TAS-428 Implement new version self-describing version of FixedEffectLMPlugin that uses Phenotype package
- TAS-428 Correct bugs in handling of Float.Nan in byteToFloatPercentage and floatToBytePercentage
- TAS-428 modify FilterGenotypeTable to handle case where the GenotypeCallTable is null
- TAS-428 add guess function line in FileLoadPlugin to recognize  header as trait file
- TAS-428 fix problems uncovered by unit test validation
- TAS-527 Corrected Linkage Disequilibrium null pointer exception when accumulate R2 is checked
- TAS-428 Converted StepwiseOLSModelFitter and StepwiseOLSModelFitterPlugin to use Phenotype package
- Fixed bug in BiparentalHaplotypeFinder and NucleotideImputationUtils that was causing FSFHap to fail if the first marker in a dataset was all missing data.
- TAS-502 Completed submitting of the discovery pipeline results to the DB
- TAS-428 Convert KinshipPlugin and Kinship class to use Phenotype package. Change KinshipPlugin to self-describing and delete KinshipMethodDialog.
- TAS-428 Convert RidgeRegression to use Phenotype package classes instead of trait package classes. Add static method to AssociationUtils to convert a primitive double array to a float array.
- TAS-514 Add support for GenotypePhenotype when adding GenotypeTable as parameter.  Also add support to limit which components of a GenotypeTable are allowed for that parameter.
- TAS-514 Add support to PluginParameter to select a component (i.e. Genotype, ReferenceProbability, Dosage, Depth, AlleleProbability) from a given GenotypeTable.
- TAS-428 Changed FilterTraitsDialog and FilterTraitsPlugin to use Phenotype package
- TAS-428 Fix bug in GenotypePhenotype. getTableColumnNames() was returning null.
- TAS-453 : added the option keepTempGenotypes.
- So you don't have to re-extract the taxa if you are working with a full set of flowcell lanes for a project (e.g., running maize landrace families from 13FL from Rare Alleles C2.1).
- TAS-428 Made changes to implement fixed effect linear models using the Phenotype package.
- TAS-428 Converted MLM to use Phenotype package and site scores. Moved public static methods from AbstractFixedEffectLM to AssociationUtils.
- TAS-501 Initial support for recording the direction of alignment in GBSdb
- TAS-428 move TableReport creation from instantiating classes to AbstractFixedEffectLM. Fix some bugs uncovered by junit tests.
- TAS-428 Convert UnionAlignmentPlugin to use Phenotype package.
- TAS-524 Updated FSFHap documentation to show how to run it from the command line.
- TAS-428 Convert DataTreePanel to use Phenotype package.
- TAS-428 Converted FilterTaxaAlignmentPlugin to Phenotype package.
- TAS -428 Converted net.maizegenetics.analysis.data package from net.maizegenetics.trait to net.maizegenetics.phenotype
- TAS-501 Adding functionality for tags that do not map, and reverse complement tags
- TAS-523 Add option to send Tassel pipeline to a file (-debug )
- TAS-428 implement PhenotypeLM a class that takes a Phenotype as input then tests whether taxa are significantly different and generates BLUEs for the taxa
- TAS-103 : No longer moves fastq files when they were not actually analyzed (for example, if the ready.txt file was not present)
- TAS-428 continued implementation of GLM classes that use the Phenotype package
- add toArray methods to DoubleMatrix
- TAS-391 Changing of AtomicLong use to LongAdder
- TAS-391 General clarifications of GBSv2 pipeline
- Fixed bug with DiscoveryPipeline minimum threshold
- TAS-511 fixed bugs in DiscreteSitesMatrixFELM, a test implementation of the method from Matrix eQTL.
- TAS-428 changes to FixedEffectLM classes to fix problems and some additions to Phenotype package classes related to help implement FixedEffectLM. FixedEffectLM and its derived classes will replace most of the code in FixedEffectLMPlugin once fully implemented and tested.
- Added code to skip markers that have no non-missing values. Also deletes a class, DiscreteSitesAddDomWithReplicationFELM, that was never actually used.
- TAS-502 Adding database close operations to GBS db
- TAS-432 Initial display of ReferenceProbability in Tassel GUI.
- TAS-334 More implementation of TOPMGenotypeTable
- TAS-500 In ReadNumericMarkerUtils, set physical position same a site number instead of 0.
- TAS-491 Redesigned GenotypeTable to make GenotypeCallTable optional.  New method hasGenotype() indicates whether the GenotypeTable has genotype values.
- TAS-513 Corrected Bug with Genetic Heat Map Colors
- TAS-511 Creates class to implement Matrix eQTL method for OLS solution to linear models
- Adding new classes and modifying LinearModelUtils and ModelEffectUtils to support new phenotype package and to improve organization of the fixed effect LM models.
- Added code to run_pipeline.pl and start_tassel.pl to fix em-dashes when setting minimum and maximum heap size.
- TAS-336 Improved formatting of Plugin error messages sent to the GUI.
- Improved error message of FastqToTagCountPlugin when no fast files found in input directory.
- TAS-502 Completed basic implementation of new SNP caller using DB and lambdas
- FIx parameter list to GenotyeCallTable.getInstanceUnknownValues used by ReadNumericeMarkerUtils.
- TAS-508 Changed Default "Min not Missing" to 0.0 for FilterTaxaPropertiesPlugin.
- TAS-508 Convert FilterTaxaPropertiesPlugin to PlugParameter Design
- Corrected Sequence Viewer to not allow reordering of columns.
- Add Files.lines() that will work with gz files, and reduce repetitive code by setting the buffer the standard size
- TAS-502 Initial redesign of DiscoverySNPCallerPluginV2
- TAS-391 Added reference tag support to Tags
- TAS-391 Fully support for variable length tags
- TAS-504 Trie based barcode processing, encapsulation of enzyme information
- TAS-500 Add handling of numerice marker files as genotypes.
- TAS-502 Copy DiscoverySNPCaller to new package
- TAS-404 Migrate one step TBT class to new package and integrate with TagDataWriter
- TAS-501 Implementation of plugins to work with Aligners - export to Fastq, import SAM to DB
- TAS-391 Improved Tag classes and builders that deal with variable length tags
- TAS-480 Commit of SQLite database adapter
- TAS-480 Initial migration of TagData and SQLite implementation
- Added ProjectPcsAndRunModelSelectionPlugin to code base for Alex to test.  Must comment out one line in FileLoadPlugin to see on GUI.
- TAS-40 Moved ALLELE_DEPTH_TYPES from AlleleDepthBuilder to AlleleDepth and made it public.  This is for easier iteration over an AlleleProbability instance.
- TAS-40 Moved ALLELE_PROBABILITY_TYPES from AlleleProbabiltyBuilder to AlleleProbability and made it public.  This is for easier iteration over an AlleleProbability instance.
- fixed GLM bug in additive only model
- TAS-40 Added ReferenceProbability to GenotypeTable and GenotypeTableBuilder
- TAS-40 Changed SiteScoreUtil to handle Float.NaN.  byte value 255 is translated to Float.NaN

## (V5.1.0) August 28, 2014 (Requires Java 1.8)

- TAS-493 implement separate phenotype on factors and joining when there are multiple observations for some taxa at the same factor levels. The primary purpose of the changes is to support table pivoting.
- Added debugging level logging to export hapmap files.
- TAS-290 Correct XML creation when flag is self describing plugin
- TAS-290 Update TasselPipelineXMLUtil.TAG_STRINGS
- TAS-498 Added option to printing usage statements in ListPlugins plugin
- TAS-498 Add getUsage() to Plugin Interface
- TAS-423 Corrected LinkageDisequilibrium calculation of total number of tests when value greater than max int.
- TAS-486 Corrections to GeneralAnnotationStorage reading / writing annotations to HDF5 when there are no annotations.
- Fastq read QC and Contiging
- TAS-486 : changed the stored annotation names to dataSetName and dataSetDescription. Used the new fxn GeneralAnnotationStorage.writeToHDF5() to write the annotations.
- TAS-486 Modified GeneralAnnotationStorage.readFromHDF5 to be more generic, so that it can be used for other data structures (i.e. Position Lists).  Also added writeToHDF5 methods for generic writing annotations to HDF5.
- TAS-486 : DataSetName & DataSetDescription are now written to and read from an HDF5GenotypeTable
- TAS-492 Improved performance of TableReportBuilder when writing to a file.
- TAS-453 : added new PluginParameter positionSourceHDF5GenoFile
- Reverted last change because I forgot the dash in TAS 453
- TAS 453 : added new PluginParameter positionSourceHDF5GenoFile
- TAS-497 Fixed bug with LD
- TAS-391 Unifying Tag classes with an AbstractTag
- Fixes a TaxaDistribution compression bug involving Unsigned byte maximum.
- TAS-473 Changes to Phenotype package classes to support table pivoting and subsetting phenotypes by factor levels
- Fixes a bug when there is no known variants at a position.
- Start code for sequence reads QC and contiging
- TAS-492 Initial implementation of TableReportBuilder
- TAS-473 Changes to Phenotype package classes to make PhenotypeBuilder more intuitive to use and to partially implement phenotype table pivoting.
- TAS-423 In TableReport...  public int getRowCount(); change to public long getRowCount(); public int getElementCount(); changed to public long getElementCount(); public Object[] getRow(int row); changed to public Object[] getRow(long row); public Object getValueAt(int row, int col); changed to public Object getValueAt(long row, int col);
- TAS-423 Removed methods getTableData() from TableReport Interface and implementations.
- TAS-423 Removed use of TableReport.getTableData() from HistogramPanel
- TAS-423 Removed use of TableReport.getTableData() from TableReportStatCategoryDataset
- TAS-423 Removed use of TableReport.getTableData() from TableReportPieDataset
- Removed use of TableReport.getTableData() from Grid2dDisplayPlugin
- TAS-423 Removed use of TableReport.getTableData() from TableReportXYDataset
- TAS-423 Removed use of TableReport.getTableData() from SynonymizerPlugin.
- TAS-423 Removed use of TableReport.getTableData() from DiversityAnalyses.
- TAS-423 Removed use of TableReport.getTableData() in TableReportBoxWhiskerCatDataset
- TAS-473 Changes to support pivoting phenotype table
- TAS-491 Add Empty GenotypeCallTable when no Genotypes (only scores)
- TAS-423 Removed unnecessary implementation of TableReport from PolymorphismDistribution.
- TAS-486 Starting code for General Annotations of Genotype Tables in HDF5.
- TAS-485 : Now sets the dataSetName and dataSetDescription via the new GenotypeTableBuilder methods
- TAS-428 changes to Phenotype package classes for use with GLM and MLM
- TAS-488 Performance and error messages improvements of Import for Hapmap
- TAS-486 Initial implementation of GeneralAnnotations for GenotypeTable.
- Removed code in TASSELMainApp that tries to set look and feel to windows.  Just letting it be default.  Throws an exception every time I can tell anyway.
- TAS-415 Refactoring improvements to PLINK import.
- TAS-336 Changed convert() method to handle Enum different, so that it compiles in Eclipse.  Apparently shouldn't have compiled under Java in the first place.  Also, checking if String to String conversion, just return input.
- TAS-344 Updated biojava libraries to 4.0.0 (snapshot) and modified forester_1034.jar (Peter fixed so prevent exiting Tassel when closed).  I will replace biojava when official 4.0.0 jars are released.
- TAS-485 : Now replaces __DATE__ with a yyyyMMdd date stamp in the dataSetDescription as well
- TAS-485 : in the output file name, __DATE__ is now replaced with _yyyyMMdd
- TAS-419 Changing title of Tassel javaDoc home page from "TASSEL 5.0: JavaDoc" to "TASSEL JavaDoc 5.0", so that searching for "Tassel Javadoc" ranks higher in the search engines.
- TAS-487 : bug fix so that HDF5Utils.isHDF5GenotypeLocked() works right
- TAS-485 : Almost ready to go, once DataSetName and DataSetDescription annotations can be added via the GenotypeTableBuilder
- TAS-103 Changed AbstractPlugin to fireDataSetReturned() even if performFunction() return null.  This will allow chaining Plugins when no input needed to start next step.
- TAS-485 : Commented out the broken line
- TAS-485 : Start of this simple plugin
- TAS-485 : added a getBuilder(String existingHDF5File) function for an unfinished HDF5 file
- TAS-480 Implementation of TaxaDistribution class that provides roughly 4-fold compression relative to HDF5 TBT
- Removed obsolete unused class net.maizegenetics.util.RightJustifiedFormat
- Removing obsolete unused classes in net.maizegenetics.util Histogram and IntLinearCongrRandom
- TAS-481 Changed Self-Describing Plugins parameter dialog to cancel when clicking X instead of running.
- TAS-483 Fixed TaxaListUtils.getAllTaxa() to return a sorted list
- TAS-481 Self-Describing Plugin: Limited number of parameters are displayed
- TAS-176 fixes a bug that occurred when converting parent calls to nucleotides when the order of taxa in the imputed GenotypeTable is different from that in the original GenotypeTable.
- TAS-103 Corrected use of lock file when running Production Pipeline
- TAS-103 Added Summary Log File to Production Pipeline.  Also changed key file to end with "_key.txt".  Also added -ko option to ProductionSNPCallerPlugin executions.
- TAS-176 debugging and improving FSFHap, specifically fixes a problem with the subpopulation check
- TAS-481 Fixed Self-Describing Plugin: Limited number of parameters are displayed
- TAS-371 : ProductionSNPCallerPlugin now adds the taxa annotation "Status" with a default value of "private"
- TAS-476 Cleaning up HDF5 experiments in preparation for GBSData interface
- TAS-453 : Fixed minor typo
- TAS-453 : Added preferredHaplotypeSize to the runFILLIN() method
- TAS-479 Remove unneeded AbstractDisplayPlugin.FileFormat
- TAS-467 Changed cache size of HDF5ByteGenotypeCallTable back to 1.5 * number of taxa.  Needs to be bigger for getting full sites.
- TAS-453 : Added updateTaxonAnnotations(Taxon origTaxon, Taxon newTaxon)
- TAS-176 FSFHap added to impute menu in GUI.  WritePopulationAlignmentPlugin now returns GenotypeTables instead of writing output to a file.
- TAS-474 Finished implementing FSFHapImputationPlugin
- TAS-336 Added AbstractPlugin.preProcessParameters() method that can be overridden to process parameters before standard handling by AbstractPlugin
- TAS-468 updated ProjectionLoadPlugin to use a GenotypeTable data set from the Data Tree instead of requiring a input filename.
- TAS-467 Changed cache size of HDF5ByteGenotypeCallTable.  It was too high especially if the number of taxa is large and the max available heap size is relatively small.

## (V5.0.9) July 24, 2014

- TAS-465 Removed jar signing from build.  Java Web Start is no longer supported.
- TAS-446 Removing unused imports from GenomeFeatureMapBuilder
- TAS-453 : replaced Taxon.getAnnotation() with Taxon.getTextAnnotation() to fix bug
- TAS-472 : Fixed a bug in addTaxon() that caused a duplicate taxon in TAXA_ORDER when TAXA_ORDER not present
- Moved to the call to createTaxaOrder(h5w) prior to the call to h5w.createGroup(path)
- TAS-176 Implementation of GUI for full sib family imputation
- TAS-474 Improvements to imputation of nucleotides
- TAS-446 Trying to merge GenomeFeatureMap back into master
- TAS-446 Added JSON library and functionality to build the directed graph holding GenomeFeature relations.
- TAS-446 Added some of the export options
- TAS-472 : added createTaxaOrder() for reverse compatibility
- TAS-453 : Added a preferredHaplotypeSize parameter
- TAS-468 Cleaned up ProjectionLoadPlugin to remove unused variables and obsolete comments.
- TAS-446 Added JSON-reading capability and fixed location lookup ranges
- removed instances of getGenotypeCopyInstance from FILLINImputationPlugin until speed fixed; fixed reporting for focus block viterbi in FILLINImputationPlugin
- TAS-176 Initial (partial) implementation of the Plugin for FSFHap Imputation
- TAS-429 debugged conversion to nucleotides and ported conversion method used for imputation paper to WritePopulationAlignment
- TAS-94 Cleaned up HDF5Utils
- TAS-94 Clean up HDF5ByteGenotypeCallTable
- TAS-415 Add Support for PLink Files
- TAS-429 Fixed bugs, made improvements, and added parameters for BiparentalHaplotypeFinder called by CallParentAlleles. This required new methods in HaplotypeCluster and HaplotypeClusterer.
- TAS-455 Corrections to Projection CLI import (-projection)
- TAS-468 Clean up ProjectionLoadPlugin
- TAS-466 Added (incomplete) FilterGenotypeTableBuilder2
- TAS-464 Fix so will output multiple phenotypes to distinct files.
- TAS-464 Added functionality to CompressedMLMusingDoubleMatrix to export residuals to file (if option chosen)
- TAS-446 Fixed location lookup on GenomeFeatureMapBuilder
- added support for FILLIN accuracy by MAF to the plugin
- TAS-453 : Complete first draft of this plugin (not tested yet)
- TAS-461 added tests of additive and dominance effects to GLM
- TAS-446 Minor tweaks to smooth out schema
- TAS-446 Clean out old comments and stuff no longer needed
- TAS-446 Changed schema for holding annotations; all now kept in a common (extensible) HashMap
- TAS-40 Continued implementation of Site Scores based on Peter's suggestions.
- Change TOGM format, adding PAV field
- Modified setVariantDef() and setVariantPosOff() so that they work with ragged arrays
- TAS-446 Added some validation for entered data
- TAS-446 Progress update (adding imports, expanding methods; not functional yet, thought)
- changes made to FSFHap to correct bug in code that filters out sites within tags and to allow processing of multiple genotype files at once so that the -separate command can be used with CallParentAlleles plugin.
- TAS-455 Add import of Projection Genotype Table to Pipeline (CLI)
- TAS-454 Changed GenotypeCallTableBuilder.getInstanceCopy() method to use Futures for multi-threading to avoid TimeoutException.
- TAS-453 : still a work in progress
- Added new restriction enzyme KpnI
- Stub created for this new plugin
- Part of the Automated Production Pipeline
- TAS-427 Added function to get Position List from Genotype Table and make it a stand alone data set on the Data Tree.
- TAS-452 Corrected DynamicBitStorage to wrap bit sets in UnmodifiableBitSets to make GenotypeTable instances unmodifiable in general.
- fixed bug. Removed cast of GenotypeTable returned by FilterGenotypeTable.getInstance to FilterGenotypeTable, since the function sometimes returns a CoreGenotypeTable.
- Removed most of the nano timers, except for the most pertinent ones
- TAS-371 : ProductionSNPCallerPlugin now adds all taxa annotations present in the key file
- fixed bug in permutation code that occurred when the model contained both factors and covariates
- TAS-446 Ported Eli's Graph structures to branch
- Added modeltype to command line parameters for StepwiseOLSModelFitterPlugin
- TAS-391 Initial stubs and libraries for using BerkeleyDB
- TAS-427 Added function to get Taxa List from Genotype Table and make it a stand alone data set on the Data Tree.
- Added depth first search traversal for Graph
- Updated size() methods for weighted edges
- implement maxHetDev option setting so that junit test runs successfully
- Adding directed graph with appropriate utilities with the builder
- Adding basic graph utilities
- Revert "Adding basic graph utilities"
- This reverts commit 75fb45b787386357bd8e4e27b3a28aeb126c9d7e.
- Adding basic graph utilities
- TAS-445 fixes row method in EJMLDoubleMatrix and makes sure that row and column methods return column vectors
- TAS-447 Created SiteNamesAvailableListModel to easily instantiate AbstractAvailableListModel for PositionList
- TAS-447 Created TaxaAvailableListModel to easily instantiate AbstractAvailableListModel for TaxaList
- TAS-448 Create HDF5 Schema Viewer
- TAS-404 : fixed a bug in readTaxaAnnotationFile() related to parsing headers for quantitative annotations
- They can start with "#" or "
- TAS-391 Testing serialization and other methods for GBS objects
- fixed a bug of transformation in PanA machine learning
- TAS-398 Added Save button to Tassel Logging Dialog.
- TAS-439 Fixed printf statement in FastqToTagCountPlugin
- TAS-446 Created the base GenomeFeature and -Map classes, and associated builders
- TAS-444 Added support for making PluginParameter (i.e. JComboBox) to be dependency for enabling/disabling another Parameter.
- TAS-444 Add Dependency on Plugin Parameter with certain value
- added CompareAlignments (deprecated) back into FILLINImputationAccuracy and changed the test to match the new class path
- TAS-103 Updated Production Pipeline to not use Enum when defining PluginParameters.
- TAS-440 Add Defaults button to Parameter Plugin GUI
- TAS-442 Added test to UTagPairToTOPMPlugin to check for empty string in output files, not just null
- TAS-409 Slight tweaks to TagPairToTOPMPlugin - added range limit to pad distance parameter and a warning if
- TAS-443 Changes currPos from int to long so going over INT.max doesn't wrap into negatives. (Source of error.)
- TAS-409 Added better GUI names for TagPairToTOPM parameters and used auto-generation of PluginParameters to create getters and setters
- TAS-441 Added a setParameter(PluginParameter, Object) wrapper to allow setting parameters without having to extract cmdLineName()
- TAS-103 Beginning adjustments to Production Pipeline to match Jeff's new design.
- TAS-439 Corrected logging statement using printf
- Removed unused src/net/maizegenetics/gui/MultiTextRowHeader.java
- readded accuracy options for FILLINImputationPlugin and removed again the FILLINAccuracyPlugin

## (V5.0.8) June 26, 2014

- TAS-398 Added option to TASSEL GUI that redirects logging messages to a window where the user can see logging messages.
- TAS-429 made changes to parameters. made BiparentalHaplotypeCaller the default method
- TAS-404 One step TBT parallelization concepts and initial implementation
- TAS-40 Added Dosage to GenotypeTableBuilder
- TAS-435 Added a few words to info message about not finding TasselBlas to clarify
- TAS-429 Haplotype clustering changes used by BiparentalHaplotypeFinder
- TAS-429 added new -cluster method to improve performance for S1's and residual heterozygosity. This replaces the old cluster method which did not perform will with few generations of selfing
- modified getUsage to provide additional info about about parameters
- TAS-398 Added -debug flag to Tassel Pipeline to set logging level to debug.
- TAS-433 Changed action of pop up menu over Sequence Viewer.  Before had to hold mouse button to keep menu visible.  Now stays until you click again.
- TAS-40 Changed Dosage to return byte since higher value integers can't be stored anyway.
- TAS-40 Added Dosage class to site score package.
- TAS-431 Remove class NucleotideGenotypeTable
- TAS-40 Continued implementation of SiteScores.  Added AlleleDepth which will be a replacement for net.maizegenetics.dna.snp.depth
- TAS-40 Created SiteScore classes for ImputeProbability
- TAS-40 Integrating AlleleProbabilityBuilder with GenotypeTableBuilder
- TAS-389 added javadoc comments for GenotypePhenotype and its builder
- TAS-389 added GenotypePhenotype and GenotypePhenotypeBuilder classes.  GenotypePhenotype will replace both MarkerPhenotype and MarkerPhenotypeAdapter.
- TAS-389 Moved TableReport implements to Phenotype interface
- added CompareAlignments (deprecated) back into FILLINImputationAccuracy
- Added restriction enzyme NspI
- TAS-424 PluginParameters Enabled/Disabled Depending on value of another PluginParameter
- TAS-40 Continued implementation of Site Scores
- TAS-105 Changed command line flag from hetThresh to mxHet in FILLINImputationPlugin.  Seems that is what it was since the test case uses it that way?
- TAS-419 Changed javadoc index.html title to TASSEL 5.0: JavaDoc
- added DEFAULT_NAME to TaxaAttribute
- Corrected HetsToUnknownPlugin.getIcon() method
- TAS-419 Changed javadoc index.html title to TASSEL JavaDoc API 5.0 (${build-date})
- TAS-221 Create public static GenotypeCallTableBuilder getHomozygousInstance(GenotypeCallTable genotype) {
- TAS-425 Refactor GenotypeTableBuilder.copyGenotypeInstance() to use GenotypeCallTableBuilder.getInstanceCopy()
- Corrected SortGenotypeFilePlugin.getIcon() method.
- TAS-378 Continued implementation of CombineGenotypeTable
- TAS-389 debug changes made to pass current version of junit tests.
- TAS-189 Handle negative values as double parameters for command line options
- Improved Data Menu icons
- Improved Load File icon
- TAS-416 Added icon for SortGenotypeFilePlugin
- TAS-426 Added SortGenotypeFilePlugin to Data Menu
- FILLINFindHaplotypes -maxOutMiss parameter typo
- FILLINFindHaplotypes functional
- TAS-40 Added SiteScoreBuilder.getFilteredInstance()
- TAS-40 Fixed byteToFloatPercentage(). It was returning 0's for all values.
- TAS-389 Finished implementing FilterPhenotype.
- finished implementing phenotype join methods in PhenotypeBuilder
- TAS-105 Changed order of FILLIN plugins to have first used to be first on menu.
- TAS-105 Created Impute menu on Tassel GUI for FILLIN and future imputation methods.
- GUI support for FILLIN
- Implementing Phenotype merging in PhenotypeBuilder. In progress.
- TAS-105 Added FILLINFindHaplotypesPlugin to Tassel Data Menu
- TAS-105 Added FILLINImputationPlugin to Tassel Data Menu
- TAS-408 Add asterisk to required Plugin Parameters in Dialog
- TAS-336 Better handling of PluginParameter verification checks.
- TAS-40 Continued implementation of Site Scores
- TAS-421 Change GenotypeCallTableBuilder.setBaseRangeForTaxon() implementation to use arraycopy()
- TAS-398 Added utilities for setting logging level to DEBUG
- TAS-163 Corrected implementation of FilterGenotypeTable.alleleDefinitions()
- TAS-418 Corrected calls to siteOfPhysicalPosition() in FilterAlignmentPlugin to pass in first chromosome.  In the past, null input meant to use first chromosome.  In these situations in FilterAlgnmentPlugin, there are only one chromosome when filtering on position.
- TAS-404 Concurrent version of DiscoveryPlugin
- TAS-404 Changed new HashMap to new HashMap
- TAS-404 Continued progress on one step TBT - parallelization
- TAS-404 Expanding the annotation so that a key file can be processed as annotated taxa list
- TAS-420 Add Glob support to DirectoryCrawler
- TAS-418 Corrected call to siteOfPhysicalPosition() in FilterAlignmentPlugin to pass in Chromosome instead of null.
- TAS-410 Improve parameter checking in DiscoverySNPCallerPlugin
- TAS-399 Added code to generate runPlugin convenience method in Plugins
- TAS-404 Improved documentation for TBT one step process
- TAS 335 : Converted to new PluginParameter design
- TAS-416 Added SortGenotypeFilePlugin (with minor changes to BuilderFromHapmap/FromVCF)
- TAS-389 Changes made to implement filtering in PhenotypeBuilder
- TAS-333 : Converted AddReferenceAlleleToHDF5Plugin to new PluginParameter design

## (V5.0.7) June 5, 2014

- TAS-389 Finished implementing file imports for PhenotypeBuilder
- TAS-403 Convert ModifyTBTHDF5Plugin to new Plugin Parameter Design
- TAS-40 Abstracted data storage from Allele Depth module and created new package net.maizegenetics.dna.snp.byte3d.  Started implementation of Site Scores using byte3d as data storage.
- TAS-103 Updated ProductionPipeline to use recent improvements to Plugin Parameter design.
- Implementing file import methods
- TAS-331 Initial implementation of FilterGenotypeTableBuilder.
- Clean up allele depth classes.
- TAS-392 : Now uses File.getName() or File.getCanonicalFile.getParent() to retrieve file names or parent directory
- Instead of using substring on the String version of the file path.
- TAS-411 Convert ProductionSNPCallerPlugin to new Plugin Parameter Design
- TAS-410 Convert DiscoverySNPCallerPlugin to new Plugin Parameter Design
- TAS-398 Convert Tassel Pipeline, Production Pipeline, and Tassel Main App to use new Logging Utilities.
- TAS-398 Added new utilities class for logging.
- TAS-412 Added HDF5SummaryPlugin to get summary data from HDF5 files (fixed typo of duplicate option flags)
- TAS-412 Minor fix to HDF5SummaryPlugin to remove performFunction() method
- TAS-412 Added HDF5SummaryPlugin to get summary data from HDF5 files
- remove try {} for debugging  in PanA
- set multiple threading  in PanA
- Position block bug fixed in PanA
- TAS-407 Adjusted layout of UNEAK and GBS menus
- TAS-397 Add Tool Tip to Plugin Parameter Dialogs
- TAS-407 Add UNEAK Menu to GUI
- TAS-406 Updated UNEAK plugins to use the new PluginParameter functionality
- TAS-405 Fix ExportUtils.WriteToVCF so it just outputs a summary of sites with no alleles.
- TAS-375 Added javadoc to BinaryToTextPlugin getters and setters
- TAS-396 Convert SeqToTBTHDF5Plugin to new Plugin Parameter Design
- TAS-391 In Utils.getBufferedReader(), set buffer size of gzip stream also to specified buffer size.
- TAS-336 Added verification that inFile() is a file (i.e. not a directory).
- TAS-399 Create Code Generator for Plugin Parameter Getters and Setters
- TAS-394 Added method for Parameter Plugin Description
- Added restriction enzymes NlaIII and SphI
- TAS-333 : Added some feedback on progress
- TAS-390 Convert SAMConverterPlugin to new Plugin Parameter Design
- TAS-336 Added Help Tab to Automatically generated dialogs for PluginParameters.
- TAS-203 Improved documentation
- TAS-393 PluginParameters now check that command-line names are unique
- TAS-336 Added getParameterFields() method to AbstractPlugin.
- TAS-389 Partial implementation of PhenotypeBuilder
- TAS-391 - Stub with todos for the revision of fast parallel fastq reader
- various bugs fixed to get the basic import test running.
- Add PanAReadToKmerPlugin and a few bug fixes in PanA
- TAS-389 partial implementation of PhenotypeBuilder with modifications to Categorical Attribute.
- TAS-387 Convert TagCountToFastqPlugin to new Plugin Parameter Design
- TAS-336 Improved AbstractPlugin verifying of Parameters using isEmpty() method.
- TAS-336 Added method AbstractPlugin.postProcessParameters() to allow extending Plugins to modify / verify Parameters after user input and before processData().
- TAS-356 Added code to AbstractPlugin that only prints out any given citation once.
- TAS-203 Fixed imports that included Java 8
- TAS-203 Adding tag direction reads and taxaList support to TBT module
- TAS-336 Added isEmpty() method to PluginParameter
- TAS-382 Changed labels on Manhattan and QQ Plots to say Log10 instead of just Log.
- TAS-203 Beginning to implement HDF5 TBT following the builders and allowing for greater depth
- Implement phenotype package classes and add PhenotypeBuilder
- Implement phenotype package classes, add new FilterPhenotype class, delete ImportPhenotypeUtils
- remove TagHDF5 interface
- new version of CorePhenotype
- make many changes to Phenotype implementation
- TAS-380 Converted MergeMultipleTagCountsPlugin to Plugin Parameter Design
- TAS-375 Minor modifications to setters/getters of BinaryToTextPlugin
- TAS-340 Removed unneeded comments in Position
- TAS-340 Removed unused import from WHICH_ALLELE
- TAS-340 Fusion of Position.Allele with WHICH_ALLELE
- TAS-340 Move Which_Allele upto DNA and add ancestral allele
- TAS-377 Added builder like methods to FastqToTagCountPlugin
- added AbstractPhenotype as the superclass of CorePhenotype and FilterPhenotype
- TAS-377 Converted FastqToTagCountPlugin to new Plugin Parameter Design
- change Phenotype to interface implemented by CorePhenotype and FilterPhenotype
- TAS-378 Implement some methods in CombineGenotypeTable
- TAS-336 Added validity check for IN_DIR and OUT_DIR to PluginParameter.  Also added File Choosers that allow on directories to be chosen.
- added TaxaAttribute and other implementation details
- TAS-336 Added IN_DIR and OUT_DIR as possible file types in PluginParameter
- TAS-336 Changed getter methods of BinaryToTextPlugin to return BinaryToTextPlugin instead of Plugin.  This is so casting not needed when stringing calls together.
- initial implementation of Phenotype objects
- TAS-375 Removed Enum from BinaryToTextPlugin and changed command line flags to original for compatibility.
- TAS-336 Moved all logic for checking Range inside PluginParameter.
- TAS-303 HapMap reader that throws error messages with the row of bad data
- TAS-160 Add ancestral, major, minor, and reference support to HDF5 position lists.
- TAS-336 Added methods to Plugin to make defining Enum of PluginParameters optional.
- Fixed PanAUsageExample.java
- TAS-346 - FILLIN GUI removed as it is being supplanted by new plugin system
- TAS-355 - Fixes pair site ordering incorrect HapMap must be ordered by position
- PanA code completed.
- TAS-132 : IBSDistanceMatrix is not using start & endword consistently - shifted to first and last word
- Tools to ignore depth of coverage when converting from VCF to HDF5
- TAS-336 Improved AbstractPlugin.getParameterInstance() to return null if not found.  This allows for a better error message.
- TAS-336 Added PluginParameter.Builder construction that takes String in addition to Enum
- TAS-368 Improvements to plugin that lists all plugins.
- TAS-374 Add UNEAK to TASSEL5 GBS - Code from Jason
- TAS-368 Initial attempt at plugin that lists all plugins.
- TAS-336 Changed AbstractPlugin.printParameterValues() to do nothing if no parameters.
- TAS-375 Convert BinaryToTextPlugin to new Plugin Parameter Design
- TAS-336 Migrating new Parameter Plugin design into the main code base.
- Pick up training data set, transform with boxcox, convert to arff. PanA Start M5 training
- TAS-367 : replaced System.out.println() calls with myLogger.info() or myLogger.error()

## (V5.0.6) May 15, 2014

- TAS-103 Added check to prevent exception in Production Pipeline in case usage statement is desired.
- TAS-103 Additional modifications to logging of Production Pipeline
- TAS-103 More adjustments to logging for the Production Pipeline.
- TAS-103 Corrected logging of Production Pipeline.
- TAS-103 Improved logging of Production Pipeline to send stdout, stderr, and log4j into same file.
- TAS-315 Final bug found in VCF to HDF5 conversion
- TAS-103 Added logging statement to ProductionPipeline to put arguments passed to ProductionSNPCallerPlugin into the log.
- TAS-103 Rewrote ProductionPipeline using new Plugin Parameter Design.  Code is much shorter and simpler.
- Adding alignments from TOPMV3 to TagGWASMap
- TAS-336 Added New BinaryToTextPlugin to GUI under new GBS Menu
- TAS-336 Make printParameterValues() protected so that plugins that extend can call that.
- TAS-103 Continue simplification of Production Pipeline
- made rawSeqFileNameRegex public static final and removed rawSeqFileNameReplaceRegex
- Modified clearVariants() and clearVariant() so that they work with ragged arrays
- Fixed getVariantOff() so that it works with ragged arrays (additional change)
- TAS-363 Added synchronized statements to addTaxa() and addSite() methods to prevent race condition.
- Fixed getVariantOff() so that it works with ragged arrays
- Added ".topm" as an acceptable suffix for a binary TOPM in writeTOPM()
- 1. Reformat tag GWAS mapping result to HDF5 2. Align tags in mapping result using Bowtie2 --very-sensitive-local -k 2 option 3. Add alignments into TOPMV3 (Multiple alignment topm), which is easier to identify unique ref tags for model training
- Plugins of tag genetic mapping
- TAS-362 Fix in GenotypeTableBuilder that sorts genotypes when it sorts the positions.
- TAS-336 Added class BinaryToTextPluginNew which I migrated BinaryToTextPlugin to use the new Parameter Plugin design.
- TAS-336 Improved Tassel Pipeline it accept either of these Plugin constructors Plugin(Frame) or Plugin(Frame, boolean).  That way if Plugin used in GUI and CLI, only Plugin(Frame, boolean) needs to be implemented.
- TAS-336 Plugin Parameters Design.  Added code to create JComboBox if parameter is type Enum.  Also added code to handle parameter default value vs required.
- TAS-336 Added code to Plugin Parameter Concept to check that directory of output files exist.
- modified -extOut in FILLINFindHaplotypes to include chromosome number and start/end positions
- added a feature for FILLINFindHaplotypes to output which taxa go into which haplotype (using flag -extOut)
- Genetic mapping code log update
- TAS-336 Added code for better validating plugin parameters.
- TAS-336 Continued development of auto-generating GUI for Plugin Parameter Concept.
- TAS-336 Added code to Parameter Plugin Concept to handle clicking "Cancel"
- TAS-336 Added code to Plugin Parameter Concept that automatically generates GUI for entering parameter values.
- TAS-315 Create memory efficient VCF to HDF5 conversion nearly functional
- Get block position from TOPM SimpleGenotypeSBit deflation set to 1
- TAS-315 Create memory efficient VCF to HDF5 insert part 2
- TAS-355 Added code to PositionListBuilder.validateOrdering() to print out two positions that are out of order.
- TAS-354 Changed Tassel Pipeline flag from -mergeAlignments to -mergeGenotypeTables
- TAS-354 Changed Menu and tool tip name from "Merge Alignments" to "Merge Genotype Tables"
- TAS-354 Renamed MergeAlignmentsPlugin to MergeGenotypeTablesPlugin
- TAS-336 Continued implementation of PluginParameter concept.
- TAS-352 Add code to FileLoadPlugin to invoke reading FASTA files.
- TAS-352 Added support for loading FASTA files
- TAS-351 Corrected problem loading Phylip files
- TAS-336 Proposed solution to one step to get parameter value for Plugins.
- TAS-336 Continue implementing ParameterConceptPlugin
- TAS-349 Fixed operation of Cancel button on Kinship Dialog
- TAS-336 In ParameterConceptPlugin changed getParameterValue() to getParameter().  This will go in AbstractPlugin once this design is finished.  Added main() to illustrate how Plugin can be setup with a Builder style method calls.
- added support for anonymous haplotype naming in FILLINFindHaplotypesPlugin
- added support for anonymous haplotype naming in FILLINFindHaplotypesPlugin
- TAS-336 Implementing getParameterValue() and setParameterValue() for PluginParameter
- TAS-326 and TAS-345. Also added an error IOException for when the output file (-o) in FILLINFindHaplotypesPlugin does not include gX
- TAS-326 and TAS-345. Also added an error IOException for when the output file (-o) in FILLINFindHaplotypesPlugin does not include gX
- TAS-336 Added Range checks for PluginParameter
- TAS-336 Continued implementation of PluginParameter concept.
- TAS-342  Added code to retain current view of TableReports which clicking between other data tree nodes.
- TAS-341 Migrate MergeAlignmentsPlugin to Tassel 5.  Written by Jason
- TAS-334 Added ability to Data Tree to show text and graphical view of TOPM
- TAS-169 Removed unused code from GenotypeSummaryPlugin
- TAS-315 Minor progress towards create memory efficient VCF to HDF5
- TAS-334 Initial implementation of Graphical TOPM Viewer
- TAS-333 : Improvements to GeneralPosition Builder
- Minor changes to the genomeVersion reporting when writePositions is true
- TAS-315 Beginning of a more efficient VCF to H5 converter.
- TAS-333 : Lots of changes so that adding the reference allele to an HDF5 genotype table will work "right"
- I ran existing unit tests where possible and they all passed. If I messed something up, we can always revert. This is an initial stab at this which could be refined in many ways. For example, population of a PositionList with reference alleles should probably be moved to a PositionListUtils class.
- moved methods in FILLINImputationAccuracy to a subclass
- TAS-336 Plugin parameter concepts highlighted for everyone
- TAS-335 : Simple plugin to split an HDF5 GenotypeTable by chromosome
- TAS-333 : Minor change: System.out.println()'s of the updated positions is now optional (hard coded to no output)
- These were for debugging/checking
- TAS-322 : changed GeneralPosition to Position when comparing PositionLists (better style?)
- TAS-333 : Now populates the global minor allele
- This new plugin is nearly ready
- TAS-324 Removed enum type FileLoadPlugin.TasselFileType.Numerical.  FileLoadPlugin.TasselFileType.Phenotype is the new and preference file format.
- TAS-322 : Added dash missing from line myArgsEngine.add("-iD", "--ignore-depth", false);
- LDResult class is expanded to permit passing around LDResults for other analyses
- Should provide a robust mechanism for only retaining the most significant results.
- TAS-322 : Added an ignoreDepth (-iD) option to the MergeHDF5GenotypesSameSitesPlugin, as well as more efficient reading of positions to confirm that they match.
- TAS-322 : Plugin to merge multiple HDF5 GenotypeTables containing the exact same sites
- TAS-323 - Added public static FilterGenotypeTable getInstance(GenotypeTable a, PositionList subPositionList)
- Added " sec" to all the nanosecond timer System.outs
- Fixed bug which slowed down writing depths to HDF5 GenotypeTable
- TAS-316: LibraryPrepID's can now be alphanumeric (i.e., Strings) OR integers. Stricter checking that the order of the columns is correct.
- Added  GenotypeTableBuilder.getTaxaIncrementalWithMerging(myTargetHDF5file, myPositionList, genoMergeRule)
- So that a new hdf5 geno file can be created to which replicate taxa (from different flowcell_lanes in the same ProductionSNPCallerPlugin run) can be added & merged

## (V5.0.5) April 17, 2014

- Added a function writeBinaryFileForChromosomalRegion()
- Moved setUpGenotypeTableBuilder() and genos.closeUnfinished() or genos.build() outside of the for loop for each fastq file
- Added nano timing code for ProductionSNPCallerPlugin.readRawSequencesAndRecordDepth() Task #TAS-103
- Add Nano timing to each method called in the performFunction() to narrow down biggest performance component.  Task #TAS-103
- Added a noDepthOutput (-ndo) option
- Fixed problem with FIle Chooser not displaying sometimes in Java 7.  Task #TAS-312
- Added BbvCI-MspI
- added method to created kinship from numeric markers so that imputed data can be used
- added option to change all selected traits to a different type to so that numeric transformed data can be used as input to GS
- Fixed bug for in memory GenotypeTableBuilder merging of taxa
- Method for directly comparing GenotypeTables - eventually this functionality needs to be moved into the core accuracy
- TAS-286 - Fix two bugs in the Viterbi section of FILLIN code
- Added 5 new enzyme pairs
- implements new kinship method in the kinship plugin
- Fixed Display of TOPM variants in Tassel GUI Task #TAS-309
- Logging corrections to GBS Production Pipeline Task: #TAS-103
- Minor improvements to GBS Production Pipeline Task: #TAS-103
- Corrected Example Run Config File output in GBS Production Pipeline.  Task #TAS-103
- modified method to improve efficiency
- Added Utils.getNumberLinesNotHashOrBlank() to count lines in a file while ignoring lines that begin with # or are blank.  Initially to determine number of sites in a VCF file.

## (V5.0.4) April 10, 2014

- added new method for calculating kinship based on Endelman and Jannink (2012) G3.
- Added code to show number of tags in info. panel when viewing TOPM.  Task #TAS-309
- Added code to put GenotypeSummaryPlugin results into Result -> Genotype Summary Task #TAS-310
- Changed Tassel Data Tree to show only nodes that have data or results.  Task #TAS-310
- Removing SBit TBit icons from Tassel 5.  Only made sense in Tassel 4.  Task #TAS-311
- Added ability to import TOPMs into Tassel GUI.  Task #TAS-309
- Fixed this by setting 0 to Double.MIN_VALUE fixed the issue in genetic mapping code that likelihood ratio equals NaN
- TAS-305 fix
- TAS-212 : ProductionSNPCaller now annotates taxa in HDF5 GenotypeTable with Flowcell_Lanes
- Something went wrong with the first commit (not quite the latest code changes), so I am trying again
- TAS-212 : ProductionSNPCaller now annotates taxa in HDF5 GenotypeTable with Flowcell_Lanes
- TAS-308 - Bug Fix VCF Reading Error with Insertions and Missing Reference Data
- TAS-305 fix
- TAS-286 - Clean up FILLIN imputation code.  Substantial renaming and reorganizing the FILLIN code that improves method names.  Second Push - some issue with first push
- TAS-286 - Clean up FILLIN imputation code.  Substantial renaming and reorganizing the FILLIN code that improves method names.
- Added icon for menu item Data -> Homozygous Genotype Task #TAS-273
- Added icon for menu item Help -> Help Manual Task #TAS-273
- Using SuperByteMatrix instead of byte[][] in BuilderFromHapMap process threads to increase performance.  Task #TAS-304
- Updated pipeline package-info Task #TAS-298
- Updated package-info.java for net.maizegenetics.pipeline Task #TAS-298
- Updated package-info.java for net.maizegenetics.prefs Task #TAS-298
- TAS-282 - Clarified the difference between packages trait (past) and phenotype (future)
- TAS-282 - Remaining package-info level stubs.
- TAS-282 - Stub out package-info for most matrix algebra and TASSEL packages
- TAS-286 - Explicitly use maximumInbredError rate for Viterbi section of FILLIN code
- TAS-301:  Bug Fix: HDF5 GenotypeTable read always creates AlleleDepth even if there is no depth
- fixed the ratios for the focus block thresholds so that they update the global variables from the plugin

## (V5.0.3) April 3, 2014

- Update Home.html Task #TAS-296
- Added code to TasselPipeline.addForFlagsIfNeeded() to not add -fork flags if there are no flags to begin with.  Task #TAS-289
- Continued refactoring of GBS ProductionPipeline Task #TAS-103
- Fixed this by setting 0 to Double.MIN_VALUE fixed the issue in genetic mapping code that likelihood ratio equals NaN
- TAS-286 - Clean up FILLIN imputation code.  Bug in Guava code with Java 7 versus 8 Ordering fixed.
- TAS-286 - Clean up FILLIN imputation code.  Moved key Virterbi methods over to FILLINImputationUtils, added unit tests of Utils, sorted the error rates with a 0.5 fudge factor so that number of sites plays a important role, developed more robust search method for best hybrids
- Fixed LinkageDiseqDisplayPlugin menu icon.  Task #TAS-273
- Fixed LinkageDisequilibriumPlugin icon.  Task #TAS-273
- Added icon for "About" menu item.  Task #TAS-273
- Fixed "Show Memory" menu icon.  Task #TAS-273
- Change chromosome to int type for tag genetic mapping result
- Continued refactoring of GBS ProductionPipeline Task #TAS-103
- Simplifying flag to run imputation on ProductionPipeline.  Task #TAS-103
- Task #TAS-220 Create public static GenotypeCallTableBuilder getInstanceCopy(GenotypeCallTable genotype) {
- Added convenience methods getInstance() and getInstanceRemoveSiteNames() to FilterGenotypeTable that accept a List.  Task #TAS-217
- Task #TAS-289 Add handling of one-fork command lines without need of -fork flags
- Small documentation fixes
- Beginning to clean up FILLIN code, part 1
- Added date to javadoc title
- Move Fei's code to appropriate analysis.gbs package
- TAS-282 - Stub out package-info for all packages.  Completion of analysis.
- Added Hets to Homozygous Plugin function to the command line.  Task #TAS-197 - Add "Hets to Unknown" (GenotypeTableBuilder.getHomozygousInstance()) to GUI
- TAS-209 - GenotypeTable.startAndEndOfChromosome() converted method name to firstLastSiteOfChromosome() to agree with inclusiveness
- Adding HetsToUnknownPlugin and to Data Menu Task #TAS-197 - Add "Hets to Unknown" (GenotypeTableBuilder.getHomozygousInstance()) to GUI
- TAS-181 : Found bug in FILLINImputation projection output code.  Added an exclusion of heterozygous sites from donors for Viterbi informative sites.  Highlighted an error with AccuracyPlugin exposing variable in ImputationPlugin.
- TAS-181:  A separate concept of a minimum distance between donor haplotypes for Viterbi
- Updated startAndEndOfChromosome() method to have first and last site returned to be inclusive.  Task #TAS-209 - GenotypeTable.startAndEndOfChromosome() returns int[]{start,end} INCLUSIVE but documentations says that end is exclusive
- Added Tassel Preferences to remember size of Tassel GUI at exit and next start of Tassel returns to that size.
- Added private support for calculating accuracy within MAF categories
- Added a FILLINImputationAccuracyPlugin that masks files, calls FILLINImputationPlugin, and calculates R2 accuracy. Depth and MAF not currently supported.  Removed the focus block error threshold user inputs and fixed them to MxInbredErr and MxHybErr
- Added restriction enzymes CviQI and Csp6I (isoschizomers of each other)
- change default PCA type to covariance
- TAS-181 : Found bug in FILLINImputation projection output code.  Added an exclusion of heterozygous sites from donors for Viterbi informative sites.
- parseSAMAlignment() now properly handles tags that don't align but are reverse complemented (i.e., flag = 20, or with 0x10 and 0x4 bits set).
- These now get reverse complemented before they are stored in the TOPM.
- Additional changes needed so it works with ragged arrays for variantOffsets[] and variantDefs[]
- Added new Tassel Logo to About Box dialog.
- TAS-181 : Cleaning up FILLINImputation code.  Removing System out and adding comments to code to make paragraphs more understandable.
- Added ability to use scientific E notation in selecting position ranges
- TAS-196 - "Site Filter" Problem with End Physical Position when Selecting chromosomes is resolved.  Issue with single chromosome GenotypeTable being accesses as the null chromosome.
- TAS-181 : Rename ProjectionAlignmentIO to ProjectionGenotypeIO to follow new convention.  Improved the JavaDoc.
- Remove obsolete code from build.xml
- modified to use DoubleMatrix instead of Colt DoubleMatrix2D and PrinComp instead of PCA in order to eliminate dependence on Colt
- fixed bug in rowSums
- add principal component class that uses DoubleMatrix
- TAS-181 - Fixed bug on on removing accuracy variables.
- TAS-181 - Improved documentation
- TAS-181:  Improving documentation of FILLIN code
- Fixed one bug on reverse direction Viterbi imputation
- Fixed a bug when chromosomes don't increase successively starting from 0 (e.g., 9 & 10 rather than 0,1,2,3â¦)
- Added getChromosomeIndex(int intChrName)
- updates to StepwiseOLSModelFitter to implement scanning
- Gets basic projection output working from FILLIN.  Still too many breakpoints with real data.
- Add a -ko (keepOpen) option to decide whether to call genos.build() or genos.closeUnfinished() once all the fastq files are read.
- Fixed a bug related to positionToSite
- Now uses a BasicGenotypeMergeRule to call genotypes for new taxa not already present in the hdf5 genotypes file
- Crediting Jeff for his code in BasicGenotypeMergeRule
- Initial reconstitution of projection genotype support to FILLIN
- Has a bug for heterozygous start regions
- Added comments explaining likelihoodRatioThreshAlleleCnt[].
- These will be deleted from the ProductionSNPCallerPlugin, so I didn't want to lose them.
- Support for empty position list
- Fixed a few strange characters in javadoc comments
- Improvements to Production Pipeline.
- Most of the changes required to make work this work with depth, and adding taxa from one fastq at a time.  Not tested yet!
- Moved accuracy out of FILLINImputationPlugin() to an independant class, ImputationAccuracy()
- shortened -errRate to -eR for consistency with ProductionSNPCallerPlugin
- Move Tassel5HDF5 to TagsHDF5
- HDF5 tags design for model training after genetic mapping tags
- Added support for target files that do not match donor files for all snp positions
- Bug fixes to FILLIN
- Fixed a bug in FILLINImputationPlugin() where all genotypes imputed as a het in two combination haplotype mode are set to missing, not only taxa below the mxHet threshold (inbreds). Also, fixed a bug in setParameters() where if mxHet or propSitesMasked not flagged then threw an error because attempting to parse a null value.
- Added a new addTaxa method that accepts int[][] depths
- Fixed typo in java doc for depthIntToByte(int[][] depth)
- Added the function depthIntToByte(int[][] depth) to handle 2D arrays
- A concept for plugin parameter that combines bounds checking, descriptions, GUIs all together
- very minimal tweaks fillin accuracy output
- slight tweaks fillin
- slight tweaks fillin
- added accuracy to FILLIN
- Added FilterSubsetPlugin created by Jason and flags for running to TasselPipeline
- Added jfxrt.jar (from jdk1.7.0_45.jdk) to lib directory for compiling javafx code.
- Additional JavaFX concepts for plugins
- JavaFX stub that works with our plugins and Swing
- Modified ProductionPipeline to send emails with subject text, that works with Basecamp.  These will be archived in the Cassava Production Pipeline Log project.
- GBSV3 alignment hypothesis test by genetic mapping
- Cleaning up of FILLIN and making non Verbose option
- Create low level HDF5 copy and subset - TAS-171
- Method for converting TableReport to Guava Table
- Add unlockHDF5TaxaModule and unlockHDF5GenotypeModule
- Fixes TAS-162 bug on single position for last chromosome
- Completed implementation of HDF5 depth export
- HDF5 with depth and merge working
- Genetic mapping code. TagAgainstAnchorPlugin

## (V5.0.2) February 20, 2014

- Fixed bug in writing depth to HDF5.  Still have bug on appending to HDF5.
- Changed javadoc of AlleleDepthBuilder.setDepth() to match new behavior
- In memory implementation of merging and recalling genotypes with depth in GenotypeTableBuilder
- Corrected -diversitySlidingWin flag implementation.
- added BasicShuffler class from Tassel ver 3 for permuting primitive arrays
- added class from Tassel v3 for permuting within levels of a factor
- Added GenotypeBuilder & HapMap sort methods
- Add chromosome to site summary
- Methods for merging taxa genotype calls
- Methods and classes to start allowing taxa merging.
- Minor changes to logging output regarding irrelevant fast files (only those in key file are reported)
- fix bugs in callParentAlleles methods introduced by the version 5 port
- fix bug in R2 computation
- Pushing HDF5 functions into HDF5Utils
- Important!  Shift HDF5 from TASSEL 4 format to TASSEL 5 format.
- Taxa annotation reading from HDF5 file.  Final commit before break from old T4 HDF5.
- Renamed TaxaList builder signatures when building just from Genotypes
- Methods for saving taxa annotations in VCF and HapMap files and for roundtripping
- cleaned up jglim and jglim.dm packages and moved all the active classes to stats.linearmodels. In addition, moved a couple of classes to analysis.association package that referenced packages outside stats and matrixalgebra. Stats should now only reference the matrixalgebra packages.
- Move analysis DistanceMatrices to analysis package.  Slightly changed the interface and renamed the interface to TASSEL 5 standards.
- Few new TaxaList annotation filter methods
- changed net.maizegenetics.gwas.modelfitter to net.maizegenetics.analysis.modelfitter
- moved AGPMap from gwas.NAM package. Then deleted the gwas.NAM  package.
- Move pd package to dna package
- Move distance and tree to taxa package
- Removed old gbs packages
- Moved BaseEncoder to net.maizegenetics.dna
- Moved net.maizegenetics.gbs.util to net.maizegenetics.dna.tag
- Moved net.maizegenetics.gbs.homology to net.maizegenetics.analysis.gbs
- Moved net.maizegenetics.gbs.pipeline to net.maizegenetics.analysis.gbs
- Migrate popgen to analysis.popgen
- Updates to ProductionPipeline for Cassava
- More updates to ProductionPipeline for Cassava Project
- Initial implementation of annotated taxa files
- Current reader will work with either a custom file or GBS keyfile.
- allow population data to be set externally instead of read from a pedigree file to enable testing with simulated data sets
- move net.maizegenetics.gwas.imputation classes to net.maizegenetics.analysis.imputation and net.maizegenetics.gwas.imputation.clustering to net.maizegenetics.gwas.analysis.clustering
- Fixes a bug in TASSEL 4 HDF5 that did not set taxa number correctly
- This is a hack to upgrade the file and fix the bug on the fly.  Once the new HDF5 directory structure is set this needs to be migrated to migration code.
- Initial Production Pipeline for Cassava Project
- Added char complementation for NucleotideAlignmentConstants.  Removed complementation from DiscoveryPipeline
- Fixing a commit and defining a char version of IUPAC lookup
- Char version of IUPAC look up
- Simplified NucleotideAlignmentConstants.getNucleotideIUPAC() implementation
- Moved SequenceDiversityPlugin to net.maizegenetics.analysis.popgen
- Moved ProjectionLoadPlugin to net.maizegenetics.analysis.data
- Moved ExtractHapmapSubsetPlugin to net.maizegenetics.analysis.data
- Moved GenotypeImputationPlugin to net.maizegenetics.analysis.imputation
- Moved KinshipPlugin, DistanceMatrixPlugin, and DistanceMatrixRangesPlugin to net.maizegenetics.analysis.distance
- Moved display plugins to net.maizegenetics.analysis.chart
- Moved net.maizegenetics.baseplugins.numericaltransform to net.maizegenetics.analysis.numericaltransform
- Move net.maizegenetics.baseplugins.chart to net.maizegenetics.analysis.chart
- Removed GenotypeImputationPlugin from GUI and CLI.
- Moved RidgeRegressionEmmaPlugin to net.maizegenetics.analysis.association
- Move ConvertAlignmentCoordinatesPlugin to net.maizegenetics.analysis.data
- Fixing icon locations after moving plugins.
- Moved MLM and GLM plugins to net.maizegenetics.analysis.association
- Migrating PE sequencing code and TOPMV3 code (annotating TOPM with multiple aligners and PE)
- Moved PlinkLoadPlugin to net.maizegenetics.analysis.data
- Moved LD Plugins to net.maizegenetics.analysis.popgen
- Moved CombineDataSetsPlugin.java, ExportMultiplePlugin.java, ExportPlugin.java, FileLoadPlugin.java, GenotypeSummaryPlugin.java, IntersectionAlignmentPlugin.java, MergeAlignmentsPlugin.java, MergeAlignmentsSameSitesPlugin.java, PassThroughPlugin.java, SeparatePlugin.java, SynonymizerPlugin.java, UnionAlignmentPlugin.java, ExtractHapmapSubsetPlugin.java, PlinkLoadPlugin.java, and ProjectionLoadPlugin.java to net.maizegenetics.analysis.data
- updated with changes from v4, deleted unused methods
- Initial migration of FILLIN
- Moved filter plugins to net.maizegenetics.analysis.filter
- Move CreateTreePlugin and ArchaeopteryxPlugin to net.maizegenetics.analysis.tree
- Deletion of some complementation code in DiscoveryPipeline.  Still need to replace more.
- Update javadoc so that it works with current dna.snp packages
- updated with version 4 changes. NucleotideImputationUtils still has many
- bugs/type conversions to fix yet.
- Creation of analysis.imputation package
- Moved DonorHaplotypes To dna.map
- Implementation of VCF
- High depth scores are calculated on fly
- Sort Taxon annotations by key
- Initialing NUCLEOTIDE_ALLELE_ARRAY to UNDEFINED_ALLELE instead of UNKNOWN
- Progress on VCF export
- Expansion of snp depth interface so that it can be accessed an used more easily
- implemented code for filtering sites by subpopulation
- slight change to transition probability so that segment length always returns a positive value
- correct problems with code so that it actually works
- changed getUsage so that it displays as info instead of error using ?
- add new clustering methods
- Fixed migration bug with VCF export and unknown diploid allele versus allele
- VCF export with sample annotations
- GeneralAnnotation with method to get all annotations as a Map
- Corrected MemoryAlleleDepth constructor to dynamically get max number of alleles.
- VCF reading support with depth, taxa and position annotations
- GeneralPosition Builder method for key-value pairs common in VCF
- Bug in AlleleDepth allele number
- Taxon with annotations to string
- Created fast array version of getNucleotideAlleleByte
- Fixed race condition with CombineDataSetsPlugin that could cause following plugins to run more than once.
- Improved reporting errors when setting parameters for Self-describing plugins.
- Corrected VCF File entry when Tassel uses major allele instead of reference.
- Removed extra ## from VCF output
- TODO added to highlight endWord feature not implemented correctly.
- Update ImputationUtils for readFromHapMap
- Improve documentation of GenotypeTable
- Simplified default import/export of HapMap files
- Begin classes for centralized read/write access to HDF5 files.
- [tassel:bugs] #172 Command-line filterAlignLocus bug: set the start and stop sites if they haven't been
- Improved reading of labels on import dialog
- Removed import / export options for BitAlignmentHDF5.  That was for Tassel 4.  Now HDF5 format stores byte values.
- Changed text on import dialog to say "Load Trait (data, covariates, or factors)"
- Removed Annotated file format from import options
- Removed "Polymorphism" format from import options.
- Removed Flapjack import / export options.

## (V5.0.1) January 16, 2014

- Moved GenotypeTableMask classes from net.maizegenetics.dna.snp to net.maizegenetics.gui
- Added logging statement to Tassel GUI and Pipeline to report current Java version being used.
- Corrected get indexOf bug in ProjectionIO
- Corrected bug with Sequence Viewer when using Genetic Distance to color sequences with heat map.
- Convenience method to projection building
- return position of marker from GenotypeTable instead of NaN
- adds Stepwise method to analysis menu
- add an icon to the GUI menu
- model fitter: added model fitting criteria (bic, etc) to the user dialog. Added an icon to the plugin. Changed the name to Stepwise in the dialog.
- class for importing phenotypes, implementation in progress
- Example use of static import to reduce boilerplate on Constants and utils classes
- Fixed bug in FilterGenotypeTable.getInstance() introduced when converted indexMatchingOf() to indexOf().  Found when running basic test suite.
- Fixed WHICH_ALLELE parameters in imputation package.
- Implemented WHICH_ALLELE design instead of SCOPE
- Added index to WHICH_ALLELE enum
- Added code to AbstractPlugin to output plugin class and citation when 100% progress reported.
- Added hdf5 output and various minor tweaks
- Added Plugin.getCitation() method that returns citation for plugins (defaults to Tassel citation).
- Added new int to byte decoding algorithm for Allele Depths designed by Robert B.
- Beginning implementation of AlleleDepth reading HDF5
- Changed all AlleleDepthBuilder methods to take depths in this orderâ¦  First dimension of depths is number of alleles (6 for Nucleotide) and second dimension is sites.
- Continued implementation of depth() methods
- Converted indexMatchingOf to the standard list indexOf
- TaxaList shift to no duplicates
- Correction to javadoc
- Initial implementation of AlleleDepthBuilder.addTaxon() for HDF5.
- Initialize Depth Module in HDF5 file.
- Implemented AlleleDepthUtil.depthByteToInt() and depthIntToByte() based on Robert's formula.
- Implemented AlleleDepthUtil.depthIntToByte().  Still need the other direction.
- Package Info for Utility class
- Elimination of Todo of TaxaArrayList
- Trivial documentation update
- ProductionPlugin updated for efficiency and TASSEL5
- Added stubs for AlleleDepth.depthForTaxon() and depthForSite()
- Initial implementation of HDF5AlleleDepth
- Removed get from method names in AlleleDepth
- Better design for AlleleDepthBuilder constructors.
- Initial AlleleDepthBuilder.addTaxon() methods added.
- Changed AlleleDepth.getDepthForAllele() to return int instead of byte.  Also added new getDepth() method.  This is to translate the byte storage to real integer value.  This had minor changes to GenotypeTable interface, VCF Export, and MergeIdenticalTaxaPlugin
- Implemented AlleleDepthBuilder.setDepth() methods.
- Corrected } syntax in AlleleDepthBuilder
- Continued structuring AlleleDepthBuilder to support HDF5 and incremental taxa.
- Structuring AlleleDepthBuilder to support HDF5 and incremental adding of taxa.
- Initial utility for translating between byte and int depth values
- Converted some syntax to Java 7
- Weigher for new BitStorage
- New implementation of BitStorage that only has one byte[] array for comparing and creating one BitSet.  This is proposal to replace current implementation.
- Adding genotypeAllTaxa(site) to GenotypeTable interface
- Removing old ProjectionAlignment
- Removed incorrectly implemented DynamicBitStorage.getInstance() method.
- Cosmetic change to logging output (fewer blank lines)
- Renamed getDepthForAlleles(taxon, site, allele) to getDepthForAllele
- Continue with initial implementation of AlleleDepth
- Made some variable names clearer (alleles -> allelesByFreq)
- changed alignment to GenotypeTable
- adds projection alignment import
- Additional documentation of dna.snp
- Class level javadoc for dna.snp package
- Removed unneeded logging statements.
- Fixed problem with DistanceMatrix copy constructor.  Thanks Peter for finding.
- Renamed PositionList.siteCount() to numberOfSites()
- Added check to GenotypeTable.getInstance() methods to ensure genotypes number of sites and taxa match site list and taxa list.
- Fixed a bug with the PositionListBuilder in getVariableSites():  An extra position was added whenever the aligned reference tag contained a gap.
- change pal.phenotype to phenotype
- adding new classes for phenotype package
- Continued documentation of dan packages
- Starting to package level documentation of the DNA package
- Removal of String taxa access method  in DistanceMatrix
- Renamed TagsAtLocus to TagLocus
- Complete rename of Alignment to GenotypeTable
- Begin rename Alignment to GenotypeTable
- Completed renaming Genotype to GenotypeCallTable
- Begin renaming Genotype to GenotypeCallTable
- Renamed net.maizegenetics.dna.snp.genotype to net.maizegenetics.dna.snp.genotypecall
- Cleaned up the DiscoverySNPCaller prior to jUnit testing
- Update jar Signatures
- Updated run_pipeline.pl and start_tassel.pl to work with older versions of perl

## (V5.0.0) (Initial Build) December 5, 2013


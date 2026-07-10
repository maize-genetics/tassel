# **TASSEL 3.0 / 4.0 Pipeline Command Line Interface:** _Guide to using Tassel Pipeline_

## **Terry Casstevens (** _tmc46@cornell.edu_ **)**

Institute for Genomic Diversity, Cornell University, Ithaca, NY 14853-2703

March 28, 2014

|Prerequisites ............................................................................................................................................................ 1 <br>|
|---|
|Source Code ............................................................................................................................................................ 1 <br>|
|Install ....................................................................................................................................................................... 1 <br>|
|Execute .................................................................................................................................................................... 1 <br>|
|Increasing Heap Size ............................................................................................................................................... 2 <br>|
|Examples ................................................................................................................................................................. 2 <br>|
|Examples_(XML Configuration Files)_.................................................................................................................... 2 <br>|
|Usage ....................................................................................................................................................................... 3 <br>|
|Pipeline Controls ................................................................................................................................................. 3 <br>|
|Data ..................................................................................................................................................................... 3 <br>|
|Analysis............................................................................................................................................................... 8 <br>|
|Results ............................................................................................................................................................... 10|

## Prerequisites

- Java SDK 6.0 or later (http://java.sun.com/javase/downloads/index.jsp).

- Tassel Standalone Build (http://www.maizegenetics.net/tassel/tassel3.0_standalone.zip)

- Tassel Standalone Build (http://www.maizegenetics.net/tassel/tassel4.0_standalone.zip)

## Source Code

- https://tassel.svn.sourceforge.net/svnroot/tassel/maizegenetics

## Install

Unzip the Tassel Standalone Build onto your file system.  Change into the root directory: . `tassel3.0_standalone or tassel4.0_standalone`

## Execute

On Windows, use `run_pipeline.bat` to execute the pipeline.

In UNIX, use `run_pipeline.pl` to execute the pipeline.  If you are using a Bash Shell on Windows, you may need to change the following line to use a ; instead of a :.

```
my $CP = join(":", @fl);
```

1

To launch the Tassel GUI that automatically executes a pipeline, use `start_tassel.bat` or `start_tassel.pl` instead of `run_pipeline.bat` or `run_pipeline.pl` respectively.

These scripts have a $top variable that can be changed to the absolute path of your installation.  That way, you can execute them any directory.

## Increasing Heap Size

To modify the initial or maximum heap size available to the Tassel Pipeline, either edit `run_pipeline.pl` or specify values via the command line.

```
./run_pipeline.pl -Xms512m -Xmx10g -fork1…
```

## Examples

```
./run_pipeline.pl -fork1 -h chr1_5000sites.txt -ld -ldd png -o
chr1_5000sites_ld.png -runfork1
```

```
./run_pipeline.pl -fork1 -h chr1_5000sites.txt -ld -ldd png -o
chr1_5000sites_ld.png -runfork1
./run_pipeline.pl -fork1 … -fork2 … -combine3 -input1 -input2 … -fork4 -
<flag> -input3 -runfork1 -runfork2
```

## **Examples** **_(XML Configuration Files)_**

This command runs the Tassel Pipeline according to the specified configuration file...  Configuration files are standard XML notation.  The tags are the same as the below documented flags although no beginning dash is used.  See the `example_pipelines` directory for some common XML configurations.

```
./run_pipeline.pl -configFile config.xml
```

This command creates the XML configuration file from the original command line flags.  Simply insert the - createXML and filename at the beginning. Only the XML is created.  It does not run the pipeline...

```
./run_pipeline.pl -createXML config.xml -fork1 ...
```

<mark>This command translates the specified XML configuration file back into the original command line flags... It does not run the pipeline...</mark>

```
./run_pipeline.pl -translateXML config.xml
```

2

## Usage

|**Pipeline Controls**||
|---|---|
|`-fork<id>`|`This flag identifies the start of a pipeline`<br>`segment that should be executed sequentially.`<br>`<id> can be numbers or characters (no spaces).`<br>`No space between -fork and <id> either.  Other`<br>`flags can reference the <id>.`|
|`-runfork<id>`|`This flag identifies a pipeline segment to`<br>`execute. This will usually be the last argument.`<br>`This explicitly executes the identified pipeline`<br>`segment. This should not be used to execute`<br>`pipeline segments that receive input from other`<br>`pipeline segments. Those will start`<br>`automatically when it receives the input.`|
|`-input<id>`|`This specifies a pipeline segment as input to`<br>`the plugin prior this flag.  That plugin must be`<br>`in the current pipeline segment. Multiple of`<br>`these can be specified after plugins that accept`<br>`multiple inputs.`<br>`./run_pipeline.pl -fork1 -h genotype.hmp.txt -`<br>`fork2 -r phenotype.txt -combine3 -input1 -input2`<br>`-intersect -runfork1 -runfork2`<br>`./run_pipeline.pl -fork1 -h genotype.hmp.txt -`<br>`fork2 -includeTaxaInFile taxaList1.txt -input1 -`<br>`export file1 -fork3 -includeTaxaInFile`<br>`taxaList2.txt -input1 -export file2 -runfork1`|
|`-inputOnce<id>`|`This specifies a pipeline segment as a one-time`<br>`input to a CombineDataSetsPlugin. As such, this`<br>`flag should follow a CombineDataSetsPlugin in`<br>`the current pipeline segment. After the`<br>`CombineDataSetsPlugin has received data from`<br>`this input, it will use it for every iteration.`<br>`Whereas CombineDataSetsPlugin waits for data`<br>`specified by -input each iteration.  Multiple of`<br>`these can be specified.`|
|`-combine<id>`|`This flag starts a new pipeline segment with a`<br>`CombineDataSetsPlugin at the beginning. The`<br>`CombineDataSetsPlugin is used to combine data`<br>`sets from multiple pipeline segments. Follow`<br>`this flag with -input<id> and/or -inputOnce<id>`<br>`flags to specify which pipeline segments should`<br>`be combined.`|
|**Data**||
||`If the filename to be imported begins with`<br>`“http”, it will be treated as an URL.`|
|`-t <trait file>`|`Loads trait file as numerical data.`|

3

|`-s <SNP file>`|`Loads SNP file as sequence alignment.`|
|---|---|
|`-p <Poly file>`|`Loads polymorphism file as polymorphism`<br>`alignment.`|
|`-a <Anno file>`|`Loads annotated alignment file.`|
|`-r <phenotype file>`|`Loads file as Table Report. This is meant to`<br>`replace numerical, polymorphism, and annotated`<br>`alignment.`|
|`-k <kinship file>`|`Loads kinship file as square matrix.`|
|`-q <population`<br>`structure file>`|`Loads population structure file as numerical`<br>`data.`|
|`-h <hapmap file>`|`Loads hapmap file(.hmp.txt or .hmp.txt.gz)`|
|`-h5 <HDF5 file>`|`Loads HDF5 Alignment file(.hmp.h5). `|
|`-b <BLOB file>`|`Loads BLOB file(.zip). Only Tassel 3.`|
|`-g <BLOB file>`|`Loads BLOB file(.gz). Only Tassel 3.`|
|`-plink -ped <ped`<br>`filename> -map <map`<br>`filename>`|`Loads Plink format given ped and map files.`|
|`-flapjack -geno <geno`<br>`file> -map <map file>`|`Loads Flapjack format given geno and map files.`|
|`-fasta <filename>`|`Loads FASTA file.`|
|`-geneticMap <filename>`|`Loads Genetic Map. `|
|`-table`|`Loads a Table(i.e. exported from LD, MLM). `|
|`-vcf <filename>`|`Loads VCF file.`|
|`-importGuess <filename>`|`Uses Tassel Guess function to load file.`|
|`-maxAllelesToRetain`<br>`<number>`|`Sets the preference maximum number of alleles to`<br>`retain. Initial default is 2. Notice this is a`<br>`preference setting and remains true on a given`<br>`machine for all import operations that follow,`<br>`even in subsequent executions of the Tassel`<br>`Pipeline or Tassel GUI.`|
|`-retainRareAllelestrue`<br>`| false`|<br>`Sets the preference whether to retain rare`<br>`alleles. Notice this is a preference setting and`<br>`remains true on a given machine for all import`<br>`operations that follow, even in subsequent`<br>`executions of the Tassel Pipeline or Tassel GUI.`|
|`-optimizeForTaxa`|`This should follow a file import, which`<br>`instructs the loader to optimize the data for`<br>`taxa based operations.  By default, loaded data`<br>`is optimized for site operations.`|
|`-convertToSiteOpt`|`Converts input Alignment to a Site Optimized`<br>`Alignment.  May return same Alignment if already`<br>`optimized for Sites.  This is mainly for testing`<br>`purposes.`|
|`-convertToTaxaOpt`|`Converts input Alignment to a Taxa Optimized`<br>`Alignment. May return same Alignment if already`<br>`optimized for Taxa. This is mainly for testing`<br>`purposes.`|

4

|`-taxaJoinStrict <type>`|`If type Strict (deprecated value true), taxa`<br>`names are only considered the same if they match`<br>`exactly.  If type NonStrict (deprecated value`<br>`false), taxa names match if all specified levels`<br>`match.  For example… B73 matches B73:XXX matches`<br>`B73:XXX:YYY. But B73:XXX does not match B73:YYY.`<br>`If type NumLevels, then taxa names match if`<br>`specified number of levels (-taxaJoinNumLevels)`<br>`match.  For example… B73:XXX:YYY matches`<br>`B73:XXX:ZZZ for 1 or 2 number of levels.`|
|---|---|
|`-taxaJoinNumLevels`<br>`<num>`|`Specifies number of levels to use when`<br>`-taxaJoinStrict is NumLevels.`|
|`-union`|`This joins (union) input datasets based taxa.`<br>`This should follow a -combine specification.`|
|`-intersect`|`This joins (intersect) input datasets based`<br>`taxa. This should follow a -combine`<br>`specification.`|
|`-separate`<br>`<chromosomes…>`|`This separates an input into its components if`<br>`possible. For example, alignments separated by`<br>`chromosome (locus). For alignments, optionally`<br>`specify list of chromosomes (separated by commas`<br>`and no spaces) to separate.  Specifying nothing`<br>`returns all chromosomes.  Example:`<br>`run_pipeline.pl -fork1 -h file.hmp.txt -separate`<br>`3,6 -export -runfork1`|
|`-mergeAlignments`|`Merges multiple Alignments regardless of taxa or`<br>`site name overlap.  Undefined taxa / sites are`<br>`set to UNKNOWN.  Duplicate taxon / site set to`<br>`last Alignment processed.  Example:`<br>`run_pipeline.pl -fork1 -h file1.hmp.txt -fork2 -`<br>`h file2.hmp.txt -combine3 -input1 -input2 -`<br>`mergeAlignments -export files_merged.hmp.txt -`<br>`runfork1 -runfork2`|
|`-mergeAlignmentsSameSites`<br>`-input <files> -output`<br>`<filename>`|<br>`Merges Alignments assuming all sites are the`<br>`same in all Hapmap files. Input files separated`<br>`by commas without spaces.  The resulting file`<br>`may have incorrect major/minor alleles, strand,`<br>`center, etc.  It uses values from first`<br>`specified input file.  Checks that Site Name,`<br>`Chromosome, and Physical Position match for each`<br>`site.  Example: run_pipeline.pl -fork1 -`<br>`mergeAlignmentsSameSites -input`<br>`file1.hmp.txt,file2.hmp.txt -output temp -`<br>`runfork1`|
|`-export`<br>`<filename1,filename2,…>`|<br>`Exports input dataset to specified filename(s).`<br>`If no -exportType follows this parameter, the`<br>`exported format will be determined by the type`<br>`of input(i.e. Alignments will default to Hapmap`|

5

||`format). Exportable datasets, other that`<br>`Alignment, only have one format option.`<br>`Therefore, there is no need to specify -`<br>`exportType. Specify none, one, or multiple`<br>`filenames matching the number of input data`<br>`sets.  If no filenames, the files will be named`<br>`the same as the input data sets. If only one`<br>`specified for multiple data sets, a count`<br>`starting with 1 will be added to each resulting`<br>`file. If multiple filenames (separated with`<br>`commas but no spaces), there should be one for`<br>`each input. When exporting Hapmap files, if the`<br>`extension is .hmp.txt.gz, the file will be`<br>`gzipped.`|
|---|---|
|`-exportType <type>`|`Defines format that previously specified -export`<br>`should use. Type can be Hapmap, HapmapDiploid,`<br>`HDF5, VCF, Plink, zipBLOB (only Tassel 3),`<br>`gzipBLOB (only Tassel 3), Flapjack, Phylip_Seq,`<br>`Phylip_Inter, Text.`|
|`-impute`|`Imputes Genotypic Data.`|
|`-imputeMethod <method>`|`Specifies the impute method to use. Method can`<br>`be Length (default), MajorAllele, SimilarWindow,`<br>`or IBDProb. This should follow the -impute flag. `|
|`-imputeMinLength <num>`|`Specifies the minimum length for the impute`<br>`Length method (default value: 30). This should`<br>`follow the -impute flag. `|
|`-imputeMaxMismatch`<br>`<num>`|`Specifies the maximum mismatch for the impute`<br>`Length method (default value: 1). This should`<br>`follow the -impute flag. `|
|`-imputeMinProb <num>`|`Specifies the minimum probability for the impute`<br>`IBDProb method (default value: 0.001). This`<br>`should follow the -impute flag. `|
|`-filterAlign `|`Filters an alignment by sites.`|
|`-filterAlignMinCount`<br>`<num>`|`Specifies the minimum count (default: 1) for the`<br>`previously specified -filterAlign.`|
|`-filterAlignMinFreq`<br>`<num>`|`Specifies the minimum frequency (default: 0.0)`<br>`for thepreviously specified -filterAlign.`|
|`-filterAlignMaxFreq`|`Specifies the maximum frequency (default 1.0)`<br>`for thepreviously specified -filterAlign.`|
|`-filterAlignStart <num>`|`Specifies the starting site index (default`<br>`value: 0) for the previously specified -`<br>`filterAlign.`|
|`-filterAlignEnd <num>`|`Specifies the end site index (default value:`<br>`last site in alignment) for the previously`<br>`specified -filterAlign.`|
|`-filterAlignLocus`<br>`<name>`|`Specifies the Locus to be used with the starting`<br>`and ending physical positions if defined.`<br>`Defaults to first Locus in the Alignment.`|

6

|`-filterAlignStartPos`<br>`<num>`|`Specifies the starting physical position`<br>`(default is first site) for the previously`<br>`specified -filterAlign.`|
|---|---|
|`-filterAlignEndPos`<br>`<num>`|`Specifies the end physical position (default is`<br>`last site) for the previously specified -`<br>`filterAlign.`|
|`-filterAlignExtInd`|`Indicates that the last specified -filterAlign`<br>`should extract indels. This is not done by`<br>`default.`|
|`-filterAlignRemMinor`|`Indicates that the last specified -filterAlign`<br>`should remove minor SNP states. This is not done`<br>`by default.`|
|`-filterAlignSliding`|`Indicates that the last specified -filterAlign`<br>`should use sliding windows. This in not done by`<br>`default.`|
|`-filterAlignHapLen`<br>`<num>`|`Specifies the haplotype length (default value:`<br>`3) if using sliding windows.`|
|`-filterAlignStepLen`<br>`<num>`|`Specifies the step length (default value: 3) if`<br>`using sliding windows.`|
|`-includeTaxa`<br>`<taxon1,taxon2,…>`|`Filters input alignment to only include`<br>`specified taxa. The taxa should be separated`<br>`with commas and no spaces.`|
|`-includeTaxaInFile`<br>`<filename>`|`Filters input alignment to only include taxa`<br>`specified in file. The taxa cannot have spaces.`<br>`Individual taxa should be separated by white`<br>`space.`|
|`-excludeTaxa`<br>`<taxon1,taxon2,…>`|`Filters input alignment to exclude specified`<br>`taxa. The taxa should be separated with commas`<br>`and no spaces.`|
|`-excludeTaxaInFile`<br>`<filename>`|`Filters input alignment to exclude taxa`<br>`specified in file. The taxa cannot have spaces.`<br>`Individual taxa should be separated by white`<br>`space.`|
|`-includeSiteNames`<br>`<siteName1,siteName2,…>`|<br>`Filters input alignment to only include`<br>`specified site names. The site names should be`<br>`separated with commas and no spaces.`|
|`-includeSiteNamesInFile`<br>`<filename>`|<br>`Filters input alignment to only include site`<br>`names specified in file. The site names cannot`<br>`have spaces. Individual site names should be`<br>`separated by white space.`|
|`-excludeSiteNames`<br>`<taxon1,taxon2,…>`|`Filters input alignment to exclude specified`<br>`site names. The site names should be separated`<br>`with commas and no spaces.`|
|`-excludeSiteNamesInFile`<br>`<filename>`|<br>`Filters input alignment to exclude site names`<br>`specified in file. The site names cannot have`<br>`spaces. Individual site names should be`<br>`separated by white space.`|

7

|`-excludeLastTrait`|`This removes last column of Phenotype data. For`<br>`example… Can be used to remove last column of`<br>`population structure for use with MLM or GLM.`|
|---|---|
|`-subsetSites <num>`|`This filters an alignment to include a random`<br>`subset of sites. If <num> is >=1, it specifies`<br>`the total number of sites to keep. If it is a`<br>`decimal, it specifies the fraction of sites to`<br>`keep. Adding the flag "-step" immediately after`<br>`<num> tells the plugin to space the selected`<br>`sites evenly instead of randomly.`|
|`-subsetTaxa <num>`|`This filters an alignment to include a random`<br>`subset of taxa. If <num> is >=1, it specifies`<br>`the total number of taxa to keep. If it is a`<br>`decimal, it specifies the fraction of taxa to`<br>`keep. Adding flag "-step" immediately after`<br>`<num> tells the plugin to space the selected`<br>`taxa evenly instead of randomly.`|
|`-step`|`This tells the previously specified -subsetTaxa`<br>`or -subsetSites plugin to select sites/taxa`<br>`evenly across the alignment instead of randomly.`|
|`-homozygous`|`This function converts all heterozygous values`<br>`to Unknown.  Example: ./run_pipeline.pl -fork1 -`<br>`h mdp_genotype.hmp.txt -homozygous -export`<br>`homozygous.hmp.txt -runfork1`|
|`-numericalGenoTransform`<br>`<type>`|<br>`Performs genotype to numerical transform. <type>`<br>`can be collapse or separated.`|
|`-newCoordinates <map`<br>`filename>`|`This converts alignment to new coordinates`<br>`specified ingiven map file.`|
|`-synonymizer`|`Runs the Synonymizer using the input dataset.`|
|**Analysis**||
|`-glm`|`This takes a Phenotype dataset as input that is`<br>`usually the intersection of sequence data, trait`<br>`data, andpopulation structure(optional). `|
|`-glmOutputFile`<br>`<filename>`|`This sends GLM results to specified filename.`|
|`-glmMaxP <number>`|`This restricts the output file to entries with P`<br>`values no larger than number specified.`|
|`-glmPermutations`<br>`<number>`|`This sets the number of permutations. Default is`<br>`to not do runpermutations.`|
|`-mlm`|`This takes a Phenotype dataset as input(usually`|

8

||`the intersection of sequence data, trait data,`<br>`and population structure (optional)) and a`<br>`Kinship matrix.`|
|---|---|
|`-mlmVarCompEst <method>`|`Defines the Variance Component Estimation for`<br>`the previously specified -mlm. Method can be P3D`<br>`(default) or EachMarker.`|
|`-mlmCompressionLevel`<br>`<level>`|`Defines the Compression Level for the previously`<br>`specified -mlm. Level can be Optimum (default),`<br>`Custom, or None.`|
|`-mlmCustomCompression`<br>`<number>`|`This specifies the compression when compression`<br>`level is Custom. Default value is 1.0.`|
|`-mlmOutputFile`<br>`<filename>`|`This sends MLM results to specified filename.`|
|`-mlmMaxP <number>`|`This restricts the output file to entries with P`<br>`values no larger than number specified.`|
|`-diversity`|`Creates a Diversity Analysis step that uses an`<br>`Alignment as input`|
|`-diversityStartBase`<br>`<number>`|`This sets start base for the previously`<br>`specified -diversity. Default is 0.`|
|`-diversityEndBase`<br>`<number>`|`This sets end base for the previously specified`<br>`-diversity. Default is last site.`|
|`-diversitySlidingWin`|`This uses sliding window analysis for the`<br>`previously specified -diversity. `|
|`-diversitySlidingWinStep`<br>`<number>`|`This sets the sliding window step size for the`<br>`previously specified -diversity. Default is 100.`|
|`-diversitySlidingWinSize`<br>`<number>`|`This sets the sliding window size for the`<br>`previously specified -diversity. Default is 500.`|
|`-ld`|`Creates LinkageDisequilibriumPlugin. Uses`<br>`Alignment from previous step to analysis linkage`<br>`disequilibrium.`|
|`-ldPermNum <number>`|`This sets permutation number for the previously`<br>`specified -ld. Default is 1000.`|
|`-ldRapidAnalysis true |`<br>`false`|<br>`Sets whether to use rapid analysis for the`<br>`previously specified -ld. Default is true.`|
|`-ldWinSize <number>`|`Sets the window size for the previously`<br>`specified -ld. Default is 50.`|
|`-ldType <type>`|`Sets the LD type for the previously specified -`<br>`ld.  Options are All, SlidingWindow (Default),`<br>`and SiteByAll.`|
|`-ldTestSite <number>`|`Sets the test site for when LD type is set to`<br>`SiteByAll.`|
|`-ldHetTreatment <type>`|`Sets the LD Heterzygous Treatment Method. Type`<br>`can be Haplotype (Default - For Inbred Lines),`<br>`Homozygous (Uses only homozygous site -`<br>`heterozygotes set to missing), or Genotype (Not`<br>`Implemented Yet). `|

9

|`-ck`|`Calculates Kinship from Marker Data.`|
|---|---|
|`-ckModelHets <type>`<br>`-ckRescale true | false`|`Sets how to model heterozygotes.  Choose default`<br>`type RelateHomo (Related to Homozygotes) or`<br>`IndepState(Independent allele state).`<br>`Set whether to rescale results between 2 and 0.`<br>`Default is true.`|
|`-tree <clustering`<br>`method>`|`This creates a tree using given clustering`<br>`method: Neighbor (default) or UPGMA.  When`<br>`exporting, use -exportType Text to get text`<br>`version.`|
|`-treeSaveDistance true`<br>`| false`|`This saves the distance matrix of a tree.`<br>`Default is true.`|
|`-distanceMatrix`|`Calculate the distance matrix of given`<br>`Alignment.`|
|`-distMatrixRanges`|`Calculates genetic distances for given taxon in`<br>`specifiedphysicalposition ranges.`|
|`-distMatrixRangesLocus`<br>`<locus>`|`Locus that specified physical positions`<br>`corresponds.`|
|`-distMatrixRangesTaxon`<br>`<taxon>`|`Taxon of interest.`|
|`-distMatrixRangesPos`<br>`<pos1,pos2,pos3,…>`|`Specified physical positions that define ranges.`<br>`A comma should separate each one with no spaces.`|
|`-distMatrixRangesPosFile`<br>`<filename>`|`File with list of physical positions that define`<br>`ranges.  Individual positions should be`<br>`separated by white space.`|
|`-gs`|`Predicts phenotypes using ridge regression for`<br>`genomic selection.`|
|`-genotypeSummary`<br>`<types>`|`This generates summaries for alignment datasets.`<br>`Types should be a comma-separated list (with no`<br>`spaces) of the following (overall, site, taxa,`<br>`all). Example -genotypeSummary overall,site`|
|**Results**||
|`-td_csv <filename>`|`Writes (comma delimited) TableReport from`<br>`previous plugin in current pipeline to specified`<br>`filename.`|
|`-td_tab <filename>`|`Writes (tab delimited) TableReport from previous`<br>`plugin in current pipeline to specified`<br>`filename.`|
|`-td_gui`|`Displays TableReport from previous plugin in`<br>`currentpipeline in GUI.`|
|`-ldd <output type>`|`Creates LinkageDiseqDisplayPlugin. If output`<br>`type isgui, thisgraphically displays results`|

10

||`from a LinkageDisequilibriumPlugin. If output`<br>`type is png, gif, bmp, jpg, or svg, then an`<br>`image of that type is written to the output file`<br>`specified with -o.`|
|---|---|
|`-ldplotsize <num>`|`Optionally specify LD plot size. Example: 1000`<br>`will produce a 1000 x 1000 plot.  Default: 500.`<br>`This should follow the -ldd flag within the`<br>`currentpipeline segment.`|
|`-ldplotlabels true |`<br>`false`|`Optionally specify whether to show the LD`<br>`Plot labels.  DEFAULT: true. This should follow`<br>`the -ldd flag within the current pipeline`<br>`segment.`|
|`-o <output file>`|`This should follow the -ldd flag within the`<br>`currentpipeline segment.`|

11

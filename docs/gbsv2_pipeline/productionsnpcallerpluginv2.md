This plugin converts data from fastq and keyfile to genotypes, then adds these to a genotype file in VCF or HDF5 format.  VCF is the default output.  An HDF5 file may be requested by using the suffix ".h5" on the file used in the output file parameter.  Merging of samples to the same HDF5 output file may be accomplished by using the –ko option described below.

The parameters to this plugin are:

*      -batchSize <Batch Size> : Number of flow cells to process simultaneously. (Default: 8)
* 	-d <Max Divergence> : Maximum divergence (edit distance) between new read and previously mapped read (Default: 0 = perfect matches only) (Default: 0)
*	-db <Input GBS Database> : Input  Database file if using SQLite (REQUIRED)
*	-e <Enzyme> : Enzyme used to create the GBS library (REGQUIRED)
*	-eR <Ave Seq Error Rate> : Average sequencing error rate per base (used to decide between heterozygous and homozygous calls) (Default: 0.01)
*	-i <Input Directory> : Input directory containing fastq AND/OR qseq files (REQUIRED)
*	-k <Key File> : Key file listing barcodes distinguishing the sample (REQUIRED)
*      -ko <true | false> : Keep HDF5 genotypes file open for future runs that add more taxa or more depth (Default: false)
*       -kmerLength <Length of Kmer> : Lemgth of kmers to grab from fastQ sequences.  This value should match the kmerLength parameter value used in the GBSSeqToTagDBPlugin step of the pipeline.  Bad values may be stored in the HDF5 file if these values are inconsistent.  (Default: 64)
*       -minPosQS < Minimum Quality Score> : The minimum quality score a SNP position must have to be output to the HDF5 file.  The quality scores are loaded into the database via the tab-delimited file read in from the UpdateSNPPositionQualityPlugin.  A value of 0 indicates all positions should be processed. (Default: 0)
*       -mnQS < Minimum Quality Score> : The minimum quality score within the barcode and read length for a position to be accepted. This filters the read values from the fastQ files. (Default: 0)
*	-do <true | false> : depth output: write depths to the output HDF5 genotypes file (Default: true)
*	-o <Output Genotypes File> : Output (target) genotypes file to which is added new genotypes.  VCF format is the default.  if the file specified has suffix ".h5" output will be to an HDF5 file. (REQUIRED)

Calling this plugin from the command line would look like as below for VCF or HDF5 output:

```
#!java

./run_pipeline.pl -fork1 -ProductionSNPCallerPluginV2
-db /Users/lcj34/development/tempDir/GBS/Chr9_10-20000000/GBSv2.db
-e ApeKI -i /Users/lcj34/git/tassel-5-test/dataFiles/GBS/Chr9_10-20000000
-k /Users/lcj34/git/tassel-5-test/dataFiles/GBS/Pipeline_Testing_key.txt
-kmerLength 64
-o /Users/lcj34/development/GBS/productioHapMap_noKO.vcf -endPlugin -runfork1

./run_pipeline.pl -fork1 -ProductionSNPCallerPluginV2
-db /Users/lcj34/development/tempDir/GBS/Chr9_10-20000000/GBSv2.db
-e ApeKI -i /Users/lcj34/git/tassel-5-test/dataFiles/GBS/Chr9_10-20000000
-k /Users/lcj34/git/tassel-5-test/dataFiles/GBS/Pipeline_Testing_key.txt
-kmerLength 64
-o /Users/lcj34/development/GBS/productioHapMap_noKO.h5 -endPlugin -runfork1

```
To call ProductionSNPCallerPluginV2 from program code :

```
#!java

        new ProductionSNPCallerPluginV2()
                .enzyme("ApeKI")
                .inputDirectory(GBSConstants.GBS_INPUT_DIR)
                .inputGBSDatabase(GBSConstants.GBS_GBS2DB_FILE)
                .keyFile(GBSConstants.GBS_TESTING_KEY_FILE)
                .kmerLength(64)
                .outputHDF5GenotypesFile("junk1.vcf")
                .performFunction(null);

```

A note on depths reported in the HDF5 output file:  TASSEL takes the integer depths values and translates them to bytes to conserve memory.  Positive values stop at 127, higher values are translated to negative numbers.  Internally, when decoding a byte depth value, TASSEL applies an algorithm that interprets each of the byte values from -128 to +127 into a specific number.  Positive values are exact, negative values are log approximations. The depths encoding byte to integer values may be seen here: [DepthEncodingTable.txt](DepthEncodingTable.txt).

HDF5 viewing utilities are not able to correctly decode the TASSEL depth translation.  Using a basic "int myInt = (int)depthByte;" calculation will return an incorrect number with positive numbers capped at 127.  To see the depths accurately from the HDF5 file, you can import your HDF5 file to TASSEL, then export it as VCF with the "Keep Depth" box checked.  VCF tools may then allow you to see the depths correctly.  Alternately, you may output directly to VDF from ProductionSNPCallerPluginV2.

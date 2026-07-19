GBSSeqToTagDBPlugin takes fastQ files as input, identifies tags and the taxa in which they appear, and stores this data to a local database. It keeps only good reads having a barcode and a cut site and no N’s in the useful part of the sequence. It trims off the barcodes and truncates sequences that (1) have a second cut site, or (2) read into the common adapter.
The parameters to this plugin are:

* -c <Min Kmer Count> : Minimum kmer count (Default: 10)  Kmers that do not appear the minimum kmer count times are removed from the database.  This is a total of X count across all taxa, not in each individual taxon.
*    -db <Output Database File> : Output Database File (REQUIRED)
*    -deleteOldData < true | false> : Indicates whether user wants to delete old database data.  If the parameter is "true" or non-existent, and the output database parameter matches a file on the system, the file is deleted and a new file by that name is created when running the plugin.  If the parameter has value "false", the plugin appends data to the existing database if the "db" parameter specifies a file that exists.  Appending requires all taxa from the key file to match taxa currently existing in the database. An error message is thrown and processing is stopped if this is not the case.  (Default: true)
*    -e <Enzyme> : Enzyme used to create the GBS library if it differs from the one listed in the key file (REQUIRED)
*    -i <Input Directory>: Input directory containing FASTQ files in text or gzipped text. NOTE: Directory will be searched recursively and should be written WITHOUT a slash after its name. (REQUIRED)
*    -k <Key File> : Key file listing barcodes distinguishing the samples (REQUIRED)
*    -mnQS <Minimum quality score> : Minimum quality score within the barcode and read length to be accepted (Default: 0)
*    -kmerLength <Kmer Length> : Length of kmer to pull from fastQ sequences (Default: 64)  This value is used when calculating the sequence read length for kmer creation.  The kmerLength value must be chosen such that the longest barcode + kmerLength < read length.
*    -minKmerL <Minimum Kmer Length> : Minimum Kmer Length (Default: 20) After all fastQ files have been processed the stored kmers are searched for second cut sites.  If a second cut site is found in the kmer sequence, the sequence is truncated at that position.  The -minKmerL parameter is used to ensure a minimum length kmer exists after the second cut site is removed.
*    -mxKmerNum <Maximum Number of Kmers> : Maximum number of Kmers stored in the database. (Default: 50000000)
*    -batchSize <Number of flow cells processed simultaneously> : Number of flow cells to process simultaneously (Default: 8)

Additionally, for this and the other plugins in the GBS pipeline, you may add the -Xms and -Xmx parameters to indicate the minimum and maximum amount of memory Java should use.  The default values are platform specific.  What you need will depend on your machine's available memory and the size of the data you intend to process.

The minimum quality score (mnQS) and maximum kmer length (kmerLength) parameters are used together to determine which reads to keep from the fastQ files. The quality score of each position in the sequence is compared to the user supplied minimum quality score.  If the first low quality position found occurs at a position in the sequence before the specified maximum kmer length, the read is discarded.  In short, all sequence positions from 0 to (barcode length + mxKmerL) must have quality scores equal to or greater than the user supplied minimum quality score to be kept.

When using the "deleteOldData" flag to retain existing data in the database:  The taxa from the key files should be identical between the old and the new runs.  When running for the second time, the previous taxa list is deleted and replaced by the new now.  This is to prevent duplicate taxa, which would become an error in later steps.

The number of sequences determined to have low quality reads can be found in the log file by searching for the line "Total number of low quality reads=".

Here is an example of the command line execution for this plugin:

```
#!bash

./run_pipeline.pl -Xms100G -Xmx200G -fork1 -GBSSeqToTagDBPlugin -e ApeKI -i /Users/lcj34/git/tassel-5-test/dataFiles/GBS/Chr9_10-20000000/
  -db /Users/lcj34/development/GBS/Chr9_10-20000000/GBSV2.db
  -k /Users/lcj34/development/GBS/Pipeline_Testing_key.txt
  -kmerLength 64 -minKmerL 20 -mnQS 20 -mxKmerNum 100000000 -endPlugin -runfork1
```

Below is an example of how this would be called from within the code:

```
#!java
new GBSSeqToTagDBPlugin()
                .enzyme("ApeKI")
                .inputDirectory(GBSConstants.GBS_INPUT_DIR)
                .outputDatabaseFile(GBSConstants.GBS_GBS2DB_FILE)
                .keyFile(GBSConstants.GBS_TESTING_KEY_FILE)
                .kmerLength(64)
                .minKmerCount(5)
                .minimumQualityScore(20)
                .maximumKmerNumber(100000000)
                .performFunction(null);

```
If the output database file does not exist prior to running this step it will be created.  If it does exist, data will be integrated with the existing file based on the value of the "deleteOldData" parameter.

This is a very memory intensive step.  The default mxKmerNum is 50,000,000.  Each tag/kmer takes about 100 bytes of memory before being loaded with large amounts of data on distribution.

The JVM should be set to a large fixed memory size roughly 2-3X larger than that needed to load the initial kmer bytes into a map.
For example, if the mxKmerNum is set to 1,000,000,000 this requires 1000000000 X 100 = 10g for the map.
Then set the JVM to a large fixed size 2.5 fold larger.
-Xms25g -Xmx25g

We are continuing to explore better ways to control this size.  Please contact Ed and Lynn about issues.
#Fastq file name format
GBSv2 expects fastq files with 3, 4, or 5 underscore-delimited values.  The file names should be represented as one of the formats below:

* flowcell_lane_fastq.txt.gz
* flowcell_s_lane_fastq.txt.gz
* code_flowcell_s_lane_fastq.txt.gz

#FAQs

*How can I add my own enzymes to the GBSv2 pipeline?*

Pull the tassel-5 source code.
Edit file lib/enzymes.ini :  add your new enzymes
Compile and run

The code will pull in the new enzymes you've defined in the lib/enzymes.ini file

*What should I do if I get the error:  java.lang.IllegalArgumentException: duplicate key: Tag*

This indicates a single tag was found multiple times in the database.  The error generally occurs when there is an existing database file with the same name as the output database file specified for the corresponding GBSSeqToTagDBPlugin parameter.   Either use a different name for the output database file, or delete the existing file.  NOTE:  This error should no longer occur with the addition of the "deleteOldData" parameter as data integration will handle duplicates.

*Plugin returns without error but has not processed any files*

This could indicate none of the files in the specified input directory are represented in the key file.  Check to ensure your file names are of the format flowcell_lane_fastq.gz OR flowcell_s_lane_fastq.gz OR code_flowcell_s_lane_fastq.gz and that the flowcell and lanes appear in the key file.  NOTE:  the files may end with one of these patterns:  .fq, fq.gz, fastq, fastq.txt, fastq.gz, fastq.txt.gz, _sequence.txt, _sequence.txt.gz

*Can I use the "deleteOldData=false" to process/store fastQ files created with different enzymes into the same database?*

This is tricky.  The key files must contain the same taxa, which may not be the case if the libraries were made with different enzymes. We don't currently support multiple enzyme processing.  Our recommendation is to process the files separately into different databases.  Run the ProductionSNPCallerV2PluginV2 to create VCF files, and then merge the SNPs from this output.

*How can I add barcodes to my fastQ files so they are usable with this pipeline?*

An R script that creates new fastQ files with barcodes is the barcode_faker.R script, found [here](https://github.com/labroo2/rtassel_supp/blob/master/barcode_faker.R).  Note that after running this script, you will need to create a keyfile from the data exported by the R script.  See the script for details.

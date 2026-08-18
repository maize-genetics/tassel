TagExportToFastqPlugin retrieves distinct tags stored in the database and reformats them to a FASTQ file that can be read by the Bowtie2 or BWA aligner program.  This output file is input to the aligner, which creates a .sam file needed for calling SNPs further down in the GBS analysis pipeline.

The parameters to this plugin are:

* -c <Min Count>: Minimum count of reads for a tag to be output (Default: 1) This is a total of X count across all taxa, not in each individual taxon.
* -db <Input DB> : Input database file with tags and taxa distribution (REQUIRED)
* -o <Output File> : Output fastq file to use as input for BWA or Bowtie2 (REQUIRED)

To help identify tags whose sequence may have been altered by the aligner, the "seqname" field of the fastQ output file is populated with the tag sequence.  This value is preserved as the QNAME in the .sam file created by the aligners and can be used for mapping tags to the database if the tag sequence from the aligner has been altered.

A sample command line execution is this:

```
#!java

 ./run_pipeline.pl -fork1 -TagExportToFastqPlugin -db /Users/lcj34/git/tassel-5-test/tempdir/GBS/Chr9_10-20000000/GBSv2.db
 -o /Users/lcj34/git/tassel-5-test/tempDir/GBS/Chr9_10-20000000/tagsForAlign.fa.gz -c 1 -endPlugin -runfork1
```

To call TagExportToFastqPlugin() from within program code:

```
#!java

	        new TagExportToFastqPlugin()
            	.inputDB(GBSConstants.GBS_GBS2DB_FILE)
            	.outputFile(outTagFasta)
                .minCount(1)
            	.performFunction(null);
```

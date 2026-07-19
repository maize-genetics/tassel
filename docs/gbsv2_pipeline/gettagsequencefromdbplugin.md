Tag data is stored in the GBSv2 SQLite database in "blob" format. GetTagSequenceFromDBPlugin was created to provide a means for the user to grab and print the tag sequences from the database.

GetTagSequenceFromDBPlugin takes an existing GBSv2 SQLite database file as input and returns a tab-delimited file containing a list of Tag Sequences stored in the specified database file.  If the user wishes to query the existence of a specific tag in the database, the -tagSequence parameter may be supplied.  If this parameter is present, the code checks for this tag sequence among the stored tags.  The tab-delimited file will contain just 1 entry showing the tag (if it is found) or a single entry showing "NOT found: <tag sequence>".

The parameters to this plugin are:

* -tagSequence <Tag Sequence>: Verify this tag sequence exists in the database (Default: null)
* -db <Input DB> : Input database file with tags and taxa distribution (REQUIRED)
* -o <Output File> : Output tab-delimited text file that can be imported to Excel (REQUIRED)

A sample command line execution requesting all tags from the database is this:

```
#!java

./run_pipeline.pl -fork1 -GetTagSequenceFromDBPlugin -db /Users/lcj34/git/tassel-5-test/tempdir/GBS/Chr9_10-20000000/GBSv2.db -o /Users/lcj34/allTagFile.txt -endPlugin -runfork1
```

A sample command line execution querying the existence of a specific tag sequence in a database is this:

```
#!java

./run_pipeline.pl -fork1 -GetTagSequenceFromDBPlugin -db /Users/lcj34/git/tassel-5-test/tempdir/GBS/Chr9_10-20000000/GBSv2.db -o /Users/lcj34/singleTagFile.txt -tagSequence "CAGCTATTTTGAAGAAAGCAGGGTGCGAACATGTTAAGACTTTGGCCCAAGCCGAAGCTG" -endPlugin -runfork1
```

To call GetTagSequenceFromDBPlugin() from within program code:

```
#!java

       new GetTagSequenceFromDBPlugin()
        .inputDB(GBSConstants.GBS_GBS2DB_FILE)
        .outputFile(dbSingleTagOutFile)
        .tagSequence(myTag.sequence())
        .performFunction(null);
```

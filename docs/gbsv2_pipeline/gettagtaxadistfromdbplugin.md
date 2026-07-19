GetTagTaxaDistFromDBPlugin takes an existing GBSV2 SQLite database file as input and returns a tab-delimited file containing a list of tags and their taxa distribution as stored in the specified database file.

The parameters to this plugin are:

* -db <Input DB>:  Full path to the GBSv2 SQLite database from which to pull the tag/taxa distribution information (REQUIRED)
* -o <Output File>: Full path to output file.  This will be a tab-delimited text file that can be imported to Excel (REQUIRED)

A sample command line exception :

```
#!java

./run_pipeline.pl -fork1 GetTagTaxaDistFromDBPlugin -db /Users/micky/test/GBSv2.db -o /Users/micky/myResults/TagTaxaDistOutput.txt -endPlugin -runfork1

```

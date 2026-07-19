The GBSv2 key file has 4 required headers:  Flowcell, Lane, Barcode and FullSampleName.    These names are case-sensitive.  They may appear in any column of the key file.  The key file may contain additional columns of data useful to the user.  These additional columns become annotations for the taxon.

Entries in the FullSampleName column may be anything meaningful to the user as a unique identifier of the final genotype.  Entries in this column are used to facilitate merging the sample if it was run multiple times.   This will be useful to users that run replicate samples.  If the same entry appears in multiple rows of the FullSampleName column, all entries for those rows will be merged and processed together.  Example entries for the FullSampleName column are shown below:

```
#!python
B73:250020986
MO17
CO5FACXX:MO328:250021060
SpecialSample2

```

An example key file may be seen here:  [Pipeline_Testing_key.txt](Pipeline_Testing_key.txt)

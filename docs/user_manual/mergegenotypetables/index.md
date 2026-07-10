# Merge Genotype Tables
## Command

```
#!shell

./run_pipeline.pl -fork1 -h group1.hmp.txt -fork2 -h group2.hmp.txt -combine3 -input1 -input2 -mergeGenotypeTables -export group1_group2_merge.hmp.txt -runfork1 -runfork2 -runfork3
```

This is the most complex merge function, and can be considered as a union join across both sites and taxa. (The actual -union join only works across taxa.) The resulting genotype table will contain all unique sites and all unique taxa from across the input datasets. If a specific site-taxon combination isn’t present in any input dataset, the value is set to missing. If a specific site-taxon combination is present in more than one input file, the output will contain the last value processed. (That is, later values overwrite earlier values even if they conflict. There are plans to change this, but they have not been implemented yet.)

## Notes

* This maps to “Data -> Merge Genotype Tables” Menu on GUI.
* Error if duplicate site names in same file. (same as with other file loadings)
* Undefined taxa / site allele values are set to UNKNOWN.
* Duplicate taxa / site set to last Alignment processed.
* Sites are identified by Locus (chromosome), Physical Position, and Site Name

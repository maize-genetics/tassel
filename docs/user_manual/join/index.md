# Intersect Join
Command

```
#!shell

./run_pipeline.pl -fork1 -h group1.hmp.txt -fork2 -h group2.hmp.txt -combine3 -input1 -input2 -intersect -export group1_group2_intersect.hmp.txt -runfork1 -runfork2 -runfork3
```

This joins multiple data sets by the intersection of their taxa. Taxa must be present in both data sets to be included. Select multiple data sets using the CTRL key in conjunction with mouse clicks, and then click on the intersection button to join the data sets. Because this function uses taxa names to join data sets, any variation in taxa names can prevent proper joining. Taxa names can be made uniform by using the “Synonymizer”.

# Union Join
Command

```
#!shell

./run_pipeline.pl -fork1 -h group1.hmp.txt -fork2 -h group2.hmp.txt -combine3 -input1 -input2 -intersect -export group1_group2_union.hmp.txt -runfork1 -runfork2 -runfork3
```

This joins multiple data sets by a union of their taxa. Missing data will be inserted if taxa are missing from one data set. Select multiple data sets using the CTRL key in conjunction with mouse clicks, and then click on the union button to join the data sets. Because this function uses taxa names to join data sets, any variation in taxa names can prevent proper joining. Taxa names can be made uniform by using the “Synonymizer”.

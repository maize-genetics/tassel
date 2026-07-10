# Estimation of Kinship using genetic markers

While PCs can be used to capture major population subdivisions, kinship can be used to capture more subtle relationships. This section shows how to create a kinship matrix based on the same SNP data used to calculate PC’s.

1. Remove monomorphic sites: Highlight the genotype and choose Filter/Sites on the menu bar. Set the threshold on MAF to 0.05, check “Remove minor SNP status,” then click Filter.

2. Estimate kinship: Highlight the filtered genotype and click Analysis/Kinship. Leave “Scaled IBS” selected in the “Choose Kinship Method” dialog and click OK. A kinship matrix will be added to the data tree under Matrix category.

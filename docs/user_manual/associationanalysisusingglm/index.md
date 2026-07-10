# Association analysis using GLM

We use three files from the tutorial data set to perform association analysis using the GLM. The first file is mdp_genotype.hmp.txt, a set of SNPs scored at 3093 sites on 281 maize inbred lines. The second one is the population structure of 282 maize inbred lines (mdp_population_structure.txt). The last one contains phenotypes for three traits measured on 282 maize inbred lines (mdp_traits.txt). The statistical model is:

Flowering time = Population structure + Marker effect + residual

1. Remove monomorphic and low coverage sites: Highlight the mdp_genotype and click Filter/Sites on the menu bar. Set “Minimum Frequency” to 0.05, “Maximum Frequency”to 1.0, and “Minimum Count” to 150. Click Filter to create a filtered genotype data set.

2. Trait selection: Highlight the phenotype and click the menu item Filter/Traits. Uncheck all the traits except flowering time (dpoll). Make sure that the Type is set to Data. Click OK to create a filtered phenotype.

3. Covariate selection: The population structure is presented as the proportion of each population. There are three populations represented as Q1, Q2, and Q3. They sum to 100%. This creates linear dependency if we use all of them as covariates. While GLM can handle that properly, it will cause MLM to complain and refuse to complete your analysis. We can eliminate the dependency by removing one of the Q variables. In this demonstration, we exclude the last one. Highlight mdp_population_structure and click Filter/Traits. Uncheck the last population (Q3). Make sure that the Type is set to Covariate. Then click OK to create a filtered population structure data.

4. Joining data: Highlight the three filtered data sets by holding the Control key while selecting the individual data sets. Then click the menu item Data/Intersect Join to create a combined data set.

5. Association analysis: Highlight the joint data set then click the menu item Analysis/GLM to perform association analysis. Two reports will be added to the data tree.

One of the reports added to data tree is labeled “GLM_Stats_” followed by the name of the joint data. The following table shows an example of the GLM Stats output as viewed with Results/Table:

![GLM output 1.png](https://bitbucket.org/repo/bRpKEG/images/686848184-GLM%20output%201.png)

In addition to the information for traits and markers, the data set contains the following statistics:

* marker_F: F value from the F test on marker;
* p: P value from the F test on marker;
* marker_Rsq: R-squared for the marker after fitting other model terms (population structure);
* add_F: F value from the F test on the additive model
* add_p: P value from the F test on the additive model
* dom_F: F value from the F test of dominance (after fitting an additive model)
* dom_p: P value from the F test of dominance
* marker_df: Degree freedom of marker;
* marker_MS: Mean square of marker;
* error_df: Degree freedom of residual error;
* error_MS: Mean square of residual error;
* model_df: Degree freedom of model;
* model_MS: Mean square of model.

Clicking “marker_p” will sort the table by P value. The smallest P value is 3.5963x10-6. A reasonable significance threshold is 1.9x10-5, which is 5% after Bonferroni multiple test correction (0.05/2559). The denominator in the Bonferroni correction is the total number of SNPs tested. The association was significant.

The other data added to the data tree is labeled “GLM_Genotypes_” followed by the name of the joint data. For the most significant SNP (highlighted in the figure below), there were two genotypes (AA and GG). There are 220 lines with genotype AA and 41 lines with allele GG. For the trait dpoll (days to pollination), the difference between the two homozygotes was 3.86 days.

![GLM output 2.png](https://bitbucket.org/repo/bRpKEG/images/2398513383-GLM%20output%202.png)

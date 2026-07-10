# GLM (General Linear Model)

This function performs association analysis using a least squares fixed effects linear model.

TASSEL utilizes a fixed effects linear model to test for association between segregating sites and phenotypes. The analysis optionally accounts for population structure using covariates that indicate degree of membership in underlying populations. A main effects only model is automatically built using all variables in the input data. A separate model is built and solved for each trait and marker combination. Any factors, covariates, reps or locations are included in every model as main effects. How the data is used must be defined either in the input data files or using the Trait Filter after the data has been imported but before it has been joined with a genotype.

General Linear Model (GLM) can be run using a numeric data set only or using numeric data joined to genotype data. If only numeric data is selected, best linear unbiased estimates (BLUEs or least square means) will be generated for the taxa for each trait. [Note: only factors and covariates intended to control field variation should be included at this stage. Population structure covariates which are intended to control for marker effects should only be included when markers are also in the analysis.] If numeric data with genotypes are analyzed, each trait by marker combination will be tested and two reports will be produced, one containing trait by marker F-tests and the other containing allele estimates.

To run GLM, select a data set and then click the GLM button. A dialog box will pop-up to allow the user to indicate that a permutation test should be run and to allow the number of permutations to be changed. The permutation test will be run using the method suggested by Anderson and Ter Braak (2003), which calculates the predicted and residual values of the reduced model (contained all terms except markers) then permutes the residuals and adds them to the predicted values. When the GLM options dialog is closed, the user is presented with a dialog allowing the output to be saved to a file rather than stored in memory and displayed by TASSEL. This option is useful when the output is expected to be very large and risks exceeding available RAM.

The following table shows an example of the Marker Test output as viewed with Results/Table:

![GLM output 1.png](https://bitbucket.org/repo/bRpKEG/images/2836584079-GLM%20output%201.png)

The table shows the F-statistics and p-values for the requested F-tests for the main and additive models, and for the F-test for dominance after fitting the additive model. It also contains marker_Rsq, mean squares (MS) and degrees of freedom (DF) for the marker effect, for the model (corrected for the mean), and for error. If taxa are replicated (across reps or environments), then the markers are tested using the taxa within marker mean square. If taxa are unreplicated, then the residual mean square is used. Marker_Rsq is the marginal R-squared for the marker calculated as SS Marker (after fitting all other model terms) / SS Total, where SS stands for sum of squares. The following table shows an example of the Allele Estimates output as viewed with Results/Table:

![GLM output 3.png](https://bitbucket.org/repo/bRpKEG/images/1366794734-GLM%20output%203.png)

For each marker and trait combination, each marker allele is listed along with the chromosome and locus position of that marker, number of observations for taxa carrying that allele (Obs), the allele, and the estimate of the effect of that allele. Because of the way that GLM codes alleles, the last allele estimate for a marker is always zero and the other allele estimates are relative to that.

# Kinship Command Line

```
#!bash

./run_pipeline.pl -fork1 -importGuess mdp_genotype.hmp.txt -FilterSiteBuilderPlugin -siteMinAlleleFreq 0.05 -endPlugin -fork2 -importGuess mdp_traits.txt -fork3 -importGuess mdp_population_structure.txt -excludeLastTrait -combine5 -input1 -input2 -input3 -intersect -FixedEffectLMPlugin -endPlugin -export glm_output
```

```
#!bash

FixedEffectLMPlugin <options>
-phenoOnly <true | false> : Should the phenotype be analyzed with no markers and BLUEs generated? (BLUE = best linear unbiased estimate) (Default: false)
-saveToFile <true | false> : Should the results be saved to a file rather than stored in memory? It true, the results will be written to a file as each SNP is analyzed in order to reduce memory requirementsand the results will NOT be saved to the data tree. Default = false. (Default: false)
-siteFile <Statistics File> : The name of the file to which these results will be saved.
-alleleFile <Genotype Effect File> : The name of the file to which these results will be saved.
-maxP <max P value> : Only results with p <= maxPvalue will be reported. Default = 1.0. [0.0..1.0] (Default: 1.0)
-permute <true | false> : Should a permutation analysis be run? The permutation analysis controls the experiment-wise error rate for individual phenotypes. (Default: false)
-nperm <Number of Permutations> : The number of permutations to be run for the permutation analysis. (Default: 0)
-genotypeComponent <Genotype Component> : If the genotype table contains more than one type of genotype data, choose the type to use for the analysis. [Genotype, ReferenceProbability, AlleleProbability] (Default: Genotype)
-minClassSize <Minimum Class Size> : The minimum acceptable genotype class size. Genotypes in a class with a smaller size will be set to missing. (Default: 0)
-biallelicOnly <true | false> : Only test sites that are bi-allelic. The alternative is to test sites with two or more alleles. (Default: false)
-siteStatsOut <true | false> : Generate an output dataset with only p-val, F statistic, and number of obs per site for all sites. (Default: false)
-siteStatFile <Site Stat File> : Site Stat File
-appendAddDom <true | false> : If true, additive and dominance effect estimates will be added to the stats report for bi-allelic sites only. The effect will only be estimated when the data source is genotype (not a probability). The additive effect will always be non-negative. (Default: false)
```

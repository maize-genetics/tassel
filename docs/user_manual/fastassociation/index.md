#Fast Association

FastAssociation implements the method described in Shabalin (2012). The method provides an OLS (ordinary least squares) solution for fixed effect models and gives p-values equivalent to those calculated by GLM for markers. It can calculate values for either an additive + dominance model or for an additive only model. It can solve models that contain covariates and/or factors. To include one or more covariates or factors, simply add those to the phenotype data set as one would for GLM. For an analysis for a single phenotype, FastAssociation is not much faster than GLM and reports only r-square and p-values. For multiple phenotypes in a single analysis, it is much faster, because it avoids repeating most of the calculations for each phenotype separately. For eQTL (gene expression data treated as phenotypes) analysis, which can involve tens of thousands of phenotypes, FastAssociation can be faster by a factor of 100 or more.

_Two requirements for FastAssociation are (1) that phenotypes have no missing data and (2) that phenotypes and genotypes have been merged using an intersect join (NOT a union join)._ Requirement (2) is actually a consequence of requirement (1).

To run a FastAssociation analysis, choose a merged genotype-phenotype data set and click "Analysis/Fast Association" on the menu bar. The following dialog box will be displayed:

![Screen Shot 2015-02-18 at 1.34.10 PM.png](https://bitbucket.org/repo/bRpKEG/images/3043304264-Screen%20Shot%202015-02-18%20at%201.34.10%20PM.png)

MaxPValue : The maximum p-value that will be output by the analysis. (Default: 0.001)

Additive Only Model : Should an additive only model be fit? If true, an additive model will be fit. If false, an additive + dominance model will be fit. Default = false. (Default: false)

Genotype Component : If the genotype table contains more than one type of genotype data, choose the type to use for the analysis. [Genotype, ReferenceProbability, AlleleProbability] (Default: Genotype)

Write to file : Should the results be saved to a file rather than stored in memory? It true, the results will be written to a file as each SNP is analyzed in order to reduce memory requirementsand the results will NOT be saved to the data tree. Default = false. (Default: false)

Output File : The name of the file to which these results will be saved.

The analysis can also be run using the EqtlAssociationPlugin.

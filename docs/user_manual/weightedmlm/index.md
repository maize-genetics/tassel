# Weighted MLM

Weighted MLM and [MLM](../mlm/index.md) are very similar. They have the same results and interface except that weighted MLM can take a list of weights for the taxa.

Association analysis of replicated phenotypes using a mixed linear model can either use a one-stage or two-stage analysis. In a two-stage analysis, the phenotype data is evaluated first to produce a single estimate for each taxon (sample). Those estimates are then combined with genotypes to test associations. A one-stage analysis does both steps at once, fitting a model that has all of the individual observations. Because the MLM scales quadratically with the number of observations, a one-stage analysis is often not computationally feasible. If the phenotype data is balanced the one and two-stage methods produce the same result. Unbalanced data can result in biased effect estimates and loss of power for the two-stage model. The weighted MLM (Xue et al. 2017) weights each sample estimate by its error variance to produce results similar to the one-stage method but with the computational cost of the two-stage method.

The user interface and results for the weighted MLM are basically the same as the non-weighted version except that the weighted MLM can take a weight matrix as an additional argurment. If a weight matrix is not provided, then the analysis defaults to the non-weighted MLM. The weight matrix needs to be input as a phenotype file with two columns. The first column is taxa, and the second column is weight.

Using the TASSEL GUI: Do NOT join the weight matrix with the phenotype data. Instead select it as a third dataset when running the analysis. In other words, select the joint phenotype-genotype data, the kinship matrix, and the weight data then click Analysis/Association/WeightedMLM.

Example of a command line for running weighted MLM:

```
./run_pipeline.pl -Xmx4G -debug \
-fork1 -importGuess mdp_genotype.hmp.txt \
-FilterSiteBuilderPlugin -siteMinAlleleFreq 0.05 -endPlugin \
-fork2 -importGuess mdp_phenotype.txt \
-fork3 -importGuess kinship.txt \
-fork4 -importGuess weights.txt \
-combine5 -input1 -input2 -intersect \
-combine6 -input3 -input4 -input5 \
-WeightedMLMPlugin -endPlugin -export mlm_output
```

Reference:
Shang Xue, Funda Ogut, Zachary Miller, Janu Verma, Peter J. Bradbury, James B. Holland (2017) Comparison of one-stage and two-stage genome-wide association studies. bioRxiv. doi: https://doi.org/10.1101/099291

[User Manual Table of Contents](../index.md)

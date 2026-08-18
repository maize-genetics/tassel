# Numerical Imputation

Some methods e.g. PCA expect the data to have no missing values. But we often have some missing values which can be of following type -

* For phenotypic data, if the value of a trait for a taxon is not known/recorded.

* For genotypic data, if the value of some attribute is missing for a SNP.

Two methods are implemented to impute the missing values of data (both phenotypic and genotypic).

* Imputation by Mean - here the missing value is replaced by mean of the values for corresponding attribute.

* Imputation by k-nearest-neighbors - If data is missing for a taxon for one of the traits, the algorithm finds other taxa (neighbors) that are most like it for the non-missing traits. It uses the average of the neighbors to impute the missing data. Similarly if the value for some attribute for a SNP is missing, we compute k other SNPs which are most similar to it.
The default value of k is 5, which can be changed.
Also there are three different choices for distance measure to be for computing nearest neighbors.

1.  Euclidean
2.  Manhattan
3.  Cosine

The Euclidean distance is chosen as default.

## Running Numerical Imputation from the command line
Use the ImputationPlugin

Usage:  ImputationPlugin <options>
-ByMean <true | false> : If imputation is performed by computing mean of the respective column (Default: false)
-nearestNeighbors <Number of nearest neighbors to be evaluated> : Choice of k in k-nearest neighbors algorithm. Default is 5. (Default: 5)
-distance <Choose Distance type> : Distance choice for computing nearest neighbors. Default choice is Euclidean distance. (Default: Euclidean)

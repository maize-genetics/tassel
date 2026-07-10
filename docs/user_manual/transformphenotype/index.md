# Transform Phenotype

The Transform Phenotype menu item under the Data menu can be used to transform phenotype data. Some analyses like GLM and MLM assume that errors, in the case of GLM, or errors and random effects, in the case of MLM, are normally distributed. While the methods are robust to moderate departures from normality, if the data is obviously not normally distributed, transforming the data before analysis is one method that is used to generate more accurate statistical tests.

![Screen Shot 2015-02-11 at 8.44.59 AM.png](https://bitbucket.org/repo/bRpKEG/images/2034455887-Screen%20Shot%202015-02-11%20at%208.44.59%20AM.png)

The Transform Data dialog provides power, natural log, log2, and log10 transformations. The power transform raises each data point to the use supplied exponent. Non-integer exponents are allowed. The standardize option subtracts the mean of all the data and divides the result by the standard deviation. If the data set contains one or more factors, the dialog will display an additional list box for the factors. Selecting a factor will standardize data within each level of that factor by calculating the mean and variance separately for each factor level.

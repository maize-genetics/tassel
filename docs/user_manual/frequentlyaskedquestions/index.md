# Frequently Asked Questions

## What do I do if TASSEL misbehaves?

TASSEL is an open source software project hosted on SourceForge and has a bug
tracking list at <https://sourceforge.net/projects/tassel> where you can notify
the developer community of problems. In order for a bug to be fixed, we must be
able to replicate the problem, so it is important to document the steps that
produced the error. If the data you are working with is not too sensitive,
please include the files that were used in the faulty procedure. If you would
rather not post your data file on SourceForge, you may email it to one of the
software developers.

## Where do I turn for more information?

If you are having difficulty with a certain aspect of TASSEL, you can either
email one of the software developers listed at
<https://www.maizegenetics.net>, or check the TASSEL forum on
[SourceForge](https://sourceforge.net/projects/tassel), as another user may have
already addressed a similar question. There is also a TASSEL discussion group at
<https://groups.google.com/group/tassel>.

## How do I join the fun: TASSEL on SourceForge?

TASSEL is an open source project distributed under the GNU General Public
License. This means that the source code is available and the user is free to
modify the code to suit their particular needs. We welcome input from developers
and those who wish to become involved in the improvement of this software. The
project is hosted on [SourceForge](https://sourceforge.net/projects/tassel),
allowing anyone to access the most recent changes to the code. This setup makes
it convenient for anyone to add special functionality to TASSEL if they so
desire. It also serves as a good platform for anyone who wishes to become
involved in a bioinformatics software development project.

## When I click on the most current version of TASSEL Web Start, a previous version appears. What should I do?

The previous version of TASSEL Web Start was cached on your machine. To replace
it with the most current version, click the **Start** button in Windows,
followed by **Run**. Type `javaws` and then click **OK**. In the window that
opens, keep the most current version of TASSEL and delete the rest.

## What should I substitute for missing values in TASSEL?

For numerical data in version 3 format, use `NA` or `NaN`. For numerical data in
version 2 format, use `-999` for missing values. For SNP data, use `N`. Kinship
does not allow missing values.

## Is it possible to change data names in the Data Tree?

Yes. Click on the desired data name in the Data Tree, wait for one second, and
then click it again or immediately press `F2`. Rename the data set and then
press `Enter` to save the change.

## How can I create a TASSEL icon on the desktop?

Click **Start** on Microsoft Windows and select **Control Panel**, then
double-click **Java** to show the **Java Control Panel**. In the **Temporary
Internet Files** section, click the **View** button to show the **Java Cache
Viewer**. Move the mouse over the TASSEL application, right-click, and select
**Install Shortcuts**.

## Why do I get empty squares in MLM association analysis?

An empty square means null information. The major reasons include
non-convergence in the estimation of variance components, or that the statistic
in question was not calculated. For example, marker F, p, and R² are not
calculated when no marker is included in the model.

## Why should I exclude one column of the population structure?

For some methods of calculating population structure, such as the software
STRUCTURE, the population proportions sum to one. This produces linear
dependence between the population co-variates. While the algorithm used by GLM
tolerates that dependency, MLM will fail because the design matrix will not be
invertible. Excluding one column eliminates linear dependence between columns.
Using PC axes to represent population structure does not result in linear
dependency because all PC columns are guaranteed to be independent.

## Can kinship replace population structure?

Sometimes. For some traits and populations, the K-only model may be as good as
or better than the Q+K model. For others, Q+K may be superior. The Q-only model
is not as effective for controlling population structure as the alternatives.
Unfortunately, no general guidelines exist for predicting which model will
perform best. As a result, an investigator may wish to fit all three models and
compare the results. If eliminating false positives is very important, then it
may make sense to accept the most conservative model. However, if the objective
is to identify candidates for further study and the cost of following up on a
false lead is low, the most liberal model may be preferred.

## Why do TASSEL and SPAGeDi give different kinship estimates?

First, many algorithms exist to calculate kinship, and their estimates will
differ from one another. Second, the algorithm in TASSEL treats each genotype as
a haplotype. It is not recommended that TASSEL be used to generate a kinship
matrix from heterozygous genotypes. In the near future, the TASSEL kinship
algorithm will be modified to handle heterozygous diploids.

## Can I get Marker R square using SAS Proc Mixed or TASSEL MLM?

SAS Proc Mixed does not produce an R² statistic. MLM in TASSEL does. The user
manual describes how it is calculated.

## Does MLM find more associations than GLM?

Sometimes. MLM has higher statistical power than GLM and may detect more true
associations. When the tested genetic markers are confounded with kinship
structure, GLM does not correct for that as effectively as MLM and may produce
more false positives.

## Do I need multiple test correction for the p value from TASSEL?

Yes.

## Can TASSEL handle diploid genotype data?

While TASSEL accepts most common sequence alignment formats that handle polyploid
genotype data, including haploid and diploid, some analyses are not appropriate
for heterozygous data. GLM and MLM fit SNPs one at a time, treating each distinct
genotype as a separate class. This has the effect of fitting an additive plus
dominance model. Separating the two effects is under consideration.

## How do I cite TASSEL?

The paper that describes TASSEL [1] as a software package and the papers that
introduce specific methods implemented in TASSEL should be cited as appropriate,
such as the unified ("Q+K") approach, EMMA, compression of the mixed linear
model, and P3D. For example:

1. Linkage disequilibrium (D', R² and P value) were calculated by TASSEL [1].
2. Association analyses were performed with the mixed linear model approach [9] implemented by TASSEL [1].
3. GWAS was performed with the compressed mixed linear model approach [4, 9] carried out by TASSEL [1], which also implemented the EMMA [3] and P3D [4] algorithms to reduce computing time.

See the [References](../appendix/references.md) appendix for the full citations.

# Frequently Asked Questions

1. What do I do if TASSEL misbehaves?

TASSEL is an open source software project hosted on SourceForge and has a bug tracking list at http://sf.net/projects/tassel where you can notify the developer community of problems. In order for a bug to be fixed, we must be able to replicate the problem. Thus, it is important to document the steps that were taken that produced the error. If the data you are working with is not too sensitive, please include the files which were used in the faulty procedure. If you would rather not post your data file on SourceForge, you may email it to one of the software developers.

1. Where do I turn for more information?

If you are having difficulty with a certain aspect of TASSEL, you can either email one of the software developers listed at www.maizegenetics.net or you may check the TASSEL forum on SourceForge http://sf.net/projects/tassel), as another user may have already addressed a similar question. There is also a TASSEL discussion group at http://groups.google.com/group/tassel.

1. How do I join the fun: TASSEL on SourceForge?

TASSEL is an open source project distributed under the GNU general public license. This means that the source code is available and the user is free to modify the code to suit their particular needs. We welcome input from developers and those who wish to become involved in the improvement of this software. The project is hosted on SourceForge (http://sf.net/projects/tassel), thereby allowing anyone to access the most recent changes to the code. This setup makes it convenient for anyone to add special functionality to TASSEL if they so desire. It also serves as a good platform for anyone who wishes to become involved in a bioinformatics software development project.

1. When I click on the most current version of TASSEL web start, a previous version appears. What should I do?
The previous version of TASSEL web start was cached in your machine. To replace it with the most current version, click the Start button in Windows, followed by Run. Type javaws and then click OK. In the window that opens, keep the most current version of TASSEL and delete the rest.

1. What should I substitute for missing values in TASSEL?

For numerical data in version 3 format, use NA or NaN. For numerical data in version 2 format, use “-999” for missing values. For SNP data, use “N”. Kinship does not allow missing values.

1. Is it possible to change data names in the Data Tree?

Yes. Click on the desired data name in the Data Tree, wait for one second, and then click it again or immediately hit the F2 key. Rename the data set and then hit Enter to save the change.

1. How can I create a TASSEL icon on desktop?

Click “Start” on Microsoft Windows and select “Control Panel”, then double click Java to show “java Control Panel”. In “Temporary Internet Files” section, click “View” button show “Java Cache Viewer”. Move mouse over TASSEL application and click right button and select “Install Shortcuts”.

1. Why do I get empty squares in MLM association analysis?
The empty square means null information. The major reasons include non-convergence in the estimation of variance components or that the statistic in question was not calculated. For example, marker F, p, and R2 are not calculated when no marker is included in the model.

1.  Why should I exclude one column of the population structure?

For some methods of calculating population structure, such as the software STRUCTURE, the population proportions sum to one. This produces linear dependence between the population co-variates. While the algorithm used by GLM tolerates that dependency, MLM will fail because the design matrix will not be invertible Excluding one column eliminates linear dependence between columns.  Using PC axes to represent population structure does not result in linear dependency because all PC columns are guaranteed to be independent.

1. Can kinship replace population structure?
Sometimes. For some traits and populations, the K-only model may be as good as or better than the Q+K model. For others, Q+K may be superior. The Q-only model is not as effective for controlling population structure as the alternatives. Unfortunately, no general guidelines exist for predicting which model will perform best. As a result, an investigator may wish to fit all three models and compare the results. If eliminating false positives is very important, then it may make sense to accept the most conservative model. However, if the objective is to identify candidates for further study and the cost of following up on a false lead is low, the most liberal model may be preferred.

1. Why do TASSEL and SPAGeDi give different kinship estimates?

First, many algorithms exist to calculate kinship and their estimates will differ from one another. Secondly, the algorithm in TASSEL treats each genotype as a haplotype. It is not recommended that TASSEL be used to generate a  kinship matrix from heterozygous genotype. In the near future, the TASSEL kinship algorithm will be modified to handle heterozygous diploids.

1. Can I get Marker R square using SAS Proc Mixed or TASSEL MLM?

SAS Proc Mixed does not produce an R2 statistic. MLM in TASSEL does. The user manual describes how it is calculated.

1. Does MLM find more associations than GLM?

Sometimes. MLM has higher statistical power than GLM and may detect more true associations.. When the tested genetic markers are confounded with kinship structure , GLM does not correct for that as effectively as MLM and may produce more false positives

1. Do I need multiple test correction for the p value from Tassel?

Yes.

1. Can TASSEL handle diploid genotype data?

While TASSEL accepts most common sequence alignment formats which handle polyploid genotype data including haploid and diploid, some analyses are not appropriate for heterozygous data. GLM or MLM fit SNPs one at a time, treating each distinct genotype as a separate class. This has the effect of fitting an additive plus dominance model. Separating the two effects is under consideration.

1. How to cite TASSEL?

The paper that describes TASSEL1 as a software package and the papers that introduce  specific methods implemented in TASSEL should be cited as appropriate, such as the unified (“Q+K”) approach, EMMA, compression of mixed linear model and P3D. For example,:
1. Linkage disequilibrium (D’, R2 and P value) were calculated by TASSEL1.
2. Association analyses were performed with the mixed linear model approach9 implemented by TASSEL1.
3. GWAS was performed with the compressed mixed linear model approach4,9 carried by TASSEL1 which also implemented the EMMA3 and P3D4 algorithms to reduce computing time.

________________

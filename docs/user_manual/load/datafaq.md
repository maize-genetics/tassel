### 1. How do you load SSR data? ###
SSRs are no longer supported specifically in TASSEL 5.  TASSEL 4 supported them.  This has to do with the memory efficiency of various data structures.  One strategy is to recode the SSRs as SNPs, note collapsing the rare alleles is generally a good idea anyhow.

Code the SSRs alleles as follows:

* Most common = A
* 2nd most common = C
* 3rd most common = G
* 4th most common  = T
* 5th most common = +
* All other alleles = - (dash)

Also any allele class with less than 10 taxa(individuals) should be included in the all other alleles class.

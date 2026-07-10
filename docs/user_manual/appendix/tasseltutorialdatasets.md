# TASSEL Tutorial Data sets

https://tassel.bitbucket.io/docs/TASSELTutorialData5.zip

File name | Type | Format | Description
:-----:|:-----:|:-----:|:-----:
d8_sequence.phy | Genotype | Phylip Alignment | Sequence of dwarf8 gene with 2466 sites on 91 maize inbred lines. The data was described by the paper on the association between Dwarf8 and flowering time.
mdp_genotype.hmp.txt | Genotype | Hapmap | 3093 SNPs on 281 maize association inbred lines. The data was created by the PANZEA project funded by NSF. Details of the data can be found at http://www.panzea.org.
mdp_genotype.plk.ped, mdp_genotype.plk.map | Genotype | Plink | Same As Above
mdp_genotype.vcf | Genotype | VCF | Same As Above
mdp_kinship.txt | Kinship | Numerical Square Matrix | Kinship created by Yu et al.
mdp_population_structure.txt | Population structure | Numerical Trait Data | Population Structure of 282 maize inbred lines.
mdp_traits.txt | Phenotype | Numerical Trait Data | Three phenotypes, including flowering time, on 282 maize inbred lines.
mdp_phenotype.txt | Phenotype | Numerical Phenotype Data | Three phenotypes, one factor (location), and three covariates (population structure) on 282 maize inbred lines
282_GBS2.7_Depth7Mask.h5 | | HDF5 | 2.7 build GBS file that genotypes have been masked at sites with physical positions divisible by 7, and have a read depth of 7. 282 maize inbred lines. They are masked by simply setting to missing. The masked genotypes are held in corresponding files with the phase "maskKey", which have all genotypes set to missing except the ones that are masked.
NAM_RILs1-5_GBS2.7_Depth7Mask.h5 | | HDF5 | Same as above but with NAM RIL families 1-5

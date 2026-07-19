GenosToABHPlugin

This plugin was designed to convert genotypes of a biparental population that came out of the TASSEL GBS pipeline in a nucleotide-based format to a parent-based format. It also removes sites in which the parent genotype is missing, ambiguous or heterozygous. The output of this plugin is a .csv file with data arrangement specific to the R package ‘qtl’ [1].

The parameters to this plugin are:

* -o <Output file>: output filename, will be written as .csv
* -parentA <Parent A>: A plain text file containing the name of all samples from parent A as they are found in the input filename. One name per line.
* -parentB <Parent B>: A plain text file containing the name of all samples from parent B as they are found in the input filename. One name per line.
* -outputFormat <output format>: if "c", output will be A,H,B for parent A, het, and parent B: if "i", output will be 0,1,2 for parent A, het, parent B, if "r", output will be 0,0.5,1 for parent A, het, parent B.  This field it not required. If absent, the default is "c" - output will be in the form of A,H,B.

Additionally, an input file must be specified.  This may be any genotype file supported by TASSEL.  It is specified on the command line with either -h or -importGuess option.  When executing from program code, the file must be imported and passed as a DataSet parameter to the plugin's "performFunction(DataSet)" method.

A sample command line execution is this:

```
#!java

 ./run_pipeline.pl -h /Users/lcj34/genosToABH.hmp.txt -GenosToABHPlugin
 -o /Users/lcj34/ABHTestOutput.csv
 -parentA parentAFile.txt -parentB parentBFile.txt -outputFormat c -endPlugin
```

To call GenosToABHPlugin from within program code:

```
#!java

        GenotypeTable result = ImportUtils.readFromHapmap("genosToABH.hmp.txt");
        GenosToABHPlugin genosToABH = new GenosToABHPlugin(null, false);
        genosToABH.outfile(actualResultsFile);
        genosToABH.parentA(parentA);
        genosToABH.parentB(parentB);
        genosToABH.outputFormat(GenosToABHPlugin.OUTPUT_CHECK.c);
        DataSet myDS = genosToABH.performFunction(new DataSet(new Datum("inputFile", result, null),null));

```

Gory Details:

The plugin was designed with the idea of using GBS data from biparental populations to map QTLs. It makes several assumptions about the input file and outputs a format specific to the R package ‘qtl’.

Your input data file can be any genotype file supported by TASSEL. Since we expect data from the TASSEL GBS pipeline it is most likely hapmap. The file should contain all mapping individuals exactly once, plus any number of replicates for the two parents used to create that population. Make sure to get the exact names of the parental samples as they appear in the input file. Also we recommend to use at least two independent samples for each parents, especially if you have non-inbred, wild species.

The plugin will first try to find consensus genotypes for both parents for each site. A consensus genotype is called at a site when:

* at least one genotype was called for that parent
* the genotypes are the same for replicate parental samples
* the genotype is not heterozygous

If any of those assumptions is wrong the genotype is set to unknown. In the next step the genotypes of all mapping individuals are recoded, if both parents show a consensus genotype. If one or both parents do not have a consensus genotype that site is removed from analysis. The idea is to get rid of sites which do not allow calling parental genotypes or which would complicate downstream analysis. In our experience only a fraction of sites is removed by those steps and since GBS produces many sites it should not decrease the power of your analysis too much.

The output format is a genotype file that can be used directly in R/qtl. The first line of the file starts with “id” followed by the marker names. The second lines starts with “NA” followed by the chromosome where each marker is supposed to be. All following lines start with the unique name of the mapping individual, followed by the genotypes coded as ABH, one individual per line. The parents themselves are not included in the output file. To read the data to R/qtl use the function read.cross with the format = “csvs” option and supply a separate file containing the phenotypes with individuals in the same order in both files. It should also be possible without much effort to use Excel, a text editor or unix tools to reformat the genotypes for use in other software that expects ABH format.

[1] Broman KW, Wu H, Sen Ś, Churchill GA (2003) R/qtl: QTL mapping in experimental crosses. Bioinformatics 19:889-890

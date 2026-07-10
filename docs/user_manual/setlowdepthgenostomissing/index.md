# SetLowDepthGenosToMissingPlugin

This plugin sets genotypes to missing that have less than the user-specified minimum read depth. The genotypes must be from a source that provides read depth (for example, being read in from a VCF file will work, but not reading from a HapMap file). This plugin does not change the read depths.

The parameters to this plugin are:

* -minDepth <Integer>: minimum read depth below which genotypes are set to missing (restricted to the range 2 through 127).

Notice in the command line example below that a file is read into TASSEL before `SetLowDepthGenosToMissingPlugin` is called and that a file is written out after this plugin. `SetLowDepthGenosToMissingPlugin` itself does not directly read or write files. It simply handles genotypes in the middle. All depth information is kept, despite any changes in genotypes.

Example command line use:

```
#!code

 ./run_pipeline.pl -vcf myAppropriatelyNamed.vcf -SetLowDepthGenosToMissingPlugin -minDepth 7 -endPlugin -export myAppropriatelyNamedMinDepth7 -exportType VCF
```

As you may notice, the file is read in before calling `SetLowDepthGenosToMissingPlugin` and the file is written out after the plugin. `SetLowDepthGenosToMissingPlugin` itself does not read or write files. It simply does what

Example use within a program:

```
#!java

GenotypeTable inputGenos = ImportUtils.read("myAppropriatelyNamed.vcf");

GenotypeTable outputGenos = new SetLowDepthGenosToMissingPlugin()
                                    .minDepth(7)
                                    .runPlugin(inputGenos);

writeToVCF(outputGenos, "myAppropriatelyNamedMinDepth7.vcf", true); // true for keepDepth
```

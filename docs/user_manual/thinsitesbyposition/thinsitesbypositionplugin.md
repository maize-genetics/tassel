ThinSitesByPositionPlugin

This plugin was designed to thin out SNP sites based on a minimum distance in bp between adjacent sites. The idea is that SNP sites very close to each other or even on the same GBS tag do not provide additional useful information.

The parameters to this plugin are:

* -o <Output file>: Output genotype name
* -minDist <iMinimum distance>: Minimum distance in basepairs between adjacent sites.

Additionally an input file must be sent to the plugin.  This may be any genotype file supported by TASSEL.  It is defined on the command line via the -h or -importGuess option.  When running from program code, it must be imported to TASSEL and translated to a GenotypeTable file.

A sample command line execution is this:

```
#!java

 ./run_pipeline.pl -importGuess /Users/lcj34/genos.h5 -ThinSitesByPositionPlugin
 -o /Users/lcj34/thin40000.vcf -minDist 40000 -endPlugin
```

To call GenosToABHPlubin from within program code:

```
#!java

       GenotypeTable result = ImportUtils.readGuessFormat("maize_chr9_10.vcf");
        ThinSitesByPositionPlugin thinSitesByPosition = new ThinSitesByPositionPlugin(null, false);
        thinSitesByPosition.outfile(actualResults100);
        thinSitesByPosition.minDist(100);

        DataSet myDS = thinSitesByPosition.performFunction(new DataSet(new Datum("inputFile", result, null),null));

```

Gory Details:

The plugin assumes sites are ordered first by chromosome and then by position. For each chromosome, it will always keep the first site. Then each chromosome is traversed and starting from the first site, each following site is evaluated based on the distance to the previous site. If a site is within the user-specified minDist it is discarded. Otherwise it is kept and used as the new “anchor” to find sites to keep.

The input format can be any genotype file supported by TASSEL. The output format is dependent on the filename the user specified. Currently, .hmp.txt, .hmp.txt.gz, .vcf, .vcf.gz and .h5 files can be written.

The GBSv2 pipeline supports both the BWA and Bowtie2 alignment programs.  Details on running these programs are described below.

## Indexing with BWA

You must have BWA installed on your computer.  Download BWA from [here](http://bio-bwa.sourceforge.net).  For additional information on BWA consult the manual page [here](http://bio-bwa.sourceforge.net/bwa.shtml).

BWA will take the *.fa.gz file created by the TagExportToFastqPlugin, align the reads against a reference genome, and return a series of support files with the name given on the command line.

The first step to running with BWA is to create an index from the reference genome.  This can be done by running the bwa index command.

```
#!java

bwa index -a bwtsw referenceGenome/Zea_maysParsed.AGPv3.23.dna.genome.fa
```

The index command creates the series of files needed for the alignment phase.  It need only be run once.  The output files for the reference genome may be stored for later calls to the BWA aligner.  The support files created by the BWA index command have the same name as the biological sequence file used as input, with the suffixes *.fa.gz.bwt, *.fa.gz.pac, *.fa.gz.ann and *.ga.gz.sa.

Once you have indexed with BWA you can now align:

```
#!java

bwa aln -t 4 referenceGenome/Zea_maysParsed.AGPv3.23.dna.genome.fa tagsForAlign.fa.gz > tagsForAlign.sai

```

In the command above "tagsForAlign.fa.gz" is the output from running TagExportToFastqPlugin.

The file created from the "bwa aln ..." call is a ".sai" file containing the suffix array coordinates of all short reads loaded to BWA.  There is one more conversion to obtain the .sam file.

```
#!java

bwa samse referenceGenome/Zea_maysParse.AGPv3.23.dna.genome.fa tagsForAlign.sai tagsForAlign.fa.gz > tagsForAlign.sam
```

The last two letters of "samse" stand for "single-ended". Paired-end read alignment is possible but not used in GBS (which produces single end reads).  The "tagsForAlign.sam" file created will be used as input for the next step in the GBSv2 pipeline.

## Indexing with Bowtie2

The Bowtie2 download and manual can be obtained [here](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)

Before running the bowtie2 aligner, the index must be built.  This is done using the bowtie2-build command with the reference genome as input parameter and a named directory as an output parameter.  Here is an example command:

```
#!java

bowtie2-build ../referenceGenome/Zea_maysParsed.AGPv3.23.dna.genome.fa ../AGPv3/ZMAGP3_9103
```

The bowtie2-build command creates the index files needed for the alignment step.  Once the indexes have been created, the alignment command can be run.

```
#!java

bowtie2 -p 15 --very-sensitive -x ../AGPv3/ZMAGP3_9103 -U tagsForAlign.fa -S tagsForAlignFullvs.sam
```

The .sam file created in the aligner step in used as input to the GBSv2 plugin SAMToGBSdbPlugin.  Please consult the bowtie2 reference manual for details on the options to each command.

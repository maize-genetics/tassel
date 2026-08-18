UpdateSNPPositionQualityPlugin reads a quality score file to obtain quality score data for positions stored in the snpposition table.  The quality score file is a user created file that supplies quality scores for SNP positions. It is up to the user to determine what values should be associated with each SNP.  SNPQualityProfilerPlugin output provides data for this analysis, or the user may base quality scores on other data of his/her choice.

The parameters to this plugin are:

* -db< Input Database> : Input database file with SNP positions stored (REQUIRED)
* -qsFile <Quality Score File> : A tab delimited txt file containing headers CHROM( String), POS (Integer) and QUALITYSCORE(Float) for filtering.  (REQUIRED)

The quality scores are float values and may be any value the user desires.  These values may be used in the ProductionSNPCallerPluginV2 plugin to pull only positions whose value exceed a specified value.

To run this command from the command line:

```
#!java

./run_pipeline.pl -fork1 -UpdateSNPPositionQualityPlugin
-db /Users/lcj34/git/tassel-5-test/tempDir/GBS/Chr9_10-20000000/GBSv2.db -qsFile /Users/lcj34/myQsFile
-endPlugin -runfork1
```

To call the UpdateSNPPositionQualityPlugin from program code:

```
#!java

                new UpdateSNPPositionQualityPlugin()
                .gBSDBFile(GBSConstants.GBS_GBS2DB_FILE)
                .qsFile(myQSFile)
                .performFunction(null);
```

An example of a quality score file can be seen here: [snpQualityScoreFile_example.txt](snpQualityScoreFile_example.txt).  The "scores" are totally arbitrary and were created by a programmer (not a biologist) so they have no particular meaning.  The relevant positions and their scores should be determined by the biologist based on information obtained from the SNPQualityProfilerPlugin data or other means.  Note:  This is a tab-delimited file which is best viewed in Excel.

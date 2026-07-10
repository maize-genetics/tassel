# Tassel and R Integration

_Terry Casstevens (tmc46@cornell.edu), Peter Bradbury September 5, 2013_

## Calling R from Tassel

First, install rJava package as described below.

Run the command > system.file('jri',package="rJava")

For Java development (i.e. Tassel), use the Java API provided (http://www.rforge.net/JRI/index.html).

Once rJava has been installed, a folder named “library/rJava/jri” can be found in the R home directory. To find the R home directory, type R.home() on the command line in R. In the jri folder is a file named “run”, which is a script that can be used to set values for environment variables and call java, so that java.library.path and classpath are set so that the jar files and native library can be found by REngine. The simplest way to use R is to always use this script to start TASSEL.

A more flexible method is to make sure that the jri library file (libjri.jnilib on a Mac, a similarly named .dll file on Windows) is in a folder on the java.library.path. On a Mac, this can be accomplished by putting a copy or a link in /Library/Java/Extensions. On Windows, the folders listed in the PATH environment variable are the java.library.path. Alternatively, if <jri> = the full  path to the jri folder then adding -Djava.library.path = <jri> to the java command line will add that folder to the library path.

In addition, the R_HOME environment variable must be set to the R home directory prior to running the java command line. From a terminal command line, the command “R RHOME” returns the R home directory. Probably the most practical way to enable the use of R without much user intervention is to add commands to the start scripts to test for the existence of R and the jri library, then set R_HOME and use -D.library.path to add the jri folder to the library path. From within TASSEL, code should test for value of the R_HOME environment variable before attempting to call R.

The examples subdirectory of jri contains rtest.java, which has examples showing how to pass commands to R and how to move data to and from R.

## Calling Tassel from R

> install.packages("rJava")

Installing package(s) into ‘/Library/Frameworks/R.framework/Versions/2.15/Resources/library’ (as ‘lib’ is unspecified)

trying URL 'http://ftp.ussg.iu.edu/CRAN/bin/macosx/leopard/contrib/2.15/rJava_0.9-3.tgz' Content type 'application/x-gzip' length 791689 bytes (773 Kb) opened URL

================================================== downloaded 773 Kb

The downloaded binary packages are in

/var/folders/Jy/JybbW3n5GcORrFKVfrqYWk+++TI/-Tmp-//RtmpWosQU6/downloaded_packages

> library(rJava)

> .jinit()

> .jaddClassPath("/Users/terry/terry/tassel4.0_standalone/sTASSEL.jar")

> .jaddClassPath("/Users/terry/terry/tassel4.0_standalone/lib/")

> .jclassPath()

[1] "/Library/Frameworks/R.framework/Versions/2.15/Resources/library/rJava/java" [2] "/Users/terry/terry/tassel4.0_standalone/sTASSEL.jar" [3] "/Users/terry/terry/tassel4.0_standalone/lib"

> obj=.jnew("net.maizegenetics.pipeline.TerryTests") > result=.jcall(obj, "[[D", "runTest") > mat=sapply(result,.jevalArray) > mat [,1] [,2] [1,]    0    1 [2,]    1    2 [3,]    2    3

/* * TerryTests */ package net.maizegenetics.pipeline;

/** * * @author terry */ public class TerryTests { public static double[][] runTest() { double[][] result = new double[2][3]; for (int x = 0; x < 2; x++) { for (int y = 0; y < 3; y++) { result[x][y] = x + y; } } return result; } }

## Other Resources

Darren Wilkinson’s Blog http://darrenjw.wordpress.com/2011/01/01/calling-java-code-from-r/ rJava Manual - http://cran.r project.org/web/packages/rJava/rJava.pdf

rJava Home Page http://www.rforge.net/rJava/

import org.jetbrains.kotlin.gradle.dsl.JvmTarget
import org.jetbrains.kotlin.gradle.tasks.KotlinCompile
import org.jreleaser.model.Active
import java.util.Locale
import org.gradle.api.tasks.bundling.Zip
import org.gradle.api.tasks.bundling.Tar
import org.gradle.api.tasks.bundling.Jar
import org.gradle.jvm.application.tasks.CreateStartScripts

plugins {
    kotlin("jvm") version "2.1.20"
    id("org.jetbrains.dokka") version "2.0.0"
    id("org.jetbrains.kotlinx.kover") version "0.9.1"
    id("org.jreleaser") version "1.18.0"

    application
    `java-library`
    `maven-publish`
    signing
}

java {
    toolchain {
        languageVersion.set(JavaLanguageVersion.of(21))
    }
    withSourcesJar()
}

kotlin {
    jvmToolchain(21)
}

group = "net.maizegenetics"
version = "5.2.96"
description = "TASSEL is a software package to evaluate traits associations, evolutionary patterns, and linkage disequilibrium."
val kotlinVersion = "2.1.21"

repositories {
    mavenCentral()
    maven {
        url = uri("https://maven.scijava.org/content/repositories/public/") // needed for JHDF5
    }
}

// openchart:openchart:1.4.2 (a transitive dependency of forester / biojava-phylo) is no
// longer resolvable in any public repo and is not referenced directly by TASSEL source.
// Exclude it globally so the runtime classpath can be resolved.
configurations.all {
    exclude(group = "openchart", module = "openchart")
}

dependencies {
    testImplementation(kotlin("test"))
    implementation("org.apache.logging.log4j:log4j-api:2.21.1")
    implementation("org.apache.logging.log4j:log4j-core:2.21.1")
    implementation("com.google.guava:guava:22.0")
    implementation("org.apache.commons:commons-math3:3.4.1")
    implementation("commons-codec:commons-codec:1.10")
    implementation("org.biojava:biojava-core:6.0.4")
    implementation("org.biojava:biojava-alignment:6.0.4")
    implementation("org.biojava:biojava-phylo:4.2.12")
    implementation("org.jfree:jfreechart:1.0.19")
    implementation("org.jfree:jfreesvg:3.2")
    implementation("net.sf.trove4j:trove4j:3.0.3")
    implementation("org.ejml:ejml-ddense:0.41")
    implementation("org.xerial.snappy:snappy-java:1.1.8.4")
    implementation("javax.mail:mail:1.4")
    implementation("com.googlecode.json-simple:json-simple:1.1.1")
    implementation("org.glassfish:javax.json:1.0.4")
    implementation("org.xerial:sqlite-jdbc:3.39.2.1")
    implementation("com.github.samtools:htsjdk:2.24.1")
    implementation("org.ahocorasick:ahocorasick:0.2.4")
    implementation("org.postgresql:postgresql:42.6.0")
    implementation("org.apache.avro:avro:1.8.1")
    implementation("colt:colt:1.2.0")
    implementation("org.biojava.thirdparty:forester:1.039")
    implementation("cisd:jhdf5:19.04.1")
    implementation("org.jetbrains.kotlin:kotlin-stdlib-jdk8:$kotlinVersion")
    testImplementation("org.jetbrains.kotlin:kotlin-test:$kotlinVersion")
    implementation("org.ini4j:ini4j:0.5.4")
    implementation("it.unimi.dsi:fastutil:8.2.2")
}

// Application configuration
application {
    mainClass.set("net.maizegenetics.tassel.TASSELMainApp")
}

tasks.named<Zip>("distZip") {
    dependsOn(tasks.named("jar"), tasks.named("statisticsTest"))
}

tasks.named<Tar>("distTar") {
    dependsOn(tasks.named("jar"), tasks.named("statisticsTest"))
}

tasks.named<CreateStartScripts>("startScripts") {
    dependsOn(tasks.named("jar"), tasks.named("sourcesJar"))
}

// General tasks
tasks {
    // Set JAR file name and add to manifest along with classes
    withType<Jar> {
        archiveFileName.set("sTASSEL.jar")

        manifest {
            attributes(
                "Main-Class" to application.mainClass.get(),
                "Class-Path" to configurations.runtimeClasspath.get()
                    .joinToString(" ") { "lib/${it.name}" }
            )
        }
    }

    // Copy runtime dependencies into build/libs/lib
    register<Copy>("copyDependencies") {
        from(configurations.runtimeClasspath)
        into(layout.buildDirectory.dir("libs/lib").get().asFile)
    }

    // Ensure dependencies are copied after build
    named("build") {
        dependsOn("copyDependencies")
    }

    // TASSEL stores a couple of non-Java resources (SQLite schema DDL) next to the
    // sources in src/main/java. They must be on the runtime classpath so that
    // TagDataSQLite / RepGenSQLite can create their schemas via getResourceAsStream().
    processResources {
        from("src/main/java") { include("**/*.sql") }
    }

    // statisticsTest is the required CI gate; wire it into `check` so
    // `./gradlew check` enforces it without running the full non-blocking suite.
    named("check") {
        dependsOn("statisticsTest")
    }

    // Compile with Java 21 bytecode target
    withType<JavaCompile> {
        sourceCompatibility = "21"
        targetCompatibility = "21"
    }

    // Optional: Kotlin compile target
    withType<KotlinCompile> {
        compilerOptions {
            jvmTarget.set(JvmTarget.JVM_21)
        }
    }

    test {
        val baseArgs = mutableListOf("-Xmx10g")

        // Detect path to BLAS native library
        val overrideDir: String? = System.getenv("BLAS_LIB_PATH")
        if (!overrideDir.isNullOrBlank()) {
            baseArgs += "-Djava.library.path=$overrideDir"
        } else {
            val os = System.getProperty("os.name").lowercase(Locale.ROOT)
            val nativeDir = when {
                "mac" in os -> {
                    val intel = "/usr/local/opt/openblas/lib"
                    val silicon = "/opt/homebrew/opt/openblas/lib"
                    when {
                        file(intel).exists() -> intel
                        file(silicon).exists() -> silicon
                        else -> {
                            logger.warn("Neither $intel nor $silicon exists. BLAS tests may fail unless BLAS_LIB_PATH is set.")
                            ""
                        }
                    }
                }
                "linux" in os -> "/usr/lib/x86_64-linux-gnu"
                else -> {
                    logger.warn("Unrecognized OS: $os. BLAS tests may fail unless BLAS_LIB_PATH is set.")
                    ""
                }
            }

            if (nativeDir.isNotBlank()) {
                baseArgs += "-Djava.library.path=$nativeDir"
            }
        }

        exclude(
            // GBS Bucket C — kept excluded by design. Bucket A/B GBS tests now self-generate their
            // data via GBSSimData (no downloads/aligners) and run in this suite.
            // -- jhdf5 native library required (read/write .h5); native lib assumed unavailable --
            "**/analysis/gbs/SeqToTBTHDF5PluginTest.class",       // writes TBT HDF5
            "**/analysis/gbs/ModifyTBTHDF5PluginTest.class",      // pivots/merges TBT HDF5
            "**/analysis/gbs/DiscoverySNPCallerPluginTest.class", // v1 caller reads/writes HDF5 TOPM
            "**/analysis/gbs/ProductionSNPCallerPluginTest.class",// v1 production writes HDF5 genos
            "**/analysis/gbs/ProductionPipelineMainTest.class",   // full v1 pipeline via HDF5
            "**/analysis/gbs/v2/ProductionSNPCallerPluginV2Test.class", // v2 production writes HDF5 genos
            "**/analysis/gbs/v2/EvaluateSNPCallQualityOfPipelineTest.class", // needs HDF5 + large golden genos
            // -- byte-exact golden fixtures / pre-made SAM from v1 legacy (.cnt/.topm/.fq.gz) --
            "**/analysis/gbs/FastqToTagCountPluginTest.class",    // compares against golden .cnt
            "**/analysis/gbs/MergeMultipleTagCountPluginTest.class", // compares against golden .cnt
            "**/analysis/gbs/TagCountToFastqPluginTest.class",    // compares against golden .fq.gz
            "**/analysis/gbs/SAMConverterPluginTest.class",       // needs pre-made SAM + golden .topm
            // -- no real assertions / heavy converted-text fixtures --
            "**/analysis/gbs/v2/GBSv2BiologyCompareTest.class",   // comparison-only, needs large fixtures
            // -- hardcoded dev-machine paths / external tools (PEAR/BWA, maize AGPv4, /Users/lcj34) --
            "**/analysis/gbs/v2/RNADeMultiPlexSeqToDBPluginTest.class", // hardcoded dev paths + external tools
            "**/analysis/gbs/repgen/RepGenLoadSeqToDBPluginTest.class", // hardcoded dev paths + external data
            "**/analysis/gbs/repgen/RepGenAlignerPluginTest.class",     // hardcoded dev paths + external data
            "**/analysis/gbs/repgen/RepGenLDAnalysisPluginTest.class",  // hardcoded dev paths + external data
            "**/analysis/gbs/repgen/RampSeqAlignFromBlastTest.class",   // needs external BLAST output
            // External database integration — need live Postgres / MonetDB instances
            "**/analysis/gobii/*Test.class",
            "**/analysis/monetdb/*Test.class",
            // HDF5 native library — need jhdf5 native lib + large fixture files
            "**/LowLevelCopyOfHDF5Test.class",
            "**/SplitHDF5ByChromosomePluginTest.class",
            "**/TagsOnPhysMapHDF5Test.class",
            "**/DistanceMatrixHDF5Test.class",
            "**/BuildUnfinishedHDF5GenotypesPluginTest.class",
            // RNA sequencing — need sequencing DB fixtures
            "**/analysis/rna/*Test.class",
            // Hardcoded paths / heavy fixtures (not in test data release)
            "**/ThinSitesByPositionPluginTest.class",
            "**/LDKNNiImputationPluginTest.class",
            "**/GenomeFeatureBuilderTest.class",
            "**/BasicGenotypeMergeRuleTest.class",
            "**/GenomeAnnosDBQueryToPositionListPluginTest.class",
            // Performance variant with hardcoded /Volumes/... path; math covered by MLMTest
            "**/FastMultithreadedAssociationPluginTest.class",
            // Dev-only utility class with hardcoded file paths, not a runnable test
            "**/CreateFastaOrFastqFiles.class",
        )

        ignoreFailures = true // broad suite: keep non-blocking while pipeline/IO tests are fixed
        jvmArgs = baseArgs

        println(jvmArgs)
    }

    // ---------------------------------------------------------------------------
    // fetchTestData — download and extract the tassel_test_data release archive
    // into dataFiles/ (git-ignored). Run once after a clean checkout before
    // executing tests locally. CI already does this in coverage.yml.
    //
    // Usage:  ./gradlew fetchTestData
    // ---------------------------------------------------------------------------
    register("fetchTestData") {
        group = "verification"
        description = "Downloads and extracts the TASSEL test data archive into dataFiles/ if absent."
        doLast {
            val dataDir = file("dataFiles")
            if (dataDir.exists() && dataDir.list()?.isNotEmpty() == true) {
                logger.lifecycle("dataFiles/ already present — skipping download.")
                return@doLast
            }
            val tarball = file("tassel_test_data_v1.tar.gz")
            logger.lifecycle("Downloading TASSEL test data archive…")
            providers.exec {
                commandLine(
                    "curl", "-L", "--fail", "-o", tarball.absolutePath,
                    "https://github.com/maize-genetics/tassel_test_data/releases/download/v1.0.0/tassel_test_data_v1.tar.gz"
                )
            }.result.get()
            logger.lifecycle("Extracting TASSEL test data archive…")
            providers.exec {
                commandLine("tar", "-xzf", tarball.absolutePath)
            }.result.get()
            tarball.delete()
            logger.lifecycle("Test data extracted to dataFiles/")
        }
    }

    // ---------------------------------------------------------------------------
    // statisticsTest — enforced CI gate for statistical-correctness tests.
    //
    // Runs only the classes that verify TASSEL's numeric results (kinship, MLM,
    // GLM, PCA, LD, distance matrices, linear models) with ignoreFailures = false,
    // so failures are visible and block CI. The broad `test` task remains
    // non-blocking while pipeline/IO tests are still being fixed.
    //
    // Usage:  ./gradlew statisticsTest
    // ---------------------------------------------------------------------------
    val statisticsClasses = listOf(
        // Association / GWAS
        "net.maizegenetics.analysis.association.MLMTest",
        "net.maizegenetics.analysis.association.ReferenceProbabilityFELMTest",
        "net.maizegenetics.analysis.association.PhenotypeLMTest",
        "net.maizegenetics.analysis.association.GenomicSelectionPluginTest",
        "net.maizegenetics.analysis.association.EqtlAssociationPluginTest",
        "net.maizegenetics.analysis.association.DiscreteSitesTest",
        // Kinship / distance
        "net.maizegenetics.analysis.distance.KinshipTest",
        "net.maizegenetics.analysis.distance.CenteredIBSTest",
        "net.maizegenetics.analysis.distance.NormalizedIBSTest",
        "net.maizegenetics.analysis.distance.DominanceCenteredIBSTest",
        "net.maizegenetics.analysis.distance.DominanceNormalizedIBSTest",
        "net.maizegenetics.analysis.distance.IBSDistanceMatrixTest",
        "net.maizegenetics.analysis.distance.AMatrixPluginTest",
        // Linear models
        "net.maizegenetics.stats.linearmodels.ModelEffectTest",
        "net.maizegenetics.stats.linearmodels.SolveByOrtholgonalizingTest",
        // PCA
        "net.maizegenetics.stats.PCA.PrinCompTest",
        // Statistics utilities
        "net.maizegenetics.stats.statistics.FisherExactTest",
        // Linkage disequilibrium
        "net.maizegenetics.popgen.LinkageDisequilibriumTest",
        // Model fitting
        "net.maizegenetics.analysis.modelfitter.StepwiseAdditiveModelFitterTest",
        "net.maizegenetics.analysis.modelfitter.AdditiveSiteTest",
        // Matrix algebra
        "net.maizegenetics.matrixalgebra.Matrix.DoubleMatrixTest",
        // Numeric transforms
        "net.maizegenetics.analysis.numericaltransform.ImputationByMeanTest",
        "net.maizegenetics.analysis.numericaltransform.kNearestNeighborsTest",
        "net.maizegenetics.analysis.numericaltransform.SubtractPhenotypeByTaxaPluginTest",
        "net.maizegenetics.analysis.numericaltransform.AvgPhenotypeByTaxaPluginTest",
        "net.maizegenetics.analysis.numericaltransform.TransformDataPluginTest",
    )

    register<Test>("statisticsTest") {
        group = "verification"
        description = "Runs the statistical-correctness test gate with ignoreFailures = false."
        dependsOn("testClasses")
        testClassesDirs = sourceSets["test"].output.classesDirs
        classpath = sourceSets["test"].runtimeClasspath

        val baseArgs = mutableListOf("-Xmx10g")
        val overrideDir: String? = System.getenv("BLAS_LIB_PATH")
        if (!overrideDir.isNullOrBlank()) {
            baseArgs += "-Djava.library.path=$overrideDir"
        } else {
            val os = System.getProperty("os.name").lowercase(Locale.ROOT)
            val nativeDir = when {
                "mac" in os -> {
                    val intel = "/usr/local/opt/openblas/lib"
                    val silicon = "/opt/homebrew/opt/openblas/lib"
                    when {
                        file(intel).exists() -> intel
                        file(silicon).exists() -> silicon
                        else -> ""
                    }
                }
                "linux" in os -> "/usr/lib/x86_64-linux-gnu"
                else -> ""
            }
            if (nativeDir.isNotBlank()) baseArgs += "-Djava.library.path=$nativeDir"
        }
        jvmArgs = baseArgs

        filter {
            statisticsClasses.forEach { includeTestsMatching(it) }
        }

        ignoreFailures = false

        reports {
            html.outputLocation.set(layout.buildDirectory.dir("reports/tests/statisticsTest"))
            junitXml.outputLocation.set(layout.buildDirectory.dir("test-results/statisticsTest"))
        }
    }

    register("printVersion") {
        group = "help"
        description = "Prints the current project.version"
        doLast {
            println(project.version)
        }
    }
}

// Kover (coverage) tasks
kover {
    reports {
        verify {
            rule {
                bound {
                    minValue = 18
                }
            }
        }
    }
}




/**
 * Generates HTML files based on Javadoc-style comments. Supports automatic insertion of Jupyter notebook tutorials,
 * (see [tutorialInjector] for details). Supports insertion of images (see [imageInjector] for details).
 *
 * This was modified from the BioKotlin project.
 */
// configure Dokka’s HTML output directory so dokkaJar can find it
tasks.named<org.jetbrains.dokka.gradle.DokkaTask>("dokkaHtml") {
    outputDirectory.set(layout.buildDirectory.dir("dokka").get().asFile)
}

val dokkaHtml by tasks.getting(org.jetbrains.dokka.gradle.DokkaTask::class)

val dokkaJar by tasks.registering(Jar::class) {
    dependsOn(dokkaHtml)
    mustRunAfter(dokkaHtml)
    group = JavaBasePlugin.DOCUMENTATION_GROUP
    description = "TASSEL: ${project.version}"
    archiveClassifier.set("javadoc")
    from(dokkaHtml.outputDirectory)
}

publishing {
    publications {

        create<MavenPublication>("maven") {
            artifactId = "tassel"
            description = "net.maizegenetics:tassel:$version"

            from(components["java"])
            artifact(dokkaJar)

            versionMapping {
                usage("java-api") {
                    fromResolutionOf("runtimeClasspath")
                }
                usage("java-runtime") {
                    fromResolutionResult()
                }
            }

            pom {
                name.set("TASSEL")
                artifactId = "tassel"
                description.set("TASSEL is a software package to evaluate traits associations, evolutionary patterns, and linkage disequilibrium. ")
                url.set("https://www.maizegenetics.net/tassel")
                licenses {
                    license {
                        name.set("The Apache License, Version 2.0")
                        url.set("http://www.apache.org/licenses/LICENSE-2.0.txt")
                    }
                }
                developers {
                    developer {
                        id.set("tmc46")
                        name.set("Terry Casstevens")
                        email.set("tmc46@cornell.edu")
                    }
                    developer {
                        name.set("Ed Buckler")
                        email.set("esb33@cornell.edu")
                    }
                    developer {
                        name.set("Zack Miller")
                        email.set("zrm22@cornell.edu")
                    }
                    developer {
                        name.set("Lynn Johnson")
                        email.set("lcj34@cornell.edu")
                    }
                    developer {
                        name.set("Peter Bradbury")
                        email.set("pjb39@cornell.edu")
                    }
                    developer {
                        name.set("Brandon Monier")
                        email.set("bm646@cornell.edu")
                    }
                }
                scm {
                    connection.set("scm:git:git://github.com/maize-genetics/tassel.git")
                    developerConnection.set("scm:git:ssh://github.com/maize-genetics/tassel.git")
                    url.set("https://github.com/maize-genetis/tassel")
                }
            }
        }
    }

    signing {
        val signingKey: String? = System.getenv("JRELEASER_GPG_SECRET_KEY")
        val signingPass: String? = System.getenv("JRELEASER_GPG_PASSPHRASE")
        useInMemoryPgpKeys(signingKey, signingPass)
        sign(publishing.publications["maven"])
    }

    repositories {
        maven {
            url = layout.buildDirectory.dir("staging-deploy").get().asFile.toURI()
        }
    }
}

tasks.named("generateMetadataFileForMavenPublication") {
    dependsOn(tasks.named("dokkaJar"))
}

signing {
    useInMemoryPgpKeys(
        System.getenv("JRELEASER_GPG_SECRET_KEY"),
        System.getenv("JRELEASER_GPG_PASSPHRASE")
    )
    sign(publishing.publications["maven"])
}

jreleaser {
    signing {
        active.set(Active.ALWAYS)
        armored.set(true)
        setMode("MEMORY")
    }
    deploy {
        active.set(Active.ALWAYS)
        release {
            github {
                skipRelease = true
                skipTag = true
            }
        }
        maven {
            active.set(Active.ALWAYS)
            mavenCentral {
                signing {
                    active.set(Active.ALWAYS)
                    armored.set(true)
                    setMode("MEMORY")
                }
                create("sonatype") {
                    active.set(Active.ALWAYS)
                    url.set("https://central.sonatype.com/api/v1/publisher")
                    stagingRepository("build/staging-deploy")
                    sign = false
                }
            }
        }
    }
}

tasks.named("publish") {
    dependsOn("dokkaJar", "sourcesJar")
}

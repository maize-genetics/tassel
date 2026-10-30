import com.github.jengelman.gradle.plugins.shadow.transformers.Log4j2PluginsCacheFileTransformer
import org.jetbrains.kotlin.gradle.dsl.JvmTarget
import org.jetbrains.kotlin.gradle.tasks.KotlinCompile
import kotlinx.kover.gradle.plugin.dsl.CoverageUnit
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
    // 8.3.x is the line that supports Gradle 8.x; Shadow 9.3+ requires Gradle 9.
    id("com.gradleup.shadow") version "8.3.11"

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
version = "5.3.0"

// Release date reported by the About box and pipeline banner. This file is the
// single source of truth for both values: `generateVersionSources` compiles
// them into TasselBuildInfo, and docs/macros.py feeds them to the MkDocs build.
// Bump `version` and `versionDate` together and nothing else needs editing.
val versionDate = "August 19, 2026"

description = "TASSEL is a software package to evaluate traits associations, evolutionary patterns, and linkage disequilibrium."
val kotlinVersion = "2.1.21"

repositories {
    mavenCentral()
    maven {
        url = uri("https://maven.scijava.org/content/repositories/public/") // needed for JHDF5
    }
    maven {
        // Gradle ignores repositories declared in dependency POMs, so declare the JBoss
        // repo here to resolve 'openchart', a transitive dependency of forester.
        url = uri("https://repository.jboss.org/maven2/")
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
    testImplementation("junit:junit:4.13.2")
    implementation("com.formdev:flatlaf:3.7.1") // modern flat Swing Look-and-Feel (light/dark, HiDPI, macOS)
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

// src/main/java is registered as both a java root and a resource root (see the sourceSets
// block below), so every icon, HTML page and XML file under it reaches `allSource` twice.
// The sources jar packs `allSource` and refuses the second copy unless told what to do with
// it. Both copies resolve to the same file on disk, so keeping the first loses nothing.
tasks.named<Jar>("sourcesJar") {
    duplicatesStrategy = DuplicatesStrategy.EXCLUDE
}

// A plain `java -cp` launch has no bundle for macOS to read branding from, so the Dock shows
// the generic Java tile named "java". Only -Xdock: can override that; the
// apple.awt.application.name property TASSELMainApp sets reaches the screen menu bar but not
// the Dock. These stay scoped to `run` because the JVM rejects -Xdock: outright on other
// platforms, which would break the start scripts shipped in distZip/distTar.
tasks.named<JavaExec>("run") {
    if ("mac" in System.getProperty("os.name").lowercase(Locale.ROOT)) {
        jvmArgs("-Xdock:name=TASSEL", "-Xdock:icon=${file("icon.png").absolutePath}")
    }
}

// General tasks
tasks {
    // Set JAR file name and add to manifest along with classes. This must stay
    // scoped to the main `jar` task: applying it to every Jar task makes
    // sourcesJar and dokkaJar write to the same file, which published three
    // identical artifacts to Maven Central in 5.2.97.
    jar {
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

    // NOTE: the SQLite schema DDL (and other non-Java resources) under src/main/java are
    // already placed on the runtime classpath by the `sourceSets.main` resources config
    // further down, so no extra processResources copy is needed here — adding one
    // duplicates net/maizegenetics/dna/tag/*.sql at build time.

    // statisticsTest is the required CI gate; wire it into `check` so
    // `./gradlew check` enforces it without running the full non-blocking suite.
    named("check") {
        dependsOn("statisticsTest", "checkIgnoreReasons")
    }

    // ---------------------------------------------------------------------------
    // checkIgnoreReasons — every @Ignore in src/test must say why the test is
    // skipped: @Ignore("#123: reason") links the issue that tracks re-enabling it,
    // and @Ignore("benchmark: reason") marks a test meant to be run by hand.
    // Annotations inside comments or string literals are not counted.
    //
    // Usage:  ./gradlew checkIgnoreReasons
    // ---------------------------------------------------------------------------
    register("checkIgnoreReasons") {
        group = "verification"
        description = "Fails if an @Ignore in src/test has no issue reference or benchmark: prefix."
        val testSources = fileTree("src/test") { include("**/*.java", "**/*.kt") }
        inputs.files(testSources)

        doLast {
            // Blanks comments and string/char literal contents, keeping every offset
            // and newline, so a match in the result points at the same place in the source.
            fun codeOnly(src: String): String {
                val out = StringBuilder(src.length)
                fun blank(from: Int, to: Int) {
                    for (j in from until to) out.append(if (src[j] == '\n') '\n' else ' ')
                }
                var i = 0
                while (i < src.length) {
                    val end = when {
                        src.startsWith("//", i) ->
                            src.indexOf('\n', i).let { if (it < 0) src.length else it }
                        src.startsWith("/*", i) ->
                            src.indexOf("*/", i + 2).let { if (it < 0) src.length else it + 2 }
                        src.startsWith("\"\"\"", i) ->
                            src.indexOf("\"\"\"", i + 3).let { if (it < 0) src.length else it + 3 }
                        src[i] == '"' || src[i] == '\'' -> {
                            var j = i + 1
                            while (j < src.length && src[j] != src[i] && src[j] != '\n') {
                                if (src[j] == '\\') j++
                                j++
                            }
                            minOf(j + 1, src.length)
                        }
                        else -> -1
                    }
                    if (end < 0) {
                        out.append(src[i])
                        i++
                    } else {
                        blank(i, end)
                        i = end
                    }
                }
                return out.toString()
            }

            val annotation = Regex("""@Ignore\b""")
            val withReason = Regex("""@Ignore\s*\(\s*"(#\d+|benchmark:)""")
            val violations = testSources.files.sorted().flatMap { file ->
                val src = file.readText()
                annotation.findAll(codeOnly(src))
                    .filterNot { withReason.matchesAt(src, it.range.first) }
                    .map { match ->
                        val line = src.substring(0, match.range.first).count { it == '\n' } + 1
                        "${file.relativeTo(projectDir)}:$line"
                    }
                    .toList()
            }

            if (violations.isNotEmpty()) {
                throw GradleException(
                    "Every @Ignore needs a reason that starts with an issue reference or \"benchmark:\", " +
                        "for example @Ignore(\"#123: needs a live database\"). Missing in:\n" +
                        violations.joinToString("\n") { "  $it" }
                )
            }
        }
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

    // Mirrors GeneralConstants.TEMP_DIR, which many tests write their outputs into.
    // It is git-ignored, so it is absent on a fresh checkout; a plugin whose outFile
    // parameter points into a missing directory calls System.exit(1) and takes the
    // whole test worker down before any report is written. Create it for every Test
    // task rather than leaving it to whichever test class happens to run first — the
    // broad `test` suite only survives because some of its classes mkdirs() it in
    // setup, and none of those classes are in the statisticsTest gate.
    withType<Test> {
        doFirst {
            file("tempDir").mkdirs()
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
                // The archive carries macOS AppleDouble `._*` entries that GNU tar would
                // extract as real files, breaking tests that list data directories.
                commandLine("tar", "-xzf", tarball.absolutePath, "--exclude=._*")
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
        "net.maizegenetics.matrixalgebra.decomposition.EJMLEigenvalueDecompositionTest",
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

        // Without this the gate's CI log shows only Gradle's "finished with non-zero
        // exit value" line, naming neither the class that was running nor the assertion
        // that failed.
        testLogging {
            events("passed", "skipped", "failed")
            exceptionFormat = org.gradle.api.tasks.testing.logging.TestExceptionFormat.SHORT
        }

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

/**
 * Compiles `version` and `versionDate` into a generated `TasselBuildInfo` class so the About box,
 * pipeline banner, and log header report exactly what Gradle built. Before this, the same values
 * were hand-copied into `TASSELMainFrame` and silently went stale whenever only Gradle was bumped.
 */
val generateVersionSources by tasks.registering {
    group = "build"
    description = "Generates TasselBuildInfo from the project version and versionDate."

    val outputDir = layout.buildDirectory.dir("generated/sources/version/java/main")
    val buildVersion = project.version.toString()

    // Declaring these as inputs is what makes the task re-run on a version bump and stay
    // up-to-date otherwise; the output directory alone carries no record of the values used.
    inputs.property("version", buildVersion)
    inputs.property("versionDate", versionDate)
    outputs.dir(outputDir)

    doLast {
        val packageDir = outputDir.get().asFile.resolve("net/maizegenetics/tassel")
        packageDir.mkdirs()
        packageDir.resolve("TasselBuildInfo.java").writeText(
            """
            package net.maizegenetics.tassel;

            /**
             * Build metadata generated from `version` and `versionDate` in build.gradle.kts.
             *
             * <p>Generated by the Gradle `generateVersionSources` task. Do not edit: this file is
             * rewritten on every build and is not checked into version control.
             */
            public final class TasselBuildInfo {

                public static final String VERSION = "$buildVersion";
                public static final String VERSION_DATE = "$versionDate";

                private TasselBuildInfo() {
                }
            }

            """.trimIndent()
        )
    }
}

// Passing the task provider (rather than the directory) is what carries the task dependency
// across to compileJava and compileKotlin.
sourceSets.main {
    java.srcDir(generateVersionSources)

    // Icons, Home.html, the workflow presets and the SQLite schemas live beside the classes
    // that load them with getResource(). Maven packaged them through an explicit <resources>
    // include; Gradle's java plugin only looks in src/main/resources, so without this every
    // such lookup returns null at runtime. The include list keeps the JNI .c/.h files out.
    resources {
        srcDir("src/main/java")
        include("**/*.gif", "**/*.GIF", "**/*.png", "**/*.html", "**/*.sql", "**/*.xml")
    }
}

/**
 * Fat ("uber") JAR published to Maven Central as `tassel-<version>-jar-with-dependencies.jar`.
 * The classifier matches what maven-assembly-plugin produced through 5.2.96 so that existing
 * consumers of that coordinate keep working.
 */
tasks.shadowJar {
    archiveBaseName.set("tassel")
    archiveClassifier.set("jar-with-dependencies")

    manifest {
        attributes(
            "Main-Class" to application.mainClass.get(),
            // Shadow inherits the `jar` manifest, whose Class-Path points at the
            // external lib/ directory of the standalone distribution. Everything is
            // already inside this archive.
            "Class-Path" to "",
            // FlatLaf and log4j2 contribute META-INF/versions/9 classes, which the JVM
            // ignores unless the archive is flagged as multi-release.
            "Multi-Release" to "true"
        )
    }

    // Several dependencies (log4j2, biojava, avro) ship provider registrations that are
    // silently lost when only the first copy of each file survives the merge.
    mergeServiceFiles()
    transform(Log4j2PluginsCacheFileTransformer())
}

// The fat JAR is a publishing artifact only. Shadow's `application` integration
// registers a "shadow" distribution whose zip/tar land in the legacy `archives`
// configuration that `assemble` builds, which would drag the 70+ MB JAR into every
// `./gradlew build` in the jDeploy and nightly workflows. Removing those artifacts
// detaches the whole shadow distribution, including `shadowJar` itself, from the
// lifecycle; `shadowDistZip` and friends still work when invoked explicitly.
configurations["archives"].artifacts.removeIf { it.name == "${project.name}-shadow" }

// Kover (coverage) tasks
//
// Run coverage through the JaCoCo engine so that method-level boilerplate
// exclusion works: JaCoCo automatically ignores any class OR method annotated
// with an annotation whose simple name contains "Generated" (retention CLASS or
// RUNTIME). Our @GeneratedGuiBoilerplate marker satisfies this, so the auto-
// generated plugin accessors and GUI hook methods (getIcon/getButtonName/
// getToolTipText) are dropped. Kover's own annotatedBy filter only excludes
// whole classes (not methods) and does not function under JaCoCo, so we rely on
// JaCoCo's built-in filter for methods and on classes()/packages() for GUI types.
kover {
    useJacoco()

    // Coverage is measured from `test` alone. Kover's report tasks depend on every
    // instrumented Test task, so leaving statisticsTest instrumented would re-run it
    // in the coverage job even though it already runs as its own required CI job.
    currentProject {
        instrumentation {
            disabledForTestTasks.add("statisticsTest")
        }
    }

    reports {
        // Exclude GUI components so coverage reflects the exercised analysis/
        // pipeline logic rather than Swing wiring.
        filters {
            excludes {
                // Pure GUI / presentation packages.
                packages(
                    "net.maizegenetics.gui",
                    "net.maizegenetics.tassel",
                    "net.maizegenetics.progress",
                    "net.maizegenetics.analysis.chart",
                )

                // Swing widgets scattered through the analysis packages.
                classes("*Dialog", "*Panel", "*Component", "*DisplayPlugin")
            }
        }

        // Measure BRANCH coverage (decision paths exercised) instead of raw line
        // counting -- a far more meaningful signal for TASSEL's numeric and
        // pipeline logic. Not wired into `check`, but CI does run koverVerify,
        // and this bound is what the coverage job gates on. Keep it one point
        // below measured branch coverage (24.46% on develop, Oct 2026) and raise
        // it as coverage improves rather than leaving slack here.
        verify {
            rule {
                bound {
                    minValue = 23
                    coverageUnits = CoverageUnit.BRANCH
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
// Under Dokka's V2 plugin mode the V1 `dokkaHtml` task is disabled, so the javadoc JAR
// has to be packed from the V2 generator's output or it ships empty.
val dokkaGeneratePublicationHtml = tasks.named("dokkaGeneratePublicationHtml")

val dokkaJar by tasks.registering(Jar::class) {
    group = JavaBasePlugin.DOCUMENTATION_GROUP
    description = "TASSEL: ${project.version}"
    archiveClassifier.set("javadoc")
    from(dokkaGeneratePublicationHtml)
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
                    url.set("https://github.com/maize-genetics/tassel")
                }
            }
        }
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

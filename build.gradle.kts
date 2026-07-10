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

configurations.all {
    // openchart:openchart:1.4.2 (transitive of biojava forester) is not published to any
    // reachable repo. It is only used by forester's charting UI, which tests do not touch.
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
    dependsOn(tasks.named("jar"))
}

tasks.named<Tar>("distTar") {
    dependsOn(tasks.named("jar"))
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
            "**/analysis/gobii/*Test.class",
            "**/analysis/gbs/*Test.class",
            "**/analysis/gbs/v2/*Test.class",
            "**/analysis/gbs/repgen/*Test.class",
            "**/analysis/monetdb/*Test.class",
            "**/LowLevelCopyOfHDF5Test.class",
            "**/SplitHDF5ByChromosomePluginTest.class",
            "**/ThinSitesByPositionPluginTest.class",
            "**/LDKNNiImputationPluginTest.class",
            "**/GenomeFeatureBuilderTest.class",
            "**/BasicGenotypeMergeRuleTest.class",
            "**/DistanceMatrixHDF5Test.class",
            "**/TagsOnPhysMapHDF5Test.class",
            "**/FastMultithreadedAssociationPluginTest.class",
            "**/BuildUnfinishedHDF5GenotypesPluginTest.class",
            "**/GenomeAnnosDBQueryToPositionListPluginTest.class",
            "**/analysis/rna/*Test.class",
            "**/CreateFastaOrFastqFiles.class", // hard coded file paths (LCJ)
        )

        ignoreFailures = true // currently setting this to 'true' until we figure out failing tests
        jvmArgs = baseArgs

        println(jvmArgs)
    }

    // Runs the GBS (legacy) and GBSv2 test suites, which the main `test` task excludes.
    // Requires a `dataFiles/` dir in the project root (see notes/gbs-tests/).
    register<Test>("gbsTest") {
        description = "Runs GBS and GBSv2 tests."
        group = "verification"
        testClassesDirs = sourceSets["test"].output.classesDirs
        classpath = sourceSets["test"].runtimeClasspath
        useJUnit()
        include("**/analysis/gbs/*Test.class", "**/analysis/gbs/v2/*Test.class")
        jvmArgs = listOf("-Xmx10g")
        ignoreFailures = true
        testLogging {
            events("passed", "skipped", "failed")
            exceptionFormat = org.gradle.api.tasks.testing.logging.TestExceptionFormat.SHORT
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
                "Minimal line coverage rate as a percentage"
                bound {
                    minValue = 15
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
    outputDirectory.set(buildDir.resolve("dokka"))
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

package net.maizegenetics.analysis.association

import com.google.common.collect.HashMultiset
import net.maizegenetics.dna.snp.GenotypeTable
import net.maizegenetics.matrixalgebra.Matrix.DoubleMatrix
import net.maizegenetics.matrixalgebra.Matrix.DoubleMatrixFactory
import net.maizegenetics.matrixalgebra.decomposition.EigenvalueDecomposition
import net.maizegenetics.phenotype.*
import net.maizegenetics.plugindef.*
import net.maizegenetics.plugindef.GeneratePluginCode.*
import net.maizegenetics.util.TableReportBuilder
import org.apache.logging.log4j.LogManager
import java.awt.Frame
import java.util.*
import java.util.concurrent.Callable
import java.util.concurrent.Executors
import java.util.concurrent.Future
import javax.swing.ImageIcon
import kotlin.math.pow
import kotlin.random.Random
import net.maizegenetics.plugindef.PluginParameter
import net.maizegenetics.stats.linearmodels.*
import net.maizegenetics.util.TableReportUtils
import java.io.File
import java.util.concurrent.TimeUnit
import kotlin.math.max
import kotlin.math.min
import kotlin.math.roundToInt


/**
 * @author Samuel B. Fernandes
 * @author Peter Bradbury
 * @author Terry Casstevens
 * Created June 24, 2019
 */

class MultivariateStepwisePlugin(parentFrame: Frame?, isInteractive: Boolean) : AbstractPlugin(parentFrame, isInteractive) {

    private val myLogger = LogManager.getLogger(MultivariateStepwisePlugin::class.java)

    private var usePermutations = PluginParameter.Builder("usePerm", false, Boolean::class.javaObjectType)
            .description("Should permutations be used to set the enter and exit limits for stepwise regression? A permutation test will be used to determine the enter limit. The exit limit will be set to 2 times the enter limit.")
            .guiName("Use permutations")
            .build()

    private var numberOfPermutations = PluginParameter.Builder("nPerm", 1000, Int::class.javaObjectType)
            .description("The number of permutations used to determine the enter limit.")
            .guiName("Number of permutations")
            .dependentOnParameter(usePermutations)
            .build()

    private var permutationAlpha = PluginParameter.Builder("permAlpha", 0.05, Double::class.javaObjectType)
            .description("Type I errors will be controlled at this level.")
            .guiName("Alpha for permutations")
            .dependentOnParameter(usePermutations)
            .build()

    private var enterLimit = PluginParameter.Builder("enterLimit", 1e-5, Double::class.javaObjectType)
            .description("When p-value is the model selection criteria, model fitting will stop when the next term chosen has a p-value greater than the enterLimit. This value will be over-ridden the permutation test, if used.")
            .guiName("enterLimit")
            .dependentOnParameter(usePermutations, false)
            .build()

    private var exitLimit = PluginParameter.Builder("exitLimit", 2e-5, Double::class.javaObjectType)
            .description("During the backward step of model fitting if p-value has been chosen as the selection criterion, if the term in model with the highest p-value has a p-value > exitLimit, it will be removed from the model.")
            .guiName("exitLimit")
            .dependentOnParameter(usePermutations, false)
            .build()

    private var isNested = PluginParameter.Builder("isNested", false, Boolean::class.javaObjectType)
            .description("Should SNPs/markers be nested within a factor, such as family?")
            .guiName("")
            .build()

    private var nestingFactor = PluginParameter.Builder("nestFactor", null, String::class.java)
            .guiName("Nesting factor")
            .description("Nest markers within this factor. This parameter cannot be set from the command line. Instead, the first factor in the data set will be used.")
            .dependentOnParameter(isNested)
            .objectListSingleSelect()
            .build()

    private var myGenotypeTable: PluginParameter<GenotypeTable.GENOTYPE_TABLE_COMPONENT> = PluginParameter.Builder("genotypeComponent", GenotypeTable.GENOTYPE_TABLE_COMPONENT.Genotype, GenotypeTable.GENOTYPE_TABLE_COMPONENT::class.java)
            .genotypeTable()
            .range(GenotypeTable.GENOTYPE_TABLE_COMPONENT.values())
            .description("If the genotype table contains more than one type of genotype data, choose the type to use for the analysis.")
            .build()

    private var createManova = PluginParameter.Builder("Manova", true, Boolean::class.javaObjectType)
            .description("Create manova reports.")
            .guiName("Create manova reports")
            .build()

/*
    private var createEffects = PluginParameter.Builder("effects", true, Boolean::class.javaObjectType)
            .description("Create a report of marker effects based on the scan results.")
            .guiName("Create effects report")
            .build()
*/

    private var createStep = PluginParameter.Builder("step", true, Boolean::class.javaObjectType)
            .description("Create a report of the which markers enter and leave the model as it is being fit.")
            .guiName("Create step report")
            .build()

/*
    private var createResiduals = PluginParameter.Builder("residuals", false, Boolean::class.javaObjectType)
            .description("Create a phenotype dataset of model residuals for each chromosome. For each chromosome, the residuals will be calculated from a model with all terms EXCEPT the markers on that chromosome.")
            .guiName("Create residuals")
            .build()
*/

    private var writeFiles = PluginParameter.Builder("saveToFile", true, Boolean::class.javaObjectType)
            .description("Should the requested output be written to files?")
            .guiName("Write to files")
            .build()

    private var outputName = PluginParameter.Builder("savePath", "", String::class.java)
            .description("The base file path for the save files. Each file saved will add a descriptive name to the base name.")
            .guiName("Base file path")
            .outFile()
            .dependentOnParameter(writeFiles)
            .build()

    private var maximumNumberOfVariantsInModel = PluginParameter.Builder("maxQTN", 100, Int::class.javaObjectType)
            .description("maximum number of QTN to be fit in the model")
            .guiName("Maximum QTN Number")
            .build()

    private var runParallel = PluginParameter.Builder("parallel", true, Boolean::class.javaObjectType)
            .description("")
            .guiName("Run Parallel")
            .build()

    private var maxThreads = PluginParameter.Builder("threads", Runtime.getRuntime().availableProcessors(), Int::class.javaObjectType)
            .description("")
            .guiName("Number of threads")
            .build()


    private lateinit var myGenoPheno: GenotypePhenotype
    private lateinit var myDatasetName: String
    private lateinit var myFactorNameList: MutableList<String>
    private lateinit var randomGenerator: Random
    private var nestingFactorModelEffect : FactorModelEffect? = null

    //TableReport builders
    private lateinit var manovaReportBuilder: TableReportBuilder
    private lateinit var permutationReportBuilder: TableReportBuilder
    private lateinit var stepsReportBuilder: TableReportBuilder

    override fun preProcessParameters(input: DataSet?) {
        myFactorNameList = ArrayList<String>()
        //input data should be a single GenotypePhenotype
        val datumList = input?.getDataOfType(GenotypePhenotype::class.java)
        if (datumList!!.size != 1)
            throw IllegalArgumentException("Choose exactly one dataset that has combined genotype and phenotype data.")
        myGenoPheno = datumList[0].data as GenotypePhenotype
        myDatasetName = datumList[0].name

        myGenoPheno.phenotype().attributeListOfType(Phenotype.ATTRIBUTE_TYPE.factor).stream()
                .map { pa -> pa.name() }
                .forEach { myFactorNameList.add(it) }

        myFactorNameList.add("None")

        if (myFactorNameList.isEmpty()) myFactorNameList.add("None")
        nestingFactor = PluginParameter<String>(nestingFactor, myFactorNameList)

    }

    override fun processData(input: DataSet): DataSet? {
        manovaReportBuilder =
                TableReportBuilder.getInstance("Manova", arrayOf("SiteID", "Chr", "Position", "approx_F", "num_df", "den_df", "probF"))
        permutationReportBuilder =
                TableReportBuilder.getInstance("Empirical Null", arrayOf("p-value"))
        stepsReportBuilder =
                TableReportBuilder.getInstance("Steps", arrayOf("SiteID", "Chr", "Position", "action", "approx_F", "num_df", "den_df", "probF"))
        randomGenerator = Random(100)

        val start = System.nanoTime()

        var xR = DoubleMatrixFactory.DEFAULT.make(myGenoPheno.numberOfObservations(), 1, 1.0)
        val Y = createY()

        //add covariate and factor design matrices to the reduced model
        val modelEffectList = ArrayList<ModelEffect>()

        val covariateAttributeList = myGenoPheno.phenotype().attributeListOfType(Phenotype.ATTRIBUTE_TYPE.covariate)
        for (covariate in covariateAttributeList) {
            val covArray = (covariate as NumericAttribute).doubleValues()
            val covModelEffect = CovariateModelEffect(covArray, covariate.name())
            modelEffectList.add(covModelEffect)

            //add the effect design matrix to xR
            xR = xR.concatenate(covModelEffect.x, false)
        }
        val factorAttributeList = myGenoPheno.phenotype().attributeListOfType(Phenotype.ATTRIBUTE_TYPE.factor)
        for (factor in factorAttributeList) {
            val factorArray = (factor as CategoricalAttribute).allIntValues()
            val factorModelEffect = FactorModelEffect(factorArray, true, factor.name())
            modelEffectList.add(factorModelEffect)

            //add the effect design matrix to xR
            xR = xR.concatenate(factorModelEffect.x, false)
        }

        val startModelSize = factorAttributeList.size + covariateAttributeList.size

        //put nesting checks here
        //if isNested, make sure there is a valid nesting factor. Otherwise, throw an error
        if (isInteractive) {
            if (isNested()) {
                if (nestingFactor == null) throw java.lang.IllegalArgumentException("If model is nested, nestingFactor must be specified")
                val factorIndex = myGenoPheno.phenotype().attributeIndexForName(nestingFactor())
                if (factorIndex < 0) throw IllegalArgumentException("The nesting factor ${nestingFactor()} does not exist")
                val nestingAttribute = myGenoPheno.phenotype().attribute(factorIndex)
                nestingFactorModelEffect = FactorModelEffect((nestingAttribute as CategoricalAttribute).allIntValues(), false, nestingAttribute.name())
            }
        } else {
            if (isNested()) {
                if (factorAttributeList.size == 0) throw IllegalArgumentException("If model is nested, data must contain a factor variable")
                nestingFactorModelEffect = FactorModelEffect((factorAttributeList[0] as CategoricalAttribute).allIntValues(), false, factorAttributeList[0].name())
            }
        }

        val numberOfBaseEffects = modelEffectList.size

        //should a permutation test be run
        if (usePermutations()) {
            runPermutationTest(Y, xR)
        }

        val snpsAddedToModel = ArrayList<Int>()
        val step1 = forwardStep(Y, xR, modelEffectList, snpsAddedToModel)
        xR = step1

        val parallel = runParallel()
        if (parallel == null || parallel) {
            while (xR != null && modelEffectList.size <= maximumNumberOfVariantsInModel()) {
                var step2 = forwardStepParallel(Y, xR, modelEffectList, snpsAddedToModel)
                xR = step2
                if (xR != null && modelEffectList.size > 1) {
                    var result = backwardStep(Y, modelEffectList, xR, numberOfBaseEffects)
                    while (result.first) result = backwardStep(Y, modelEffectList, result.second, numberOfBaseEffects)
                    xR = result.second
                }
            }

        } else {
            while (xR != null && modelEffectList.size <= maximumNumberOfVariantsInModel()) {
                var step2 = forwardStep(Y, xR, modelEffectList, snpsAddedToModel)
                xR = step2
                if (xR != null && modelEffectList.size > 1) {
                    var result = backwardStep(Y, modelEffectList, xR, numberOfBaseEffects)
                    while (result.first) result = backwardStep(Y, modelEffectList, result.second, numberOfBaseEffects)
                    xR = result.second
                }
            }
        }

        myLogger.debug(String.format("ran analysis in %d ms.", (System.nanoTime() - start) / 1000000))

        calculateModelForManovaReport(Y, modelEffectList, startModelSize)

        val datumList = ArrayList<Datum>()

        if (createStep.value()) {
            val outputStep = stepsReportBuilder.build()
            datumList.add(Datum("Steps", outputStep, "Multivariate Stepwise regression results:\n" + "Model fitting steps\n"))

            if (writeFiles.value()) {
                val filenameStep = outputName.value() + "_steps.txt"
                TableReportUtils.saveDelimitedTableReport(outputStep, File(filenameStep))
            }
        }

        if (createManova.value()) {
            val outputMan = manovaReportBuilder.build()
            datumList.add(Datum("Manova", outputMan, "Multivariate Stepwise regression results:\\n Manova for the final model \\n"))

            if (writeFiles.value()) {
                val filenameMan = outputName.value() + "_Manova.txt"
                TableReportUtils.saveDelimitedTableReport(outputMan, File(filenameMan))
            }
        }

        if (usePermutations.value()) {
            val outputPer = permutationReportBuilder.build()
            datumList.add(Datum("Permutation", outputPer, "permutation report"))

            if (writeFiles.value()) {
                val filename = outputName.value() + "_permutations.txt"
                TableReportUtils.saveDelimitedTableReport(outputPer, File(filename))
            }
        }

        return DataSet(datumList, this)
    }

    fun forwardStep(Y: DoubleMatrix, xR: DoubleMatrix, modelEffectList: MutableList<ModelEffect>, snpsAdded: MutableList<Int>): DoubleMatrix? {
        val nSites = myGenoPheno.genotypeTable().numberOfSites()
        var minPval = 1.0
        var bestModelEffect: ModelEffect? = null
        lateinit var bestResult: List<Double>

        for (sitenum in 0 until nSites) {
            if (!snpsAdded.contains(sitenum)) {
                val genotypesForSite = imputeNsInGenotype(myGenoPheno.getStringGenotype(sitenum), randomGenerator)

                val modelEffect = if (isNested()) {
                    val snpEffect = FactorModelEffect(ModelEffectUtils.getIntegerLevels(genotypesForSite),
                            false, SnpData(myGenoPheno.genotypeTable().siteName(sitenum),
                            sitenum, myGenoPheno.genotypeTable().chromosomeName(sitenum),
                            myGenoPheno.genotypeTable().chromosomalPosition(sitenum)))
                    if (nestingFactorModelEffect == null) throw java.lang.IllegalArgumentException("")
                    NestedFactorModelEffect(snpEffect, nestingFactorModelEffect!!, snpEffect.id)
                } else {
                    FactorModelEffect(ModelEffectUtils.getIntegerLevels(genotypesForSite),
                            true, SnpData(myGenoPheno.genotypeTable().siteName(sitenum),
                            sitenum, myGenoPheno.genotypeTable().chromosomeName(sitenum),
                            myGenoPheno.genotypeTable().chromosomalPosition(sitenum)))
                }
                val siteDesignMatrix = xR.concatenate(modelEffect.x, false)
                val result = manovaPvalue(Y, siteDesignMatrix, xR)
                val pval = result[3]
                if (pval < minPval) {
                    minPval = pval
                    bestModelEffect = modelEffect
                    bestResult = result
                }
            }
        }

        if (minPval <= enterLimit() && bestModelEffect != null) {
            modelEffectList.add(bestModelEffect)
            snpsAdded.add((bestModelEffect.id as SnpData).index)
            println("${(bestModelEffect.id as SnpData).name}, pval = $minPval")

            //add a row to the step report
            //"SiteID", "Chr", "Position", "action", "approx_F", "num_df", "den_df", "probF"
            val snpdataForSite = bestModelEffect.id as SnpData


            stepsReportBuilder.add(arrayOf(snpdataForSite.name, snpdataForSite.chromosome, snpdataForSite.position,
                    "add", bestResult[0], bestResult[1], bestResult[2], bestResult[3]))
            return xR.concatenate(bestModelEffect.x, false)
        }

        return null
    }

    fun forwardStepParallel(Y: DoubleMatrix, xR: DoubleMatrix, modelEffectList: MutableList<ModelEffect>, snpsAdded: MutableList<Int>): DoubleMatrix? {
        //submit a range of snps to ManovaTester, which will return the SnpData and List<Double> result from manovaPvalue
        //for the snp with the lowest p value
        //repeat for all batches until all have been processed then pick the lowest p-value of those

        //set up an ExecutorService
        val nAvailableProcessors = Runtime.getRuntime().availableProcessors()
        val maxNumberThreads = maxThreads()
        val nThreads = if (maxNumberThreads != null) min(maxNumberThreads, nAvailableProcessors) else nAvailableProcessors


        val myExecutor = Executors.newFixedThreadPool(nThreads)


        val nSites = myGenoPheno.genotypeTable().numberOfSites()
        val batchSize = Math.max(1, nSites / nThreads / 100)

        var siteIndex = 0
        val futureList = ArrayList<Future<Pair<ModelEffect, List<Double>>>>()
        myLogger.debug("parallel execution using ${nThreads} threads, batch size = $batchSize")

        while (siteIndex < nSites) {
            val start = siteIndex
            siteIndex = Math.min(nSites, siteIndex + batchSize)
            futureList.add(myExecutor.submit(ManovaTester(start, siteIndex, Y, xR, snpsAdded, myGenoPheno, isNested(), nestingFactorModelEffect)))
        }

        lateinit var bestStatList: List<Double>
        lateinit var bestModelEffect: ModelEffect
        var minPval = 1.0
        futureList.forEach {
            val result = it.get(25, TimeUnit.SECONDS)
            if (result != null) {
                val statistics = result.second
                val probability = statistics[3]
                if (probability <= minPval) {
                    minPval = probability
                    bestModelEffect = result.first
                    bestStatList = statistics
                }
            }
        }

        myExecutor.shutdown()

        if (minPval <= enterLimit()) {
            modelEffectList.add(bestModelEffect)
            snpsAdded.add((bestModelEffect.id as SnpData).index)
            myLogger.debug("${(bestModelEffect.id as SnpData).name}, pval = $minPval")

            //add a row to the step report
            //"SiteID", "Chr", "Position", "action", "approx_F", "num_df", "den_df", "probF"
            val snpdataForSite = bestModelEffect.id as SnpData


            stepsReportBuilder.add(arrayOf(snpdataForSite.name, snpdataForSite.chromosome, snpdataForSite.position,
                    "add", bestStatList[0], bestStatList[1], bestStatList[2], bestStatList[3]))
            return xR.concatenate(bestModelEffect.x, false)
        }

        return null
    }

    fun backwardStep(Y: DoubleMatrix, modelEffectList: MutableList<ModelEffect>, originalXR: DoubleMatrix, numberOfBaseEffects: Int): Pair<Boolean, DoubleMatrix> {
        //test each snp one at a time, except for the last snp added
        //record the maximum pvalue
        //if max pvalue > exit limit remove the snp and return true
        //if max pvalue <= exit limit return false

        var maxPValue = 0.0
        var maxSnpIndex = 0
        var bestResult: List<Double>? = null

        val nEffectsInModel = modelEffectList.size
        val nObs = myGenoPheno.numberOfObservations()

        val intercept = DoubleMatrixFactory.DEFAULT.make(nObs, 1, 1.0)
        lateinit var maxXR: DoubleMatrix

        //do not reanalyze the last snp added
        for (snp in numberOfBaseEffects until nEffectsInModel - 1) {
            //xR includes all effects except snp
            var xR = intercept
            for (ndx in numberOfBaseEffects until nEffectsInModel) {
                if (ndx != snp) xR = xR.concatenate(modelEffectList[ndx].x, false)
            }
            val siteDesignMatrix = xR.concatenate(modelEffectList[snp].x, false)
            val result = manovaPvalue(Y, siteDesignMatrix, xR)
            val pval = result[3]
            if (pval > maxPValue) {
                maxPValue = pval
                maxSnpIndex = snp
                maxXR = xR
                bestResult = result
            }
        }

        if (maxPValue > exitLimit()) {
            //remove the effect and return the new xR and true
            val removedEffect = modelEffectList.removeAt(maxSnpIndex)

            if (bestResult != null) {
                //add a row to the step report
                //"SiteID", "Chr", "Position", "action", "approx_F", "num_df", "den_df", "probF"
                val snpdataForSite = removedEffect.id as SnpData

                stepsReportBuilder.add(arrayOf(snpdataForSite.name, snpdataForSite.chromosome, snpdataForSite.position,
                        "remove", bestResult[0], bestResult[1], bestResult[2], bestResult[3]))

            }
            return Pair(true, maxXR)
        } else {
            return Pair(false, originalXR)
        }

    }

    fun calculateModelForManovaReport(Y: DoubleMatrix, modelEffectList: MutableList<ModelEffect>, startModelSize: Int) {

        val nEffects = modelEffectList.size
        val nObs = myGenoPheno.numberOfObservations()

        for (snp in startModelSize until nEffects) {
            var xR = DoubleMatrixFactory.DEFAULT.make(nObs, 1, 1.0)
            for (effect in startModelSize until nEffects) {
                if (effect != snp) {
                    xR = xR.concatenate(modelEffectList[effect].x, false)
                }
            }

            val xF = xR.concatenate(modelEffectList[snp].x, false)
            val result = manovaPvalue(Y, xF, xR)
            modelEffectList[snp].id
            val snpdataForSite = modelEffectList[snp].id as SnpData
            //                "SiteID", "Chr", "Position", "approx_F", "num_df", "den_df", "probF"
            manovaReportBuilder.add(arrayOf(snpdataForSite.name, snpdataForSite.chromosome, snpdataForSite.position,
                    result[0], result[1], result[2], result[3]))

        }
    }

    fun runPermutationTest(Y: DoubleMatrix, xR: DoubleMatrix) {
        val nObs = Y.numberOfRows()
        val permutationPvalues = ArrayList<Double>()
        val index = (0 until nObs).toMutableList()
        val nPerm = numberOfPermutations()
        val permAlpha: Double = permutationAlpha() ?: throw IllegalArgumentException("permutation alpha level is null.")
        val permPercentile = max(0,(nPerm * permAlpha).roundToInt() - 1)

        for (i in  0 until numberOfPermutations()) {
            if (i % 10 == 0) myLogger.debug("Running permutation $i")
            //shuffle rows or Y
            index.shuffle(randomGenerator)
            val shuffledY = Y.getSelection(index.toIntArray(), null)

            //fit snps using a modification of forward step and return lowest p-value
            //add p-value to list of p-values
            val executeParallel = runParallel() ?: true
            if (executeParallel)permutationPvalues.add(permutationParallel(shuffledY, xR))
            else permutationPvalues.add(permutationTestNotParallel(shuffledY, xR))
        }

        //set the enter and exit limits based on the alpha percentile of pvalues
        permutationPvalues.sort()
        val permutationLimit = permutationPvalues[permPercentile]
        println("permutation limit = ${permutationLimit}")
        enterLimit(permutationLimit)
        exitLimit(permutationLimit)
        permutationPvalues.forEach { permutationReportBuilder.addElements(it); println(it) }
    }

    fun permutationParallel(Y: DoubleMatrix, xR: DoubleMatrix): Double {
        //submit a range of snps to ManovaTester, which will return the SnpData and List<Double> result from manovaPvalue
        //for the snp with the lowest p value
        //repeat for all batches until all have been processed then pick the lowest p-value of those

        //set up an ExecutorService
        val nAvailableProcessors = Runtime.getRuntime().availableProcessors()
        val maxNumberThreads = maxThreads()
        val nThreads = if (maxNumberThreads != null) min(maxNumberThreads, nAvailableProcessors) else nAvailableProcessors


        val myExecutor = Executors.newFixedThreadPool(nThreads)


        val nSites = myGenoPheno.genotypeTable().numberOfSites()
        val batchSize = Math.max(1, nSites / nThreads / 100)

        var siteIndex = 0
        val futureList = ArrayList<Future<Pair<ModelEffect, List<Double>>>>()

        val snpsAdded = ArrayList<Int>() //an empty list, because no snps have been added
        while (siteIndex < nSites) {
            val start = siteIndex
            siteIndex = Math.min(nSites, siteIndex + batchSize)
            futureList.add(myExecutor.submit(ManovaTester(start, siteIndex, Y, xR, snpsAdded, myGenoPheno, isNested(), nestingFactorModelEffect)))
        }

        lateinit var bestStatList: List<Double>
        var minPval = 1.0
        futureList.forEach {
            val result = it.get(25, TimeUnit.MINUTES)
            if (result != null) {
                val statistics = result.second
                val probability = statistics[3]
                if (probability <= minPval) {
                    minPval = probability
                    bestStatList = statistics
                }
            }
        }

        myExecutor.shutdown()


        return minPval;
    }

    fun permutationTestNotParallel(Y: DoubleMatrix, xR: DoubleMatrix): Double {
        val nSites = myGenoPheno.genotypeTable().numberOfSites()
        var minPval = 1.0

        for (sitenum in 0 until nSites) {
            val genotypesForSite = imputeNsInGenotype(myGenoPheno.getStringGenotype(sitenum), randomGenerator)

            val modelEffect = if (isNested()) {
                val snpEffect = FactorModelEffect(ModelEffectUtils.getIntegerLevels(genotypesForSite),
                        false, SnpData(myGenoPheno.genotypeTable().siteName(sitenum),
                        sitenum, myGenoPheno.genotypeTable().chromosomeName(sitenum),
                        myGenoPheno.genotypeTable().chromosomalPosition(sitenum)))
                if (nestingFactorModelEffect == null) throw java.lang.IllegalArgumentException("")
                NestedFactorModelEffect(snpEffect, nestingFactorModelEffect!!, snpEffect.id)
            } else {
                FactorModelEffect(ModelEffectUtils.getIntegerLevels(genotypesForSite),
                        true, SnpData(myGenoPheno.genotypeTable().siteName(sitenum),
                        sitenum, myGenoPheno.genotypeTable().chromosomeName(sitenum),
                        myGenoPheno.genotypeTable().chromosomalPosition(sitenum)))
            }
            val siteDesignMatrix = xR.concatenate(modelEffect.x, false)
            val result = manovaPvalue(Y, siteDesignMatrix, xR)
            val pval = result[3]
            if (pval < minPval) {
                minPval = pval
            }

        }

        return minPval
    }

    fun createY(): DoubleMatrix {
        val dataAttributeList = myGenoPheno.phenotype().attributeListOfType(Phenotype.ATTRIBUTE_TYPE.data)
        val nTraits = dataAttributeList.size
        val nObs = myGenoPheno.phenotype().numberOfObservations()
        val Y = DoubleMatrixFactory.DEFAULT.make(nObs, nTraits)
        for (t in 0 until nTraits) {
            val traitvals = dataAttributeList[t].allValues() as FloatArray
            traitvals.forEachIndexed { index, fl ->
                val traitvalue = fl.toDouble()
                if (traitvalue.isNaN()) throw java.lang.IllegalArgumentException("A trait value is missing. No missing trait values are allowed.")
                Y.set(index, t, traitvalue)
            }
        }
        return Y
    }


    override fun getIcon(): ImageIcon? {
        return null
    }

    override fun getButtonName(): String {
        return "Multivariate Stepwise"
    }

    override fun getToolTipText(): String {
        return "Multivariate Stepwise"
    }

    override fun getCitation(): String {
        return "Fernandes, SB; Casstevens, TM; Bradbury, PJ; and Lipka, AE. ... 2021"
    }

    override fun pluginUserManualURL(): String {
        return "https://bitbucket.org/tasseladmin/tassel­5­source/wiki/UserManual/..."
    }


    /**
     * Should permutations be used to set the enter and exit
     * limits for stepwise regression? A permutation test
     * will be used to determine the enter limit. The exit
     * limit will be set to 2 times the enter limit.
     *
     * @return Use permutations
     */
    fun usePermutations(): Boolean {
        return usePermutations.value()
    }

    /**
     * Set Use permutations. Should permutations be used to
     * set the enter and exit limits for stepwise regression?
     * A permutation test will be used to determine the enter
     * limit. The exit limit will be set to 2 times the enter
     * limit.
     *
     * @param value Use permutations
     *
     * @return this plugin
     */
    fun usePermutations(value: Boolean): MultivariateStepwisePlugin {
        usePermutations = PluginParameter(usePermutations, value)
        return this
    }

    /**
     * The number of permutations used to determine the enter
     * limit.
     *
     * @return Number of permutations
     */
    fun numberOfPermutations(): Int {
        return numberOfPermutations.value()
    }

    /**
     * Set Number of permutations. The number of permutations
     * used to determine the enter limit.
     *
     * @param value Number of permutations
     *
     * @return this plugin
     */
    fun numberOfPermutations(value: Int): MultivariateStepwisePlugin {
        numberOfPermutations = PluginParameter(numberOfPermutations, value)
        return this
    }

    /**
     * Type I errors will be controlled at this level.
     *
     * @return Alpha for permutations
     */
    fun permutationAlpha(): Double? {
        return permutationAlpha.value()
    }

    /**
     * Set Alpha for permutations. Type I errors will be controlled
     * at this level.
     *
     * @param value Alpha for permutations
     *
     * @return this plugin
     */
    fun permutationAlpha(value: Double?): MultivariateStepwisePlugin {
        permutationAlpha = PluginParameter(permutationAlpha, value)
        return this
    }

    /**
     * When p-value is the model selection criteria, model
     * fitting will stop when the next term chosen has a p-value
     * greater than the enterLimit. This value will be over-ridden
     * the permutation test, if used.
     *
     * @return enterLimit
     */
    fun enterLimit(): Double {
        return enterLimit.value()
    }

    /**
     * Set enterLimit. When p-value is the model selection
     * criteria, model fitting will stop when the next term
     * chosen has a p-value greater than the enterLimit. This
     * value will be over-ridden the permutation test, if
     * used.
     *
     * @param value enterLimit
     *
     * @return this plugin
     */
    fun enterLimit(value: Double): MultivariateStepwisePlugin {
        enterLimit = PluginParameter(enterLimit, value)
        return this
    }

    /**
     * During the backward step of model fitting if p-value
     * has been chosen as the selection criterion, if the
     * term in model with the highest p-value has a p-value
     * > exitLimit, it will be removed from the model.
     *
     * @return exitLimit
     */
    fun exitLimit(): Double {
        return exitLimit.value()
    }

    /**
     * Set exitLimit. During the backward step of model fitting
     * if p-value has been chosen as the selection criterion,
     * if the term in model with the highest p-value has a
     * p-value > exitLimit, it will be removed from the model.
     *
     * @param value exitLimit
     *
     * @return this plugin
     */
    fun exitLimit(value: Double): MultivariateStepwisePlugin {
        exitLimit = PluginParameter(exitLimit, value)
        return this
    }

    /**
     * Should SNPs/markers be nested within a factor, such
     * as family?
     *
     * @return Is Nested
     */
    fun isNested(): Boolean {
        return isNested.value()
    }

    /**
     * Set Is Nested. Should SNPs/markers be nested within
     * a factor, such as family?
     *
     * @param value Is Nested
     *
     * @return this plugin
     */
    fun isNested(value: Boolean): MultivariateStepwisePlugin {
        isNested = PluginParameter(isNested, value)
        return this
    }

    /**
     * Nest markers within this factor. This parameter cannot
     * be set from the command line. Instead, the first factor
     * in the data set will be used.
     *
     * @return Nesting factor
     */
    fun nestingFactor(): String {
        return nestingFactor.value()
    }

    /**
     * Set Nesting factor. Nest markers within this factor.
     * This parameter cannot be set from the command line.
     * Instead, the first factor in the data set will be used.
     *
     * @param value Nesting factor
     *
     * @return this plugin
     */
    fun nestingFactor(value: String): MultivariateStepwisePlugin {
        nestingFactor = PluginParameter(nestingFactor, value)
        return this
    }

    /**
     * If the genotype table contains more than one type of
     * genotype data, choose the type to use for the analysis.
     *
     * @return Genotype Component
     */
    fun genotypeTable(): GenotypeTable.GENOTYPE_TABLE_COMPONENT {
        return myGenotypeTable.value()
    }

    /**
     * Set Genotype Component. If the genotype table contains
     * more than one type of genotype data, choose the type
     * to use for the analysis.
     *
     * @param value Genotype Component
     *
     * @return this plugin
     */
    fun genotypeTable(value: GenotypeTable.GENOTYPE_TABLE_COMPONENT): MultivariateStepwisePlugin {
        myGenotypeTable = PluginParameter(myGenotypeTable, value)
        return this
    }

    /**
     * Create pre- and post-scan anova reports.
     *
     * @return Create anova reports
     */
    fun createManova(): Boolean {
        return createManova.value()
    }

    /**
     * Set Create anova reports. Create pre- and post-scan
     * anova reports.
     *
     * @param value Create anova reports
     *
     * @return this plugin
     */
    fun createManova(value: Boolean): MultivariateStepwisePlugin {
        createManova = PluginParameter(createManova, value)
        return this
    }

    /**
     * Create a report of marker effects based on the scan
     * results.
     *
     * @return Create effects report
     */
/*
    fun createEffects(): Boolean {
        return createEffects.value()
    }
*/

    /**
     * Set Create effects report. Create a report of marker
     * effects based on the scan results.
     *
     * @param value Create effects report
     *
     * @return this plugin
     */
/*
    fun createEffects(value: Boolean): MultivariateStepwisePlugin {
        createEffects = PluginParameter(createEffects, value)
        return this
    }
*/

    /**
     * Create a report of the which markers enter and leave
     * the model as it is being fit.
     *
     * @return Create step report
     */
    fun createStep(): Boolean {
        return createStep.value()
    }

    /**
     * Set Create step report. Create a report of the which
     * markers enter and leave the model as it is being fit.
     *
     * @param value Create step report
     *
     * @return this plugin
     */
    fun createStep(value: Boolean): MultivariateStepwisePlugin {
        createStep = PluginParameter(createStep, value)
        return this
    }

    /**
     * Create a phenotype dataset of model residuals for each
     * chromosome. For each chromosome, the residuals will
     * be calculated from a model with all terms EXCEPT the
     * markers on that chromosome.
     *
     * @return Create residuals
     */
/*
    fun createResiduals(): Boolean {
        return createResiduals.value()
    }
*/

    /**
     * Set Create residuals. Create a phenotype dataset of
     * model residuals for each chromosome. For each chromosome,
     * the residuals will be calculated from a model with
     * all terms EXCEPT the markers on that chromosome.
     *
     * @param value Create residuals
     *
     * @return this plugin
     */
/*
    fun createResiduals(value: Boolean): MultivariateStepwisePlugin {
        createResiduals = PluginParameter(createResiduals, value)
        return this
    }
*/

    /**
     * Should the requested output be written to files?
     *
     * @return Write to files
     */
    fun writeFiles(): Boolean {
        return writeFiles.value()
    }

    /**
     * Set Write to files. Should the requested output be
     * written to files?
     *
     * @param value Write to files
     *
     * @return this plugin
     */
    fun writeFiles(value: Boolean): MultivariateStepwisePlugin {
        writeFiles = PluginParameter(writeFiles, value)
        return this
    }

    /**
     * The base file path for the save files. Each file saved
     * will add a descriptive name to the base name.
     *
     * @return Base file path
     */
    fun outputName(): String {
        return outputName.value()
    }

    /**
     * Set Base file path. The base file path for the save
     * files. Each file saved will add a descriptive name
     * to the base name.
     *
     * @param value Base file path
     *
     * @return this plugin
     */
    fun outputName(value: String): MultivariateStepwisePlugin {
        outputName = PluginParameter(outputName, value)
        return this
    }

    /**
     * maximum number of QTN to be fit in the model
     *
     * @return Maximum QTN Number
     */
    fun maximumNumberOfVariantsInModel(): Int {
        return maximumNumberOfVariantsInModel.value()
    }

    /**
     * Set Maximum QTN Number. maximum number of QTN to be
     * fit in the model
     *
     * @param value Maximum QTN Number
     *
     * @return this plugin
     */
    fun maximumNumberOfVariantsInModel(value: Int): MultivariateStepwisePlugin {
        maximumNumberOfVariantsInModel = PluginParameter(maximumNumberOfVariantsInModel, value)
        return this
    }

    /**
     * Run Parallel
     *
     * @return Run Parallel
     */
    fun runParallel(): Boolean? {
        return runParallel.value()
    }

    /**
     * Set Run Parallel. Run Parallel
     *
     * @param value Run Parallel
     *
     * @return this plugin
     */
    fun runParallel(value: Boolean?): MultivariateStepwisePlugin {
        runParallel = PluginParameter(runParallel, value)
        return this
    }

    /**
     * Number of threads
     *
     * @return Number of threads
     */
    fun maxThreads(): Int? {
        return maxThreads.value()
    }

    /**
     * Set Number of threads. Number of threads
     *
     * @param value Number of threads
     *
     * @return this plugin
     */
    fun maxThreads(value: Int?): MultivariateStepwisePlugin {
        maxThreads = PluginParameter(maxThreads, value)
        return this
    }

}

/*fun main(args : Array<String>) {
    generate(MultivariateStepwisePlugin::class.java)
}*/

data class BetaValue(val B: DoubleMatrix, val H: DoubleMatrix)
data class SnpData(val name: String, val index: Int, val chromosome: String, val position: Int)

class ManovaTester(val start: Int, val end: Int, val Y: DoubleMatrix, val xR: DoubleMatrix,
                   val snpsAdded: MutableList<Int>, val genoPheno: GenotypePhenotype, val isNested: Boolean,
                   val nestingFactorModelEffect : FactorModelEffect?) : Callable<Pair<ModelEffect, List<Double>>> {

    val randomGenerator = Random(start)

    override fun call(): Pair<ModelEffect, List<Double>>? {
        var minPval = 1.0
        var bestModelEffect: ModelEffect? = null
        var bestResult: List<Double>? = null

        for (sitenum in start until end) {
             if (!snpsAdded.contains(sitenum)) {
                val genotypesForSite = imputeNsInGenotype(genoPheno.getStringGenotype(sitenum), randomGenerator)

                val modelEffect = if (isNested) {
                    val snpEffect = FactorModelEffect(ModelEffectUtils.getIntegerLevels(genotypesForSite),
                            false, SnpData(genoPheno.genotypeTable().siteName(sitenum),
                            sitenum, genoPheno.genotypeTable().chromosomeName(sitenum),
                            genoPheno.genotypeTable().chromosomalPosition(sitenum)))
                    nestingFactorModelEffect
                    NestedFactorModelEffect(snpEffect, nestingFactorModelEffect!!, snpEffect.id)
                } else {
                    FactorModelEffect(ModelEffectUtils.getIntegerLevels(genotypesForSite),
                            true, SnpData(genoPheno.genotypeTable().siteName(sitenum),
                            sitenum, genoPheno.genotypeTable().chromosomeName(sitenum),
                            genoPheno.genotypeTable().chromosomalPosition(sitenum)))
                }
                val siteDesignMatrix = xR.concatenate(modelEffect.x, false)
                val result = manovaPvalue(Y, siteDesignMatrix, xR)
                val pval = result[3]
                if (pval <= minPval) {
                    minPval = pval
                    bestModelEffect = modelEffect
                    bestResult = result
                }
            }
        }

        if (bestModelEffect == null || bestResult == null) return null
        return Pair(bestModelEffect, bestResult)
    }
}

/**
 * Y = Phenotypic data
 * xF = Full model
 * xR = Reduced model
 * returns Wilk's lambda pvalue
 */
fun manovaPvalue(Y: DoubleMatrix, xF: DoubleMatrix, xR: DoubleMatrix): List<Double> {
    //total SQ
    val YtY: DoubleMatrix = Y.crossproduct()
    //number of variables
    val p = Y.numberOfColumns().toDouble()
    //full model
    val hF = calcBeta(Y, xF).H
    //Residual (Error) matrix
    val E: DoubleMatrix = YtY.minus(hF)
    //reduced model
    val hR = calcBeta(Y, xR).H
    //adjusted hypothesis matrix
    val hA: DoubleMatrix = hF.minus(hR)
    //Wilks lambda
    val eigen: EigenvalueDecomposition = E.inverse().mult(hA).eigenvalueDecomposition
    var lambda = 1.0
    for (i in 0..(eigen.eigenvalues.size - 1)) {
        lambda *= 1 / (1 + eigen.eigenvalues[i])
    }
    //Degree of Freedom
    val df = xF.columnRank().toDouble() - xR.columnRank().toDouble()//data[data.names[0]].asStrings().distinctBy {it.hashCode()}.size.toDouble()
    if (df == 0.0) {
        return listOf(1.0, 0.0, 0.0, 1.0)
    }
    val t: Double
    if ((p.pow(2) + df.pow(2) - 5) > 0) {
        t = Math.sqrt((p.pow(2) * df.pow(2) - 4) / (p.pow(2) + df.pow(2) - 5))
    } else {
        t = 1.0
    }

    val ve = Y.numberOfRows().toDouble() - xF.columnRank().toDouble()

    val r = ve - (p - df + 1) / 2
    val f = (p * df - 2) / 4
    val fCalc = ((1 - lambda.pow(1 / t)) / (lambda.pow(1 / t))) * ((r * t - 2 * f) / (p * df))
    val num: Double = p * df
    val den: Double = (r * t) - (2 * f)
    val pvalue: Double = LinearModelUtils.Ftest(fCalc, num, den)
//        return pvalue
    return listOf(fCalc, num, den, pvalue)
}

fun calcBeta(Y: DoubleMatrix, X: DoubleMatrix): BetaValue {
    val XtX: DoubleMatrix = X.crossproduct()
    val XtY: DoubleMatrix = X.crossproduct(Y)
    val XtXinv: DoubleMatrix = XtX.generalizedInverse()
    val B: DoubleMatrix = XtXinv.mult(XtY)
    val H: DoubleMatrix = B.transpose().mult(XtY)
    return BetaValue(B, H)
}

fun imputeNsInGenotype(genotypes: Array<String>, randomGenerator: Random): Array<String> {
    //TODO write test to make sure this works
    val alleleCounter = HashMultiset.create<String>()
    for (allele in genotypes) if (!allele.equals("N")) alleleCounter.add(allele)
    val totalCount = alleleCounter.count().toDouble()
    val alleleProportionList = ArrayList<Pair<Double, String>>()
    var cumulativeSum = 0
    for (entry in alleleCounter.entrySet()) {
        cumulativeSum += entry.count
        alleleProportionList.add(Pair(cumulativeSum / totalCount, entry.element))
    }

    return genotypes.map {
        if (it.equals("N")) {
            val ran = randomGenerator.nextDouble()
            var ndx = 0
            while (ran >= alleleProportionList[ndx].first) ndx++
            alleleProportionList[ndx].second
        } else it
    }.toTypedArray()
}


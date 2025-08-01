package net.maizegenetics.analysis.data

import com.google.common.collect.ImmutableRangeMap
import com.google.common.collect.Range
import com.google.common.collect.RangeMap
import net.maizegenetics.dna.map.Position
import net.maizegenetics.plugindef.AbstractPlugin
import net.maizegenetics.plugindef.DataSet
import net.maizegenetics.plugindef.GeneratePluginCode
import net.maizegenetics.plugindef.PluginParameter
import net.maizegenetics.util.Utils
import org.apache.logging.log4j.LogManager
import java.awt.Frame
import javax.swing.ImageIcon

class CountAssociationsPlugin(parentFrame: Frame? = null, isInteractive: Boolean = false) : AbstractPlugin(parentFrame, isInteractive) {

    private val myLogger = LogManager.getLogger(CountAssociationsPlugin::class.java)

    private var intervalFile = PluginParameter.Builder("intervals", null, String::class.java)
            .description("Interval File.")
            .inFile()
            .required(true)
            .build()

    private var gwasResults = PluginParameter.Builder("gwasResults", null, String::class.java)
            .description("GWAS Result File.")
            .inFile()
            .required(true)
            .build()

    private var outputFile = PluginParameter.Builder("outputCountFile", null, String::class.java)
            .description("Output Count File.")
            .outFile()
            .required(true)
            .build()


    override fun processData(input: DataSet?): DataSet? {

        //Import the intervals into a rangeSet
        val intervalMap = loadInIntervalFile(intervalFile())

        //Map<IntervalName, Map<TraitName, count>>
        val outputMap = mutableMapOf<String, MutableMap<String, Int>>()

        val traitSet = mutableSetOf<String>()
        val intervalSet = intervalMap.asMapOfRanges().map { it.value }.toSortedSet()
        var traitCount = 0

        Utils.getBufferedReader(gwasResults()).use { reader ->
            var currLine: String? = reader.readLine()

            while (currLine != null) {
//                    if (currLine.startsWith("Trait")) {
                if (currLine.startsWith("trait")) {
                    currLine = reader.readLine()
                    continue
                }

                val currLineSplit = currLine.split(",")

                val traitName = currLineSplit[0]


                traitSet.add(traitName)
//                    val pos = Position.of(currLineSplit[2], currLineSplit[3].toInt())
                val pos = Position.of(currLineSplit[1], currLineSplit[2].toInt())
                val intervalHit = intervalMap.get(pos)


                if (intervalHit == null) {
                    currLine = reader.readLine()
                    continue
                }

                if (!outputMap.containsKey(intervalHit)) {
                    outputMap[intervalHit] = mutableMapOf()
                }
                if (!outputMap[intervalHit]!!.containsKey(traitName)) {
                    outputMap[intervalHit]?.put(traitName, 1)
                } else {
                    val currentCount = outputMap[intervalHit]?.get(traitName) ?: 0
                    outputMap[intervalHit]?.put(traitName, currentCount + 1)
                }


                currLine = reader.readLine()
            }
        }


        Utils.getBufferedWriter(outputFile()).use { output ->

            val sortedTraits = traitSet.sorted()
            output.write("Interval\t${sortedTraits.joinToString("\t")}\n")

            for (interval in intervalSet) {
                val countMap = outputMap[interval]
                val countList = mutableListOf<Int>()
                for (trait in sortedTraits) {
                    if (countMap == null) {
                        countList.add(0)
                    } else {
                        val count = countMap[trait] ?: 0
                        countList.add(count)
                    }
                }

                output.write("${interval}\t${countList.joinToString("\t")}\n")
            }
        }
        return null
    }

    fun loadInIntervalFile(intervalFile: String): RangeMap<Position, String> {
        val intervalMapBuilder = ImmutableRangeMap.Builder<Position, String>()
        Utils.getBufferedReader(intervalFile).readLines()
                .filter { !it.startsWith("seqnames") }
                .map { it.split(",") }
                .map { Pair(Range.closed(Position.of(it[0], it[1].toInt()), Position.of(it[0], it[2].toInt())), it[3]) }
                .forEach { intervalMapBuilder.put(it.first, it.second) }
        return intervalMapBuilder.build()
    }

    override fun getIcon(): ImageIcon? {
        return null
    }

    override fun getButtonName(): String {
        return ("Count Associations")
    }

    override fun getToolTipText(): String {
        return ("Create a count file For Associations")
    }

    /**
     * Interval File.
     *
     * @return Intervals
     */
    fun intervalFile(): String {
        return intervalFile.value()
    }

    /**
     * Set Intervals. Interval File.
     *
     * @param value Intervals
     *
     * @return this plugin
     */
    fun intervalFile(value: String): CountAssociationsPlugin {
        intervalFile = PluginParameter<String>(intervalFile, value)
        return this
    }

    /**
     * GWAS Result File.
     *
     * @return Gwas Results
     */
    fun gwasResults(): String {
        return gwasResults.value()
    }

    /**
     * Set Gwas Results. GWAS Result File.
     *
     * @param value Gwas Results
     *
     * @return this plugin
     */
    fun gwasResults(value: String): CountAssociationsPlugin {
        gwasResults = PluginParameter<String>(gwasResults, value)
        return this
    }

    /**
     * Output Count File.
     *
     * @return Output Count File
     */
    fun outputFile(): String {
        return outputFile.value()
    }

    /**
     * Set Output Count File. Output Count File.
     *
     * @param value Output Count File
     *
     * @return this plugin
     */
    fun outputFile(value: String): CountAssociationsPlugin {
        outputFile = PluginParameter<String>(outputFile, value)
        return this
    }


}

fun main(args: Array<String>) {
    GeneratePluginCode.generateKotlin(CountAssociationsPlugin::class.java)
}
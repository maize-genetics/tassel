package net.maizegenetics.analysis.data

import com.google.common.collect.ImmutableRangeMap
import com.google.common.collect.Range
import com.google.common.collect.RangeMap
import net.maizegenetics.dna.map.Position
import net.maizegenetics.plugindef.AbstractPlugin
import net.maizegenetics.plugindef.DataSet
import net.maizegenetics.plugindef.PluginParameter
import net.maizegenetics.util.Utils
import org.apache.logging.log4j.LogManager
import java.awt.Frame
import javax.swing.ImageIcon

class MeanR2FromLDPlugin(parentFrame: Frame? = null, isInteractive: Boolean = false) : AbstractPlugin(parentFrame, isInteractive) {

    private val myLogger = LogManager.getLogger(MeanR2FromLDPlugin::class.java)

    private var intervalFile = PluginParameter.Builder("intervals", null, String::class.java)
            .description("Interval File.")
            .inFile()
            .required(true)
            .build()

    private var ldResultFile = PluginParameter.Builder("ldResultFile", null, String::class.java)
            .description("LD Output File.")
            .inFile()
            .required(true)
            .build()

    private var output = PluginParameter.Builder("output", null, String::class.java)
            .description("Output Mean R2 Summary.")
            .outFile()
            .required(true)
            .build()

    override fun processData(input: DataSet?): DataSet? {

        val intervalMap = intervalMap(intervalFile())

        Utils.getBufferedReader(ldResultFile()).readLines()
            .drop(1)
            .map { it.split("\t") }
            .map { Pair(Position.of(it[0], it[1].toInt()), it[13].toDouble()) }
            .filter { !it.second.isNaN() }
                .forEach {
                    intervalMap[it.first]?.let { accumulateR2 ->
                        accumulateR2.r2Total += it.second
                        accumulateR2.r2Count++
                    }
                }

        Utils.getBufferedWriter(output()).use { writer ->
            writer.write("seqid,start,end,rr_id,average_r2_ld\n")
            intervalMap.asMapOfRanges()
                    .map { it.value }
                    .forEach {
                        writer.write(it.chr)
                        writer.write(",")
                        writer.write(it.start.toString())
                        writer.write(",")
                        writer.write(it.end.toString())
                        writer.write(",")
                        writer.write(it.id)
                        writer.write(",")
                        writer.write(if (it.r2Count == 0) "0" else (it.r2Total / it.r2Count.toDouble()).toString())
                        writer.write("\n")
                    }
        }

        return null
    }

    data class AccumulateR2(val id: String, val chr: String, val start: Int, val end: Int, var r2Total: Double = 0.0, var r2Count: Int = 0)

    fun intervalMap(intervalFile: String): RangeMap<Position, AccumulateR2> {
        val intervalMapBuilder = ImmutableRangeMap.Builder<Position, AccumulateR2>()
        Utils.getBufferedReader(intervalFile).readLines()
                .filter { !it.startsWith("seqnames") }
                .map { it.split(",") }
                .map { Pair(Range.closed(Position.of(it[0], it[1].toInt()), Position.of(it[0], it[2].toInt())), AccumulateR2(it[3], it[0], it[1].toInt(), it[2].toInt())) }
                .forEach { intervalMapBuilder.put(it.first, it.second) }
        return intervalMapBuilder.build()
    }

    override fun getIcon(): ImageIcon? {
        return null
    }

    override fun getButtonName(): String {
        return "Mean R2 From LD"
    }

    override fun getToolTipText(): String {
        return "Mean R2 From LD"
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
    fun intervalFile(value: String): MeanR2FromLDPlugin {
        intervalFile = PluginParameter(intervalFile, value)
        return this
    }

    /**
     * LD Output File.
     *
     * @return Ld Result File
     */
    fun ldResultFile(): String {
        return ldResultFile.value()
    }

    /**
     * Set Ld Result File. LD Output File.
     *
     * @param value Ld Result File
     *
     * @return this plugin
     */
    fun ldResultFile(value: String): MeanR2FromLDPlugin {
        ldResultFile = PluginParameter(ldResultFile, value)
        return this
    }

    /**
     * Output Mean R2 Summary.
     *
     * @return Output
     */
    fun output(): String {
        return output.value()
    }

    /**
     * Set Output. Output Mean R2 Summary.
     *
     * @param value Output
     *
     * @return this plugin
     */
    fun output(value: String): MeanR2FromLDPlugin {
        output = PluginParameter(output, value)
        return this
    }

}

/*
* ListPluginParameters
*/

package net.maizegenetics.analysis

import net.maizegenetics.plugindef.AbstractPlugin
import net.maizegenetics.plugindef.DataSet
import net.maizegenetics.plugindef.Plugin
import net.maizegenetics.plugindef.PluginParameter
import net.maizegenetics.util.TableReport
import net.maizegenetics.util.TableReportBuilder
import net.maizegenetics.util.Utils
import org.apache.logging.log4j.Level
import org.apache.logging.log4j.LogManager
import java.awt.Frame
import java.io.File
import java.util.*
import java.util.zip.ZipFile
import javax.swing.ImageIcon

/**
 *
 * @author Terry Casstevens
 */

private val myLogger = LogManager.getLogger(ListPluginParameters::class.java)

class ListPluginParameters(parentFrame: Frame? = null, isInteractive: Boolean = false) : AbstractPlugin(parentFrame, isInteractive) {

    private var jarFiles = PluginParameter.Builder<List<*>>("jarFiles", listOf("phg.jar"), List::class.java)
            .description("List of jar files.")
            .build()

    override fun processData(input: DataSet?): DataSet? {

        val temp: MutableList<String> = ArrayList()

        for (current in jarFiles()) {
            val classes = classes(current)
            for (current in classes) {
                if (Plugin.isPlugin(current)) {
                    temp.add(current)
                }
            }
        }

        val result = TableReportBuilder.getInstance("PluginParameters", arrayOf(Plugin.PARAMETER_PROPERTIES.Plugin, Plugin.PARAMETER_PROPERTIES.Parameter, Plugin.PARAMETER_PROPERTIES.Required, Plugin.PARAMETER_PROPERTIES.Default))
        temp.sort()
        for (current in temp) {
            val pluginName = Utils.getBasename(current)
            val plugin = Plugin.getPluginInstance(current, null)
            if (plugin != null) {
                val currentPluginParameters = (plugin as AbstractPlugin).usageParameters()
                for (parameterName in currentPluginParameters.keys) {
                    val parameterProperties = currentPluginParameters[parameterName]!!
                    result.add(arrayOf(pluginName, parameterName,
                            parameterProperties[Plugin.PARAMETER_PROPERTIES.Required],
                            parameterProperties[Plugin.PARAMETER_PROPERTIES.Default]))
                }
            }
        }

        return DataSet.getDataSet(result.build());

    }

    private fun classes(jarFile: String): Set<String> {
        val classpath = System.getProperty("java.class.path")
        val paths = classpath.split(File.pathSeparator.toRegex()).toTypedArray()
        var tasselPath: String? = null
        for (path in paths) {
            if (path.trim { it <= ' ' }.isNotEmpty()) {
                val file = File(path)
                if (file.exists()) {
                    tasselPath = file.absolutePath
                    if (tasselPath.endsWith(jarFile)) {
                        break
                    }
                }
            }
        }
        val classes: MutableSet<String> = LinkedHashSet()
        try {
            ZipFile(tasselPath).use { zFile ->
                val entries = zFile.entries()
                while (entries.hasMoreElements()) {
                    val entry = entries.nextElement()
                    if (!entry.isDirectory) {
                        var name = entry.name.replace(File.separator, ".")
                        if (name.endsWith(".class") && !name.contains("$")) {
                            name = name.substring(0, name.lastIndexOf(".class"))
                            classes.add(name)
                        }
                    }
                }
            }
        } catch (e: Exception) {
            e.printStackTrace()
        }
        return classes
    }

    /**
     * List of jar files.
     *
     * @return Name List
     */
    fun jarFiles(): List<String> {
        return jarFiles.value() as List<String>
    }

    /**
     * Set Name List. List of jar files.
     *
     * @param value Name List
     *
     * @return this plugin
     */
    fun jarFiles(value: List<String>): ListPluginParameters {
        jarFiles = PluginParameter<List<*>>(jarFiles, value)
        return this
    }

    override fun getIcon(): ImageIcon? {
        return null
    }

    override fun getButtonName(): String {
        return "List Plugins"
    }

    override fun getToolTipText(): String {
        return "List Plugins"
    }

    init {
        LogManager.getLogger("net.maizegenetics").atLevel(Level.OFF)
        LogManager.getLogger("net.maizegenetics.plugindef").atLevel(Level.INFO)
    }

    companion object {
        @JvmStatic
        fun forPHG() = ListPluginParameters().performFunction(null).getData(0).data as TableReport

        @JvmStatic
        fun forTASSEL() = ListPluginParameters().jarFiles(arrayListOf("sTASSEL.jar")).performFunction(null).getData(0).data as TableReport
    }

}

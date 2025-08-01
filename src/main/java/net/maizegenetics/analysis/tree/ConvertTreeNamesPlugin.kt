package net.maizegenetics.analysis.tree

import net.maizegenetics.plugindef.AbstractPlugin
import net.maizegenetics.plugindef.DataSet
import net.maizegenetics.plugindef.Datum
import net.maizegenetics.plugindef.PluginParameter
import net.maizegenetics.taxa.tree.Tree
import net.maizegenetics.taxa.tree.convertNames
import org.apache.logging.log4j.LogManager
import java.awt.Frame
import javax.swing.ImageIcon

private val myLogger = LogManager.getLogger(ConvertTreeNamesPlugin::class.java)

class ConvertTreeNamesPlugin(parentFrame: Frame? = null, isInteractive: Boolean = false) : AbstractPlugin(parentFrame, isInteractive) {

    private var filename = PluginParameter.Builder<String>("nameList", null, String::class.java)
            .description("Conversion filename. Each line of the file should have two names separated with a tab. The first name should be converted to the second.")
            .inFile()
            .required(true)
            .build()

    override fun preProcessParameters(input: DataSet?) {
        val temp = input?.getDataOfType(Tree::class.java)
        val numTrees = temp?.size ?: 0
        if (numTrees != 1) {
            throw IllegalArgumentException("ConvertTreeNamesPlugin: processData: must input 1 tree")
        }
    }

    override fun processData(input: DataSet?): DataSet? {
        val tree = input!!.getDataOfType(Tree::class.java).get(0).data as Tree
        return DataSet(Datum("Convert Tree Names", convertNames(tree, filename()), null), this)
    }

    /**
     * Conversion filename. Each line of the file should have
     * two names separated with a tab. The first name should
     * be converted to the second.
     *
     * @return Name List
     */
    fun filename(): String {
        return filename.value()
    }

    /**
     * Set Name List. Conversion filename. Each line of the
     * file should have two names separated with a tab. The
     * first name should be converted to the second.
     *
     * @param value Name List
     *
     * @return this plugin
     */
    fun filename(value: String): ConvertTreeNamesPlugin {
        filename = PluginParameter<String>(filename, value)
        return this
    }

    override fun getToolTipText(): String {
        return "Convert Tree Names"
    }

    override fun getButtonName(): String {
        return "Convert Tree Names"
    }

    override fun getIcon(): ImageIcon? {
        return null
    }

}
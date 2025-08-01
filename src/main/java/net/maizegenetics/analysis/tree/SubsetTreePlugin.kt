package net.maizegenetics.analysis.tree

import net.maizegenetics.analysis.filter.FilterTaxaBuilderPlugin
import net.maizegenetics.plugindef.AbstractPlugin
import net.maizegenetics.plugindef.DataSet
import net.maizegenetics.plugindef.Datum
import net.maizegenetics.plugindef.PluginParameter
import net.maizegenetics.taxa.TaxaList
import net.maizegenetics.taxa.tree.Tree
import net.maizegenetics.taxa.tree.subsetTree
import org.apache.logging.log4j.LogManager
import java.awt.Frame
import javax.swing.ImageIcon

private val myLogger = LogManager.getLogger(SubsetTreePlugin::class.java)

class SubsetTreePlugin(parentFrame: Frame? = null, isInteractive: Boolean = false) : AbstractPlugin(parentFrame, isInteractive) {

    private var nameList = PluginParameter.Builder<TaxaList>("nameList", null, TaxaList::class.java)
            .description("List of names to include. This can be a comma separated list of names (no spaces unless surrounded by quotes), file (.txt) with list of names to include, or a taxa list file (.json or .json.gz).")
            .required(true)
            .build()

    override fun preProcessParameters(input: DataSet?) {
        val temp = input?.getDataOfType(Tree::class.java)
        val numTrees = temp?.size ?: 0
        if (numTrees != 1) {
            throw IllegalArgumentException("SubsetTreePlugin: processData: must input 1 tree")
        }
    }

    override fun processData(input: DataSet?): DataSet? {
        val tree = input!!.getDataOfType(Tree::class.java).get(0).data as Tree
        return DataSet(Datum("Subset Tree", subsetTree(tree, nameList()), null), this)
    }

    /**
     * List of names to include. This can be a comma separated
     * list of names (no spaces unless surrounded by quotes),
     * file (.txt) with list of names to include, or a taxa
     * list file (.json or .json.gz).
     *
     * @return Name List
     */
    fun nameList(): TaxaList {
        return nameList.value()
    }

    /**
     * Set Name List. List of names to include. This can be
     * a comma separated list of names (no spaces unless surrounded
     * by quotes), file (.txt) with list of names to include,
     * or a taxa list file (.json or .json.gz).
     *
     * @param value Name List
     *
     * @return this plugin
     */
    fun nameList(value: TaxaList): SubsetTreePlugin {
        nameList = PluginParameter<TaxaList>(nameList, value)
        return this
    }

    override fun getToolTipText(): String {
        return "Subset Tree"
    }

    override fun getButtonName(): String {
        return "Subset Tree"
    }

    override fun getIcon(): ImageIcon? {
        val imageURL = FilterTaxaBuilderPlugin::class.java.getResource("/net/maizegenetics/analysis/images/FilterNew.gif")
        return imageURL?.let { ImageIcon(it) }
    }

}
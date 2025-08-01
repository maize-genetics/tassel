package net.maizegenetics.analysis.tree

import net.maizegenetics.plugindef.AbstractPlugin
import net.maizegenetics.plugindef.DataSet
import net.maizegenetics.plugindef.Datum
import net.maizegenetics.taxa.tree.Tree
import net.maizegenetics.taxa.tree.mergeTrees
import org.apache.logging.log4j.LogManager
import java.awt.Frame
import javax.swing.ImageIcon

private val myLogger = LogManager.getLogger(MergeTreesPlugin::class.java)

class MergeTreesPlugin(parentFrame: Frame? = null, isInteractive: Boolean = false) : AbstractPlugin(parentFrame, isInteractive) {

    override fun processData(input: DataSet?): DataSet {

        val temp = input?.getDataOfType(Tree::class.java)
        val numTrees = temp?.size ?: 0
        if (numTrees < 2) {
            throw IllegalArgumentException("MergeTreesPlugin: processData: must input at least 2 trees")
        }
        val trees = temp!!.map { it.data as Tree }

        return DataSet(Datum("Merged Tree", mergeTrees(trees), null), this)

    }

    override fun getToolTipText(): String {
        return "Merge Trees"
    }

    override fun getButtonName(): String {
        return "Merge Trees"
    }

    override fun getIcon(): ImageIcon? {
        val imageURL = MergeTreesPlugin::class.java.getResource("/net/maizegenetics/analysis/images/Merge.gif")
        return imageURL?.let { ImageIcon(it) }
    }

}
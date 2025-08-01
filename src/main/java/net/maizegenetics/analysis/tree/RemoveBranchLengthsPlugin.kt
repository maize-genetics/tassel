package net.maizegenetics.analysis.tree

import net.maizegenetics.plugindef.AbstractPlugin
import net.maizegenetics.plugindef.DataSet
import net.maizegenetics.plugindef.Datum
import net.maizegenetics.taxa.tree.Tree
import net.maizegenetics.taxa.tree.removeBranchLengths
import org.apache.logging.log4j.LogManager
import java.awt.Frame
import javax.swing.ImageIcon

private val myLogger = LogManager.getLogger(RemoveBranchLengthsPlugin::class.java)

class RemoveBranchLengthsPlugin(parentFrame: Frame? = null, isInteractive: Boolean = false) : AbstractPlugin(parentFrame, isInteractive) {

    override fun processData(input: DataSet?): DataSet {

        val temp = input?.getDataOfType(Tree::class.java)
        val numTrees = temp?.size ?: 0
        if (numTrees != 1) {
            throw IllegalArgumentException("RemoveBranchLengthsPlugin: processData: must input 1 tree")
        }
        val tree = temp!![0].data as Tree

        return DataSet(Datum("Remove Tree Branch Lengths", removeBranchLengths(tree), null), this)

    }

    override fun getToolTipText(): String {
        return "Remove Tree Branch Lengths"
    }

    override fun getButtonName(): String {
        return "Remove Tree Branch Lengths"
    }

    override fun getIcon(): ImageIcon? {
        return null
    }

}
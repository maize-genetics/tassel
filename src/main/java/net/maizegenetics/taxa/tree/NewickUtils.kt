@file:JvmName("NewickUtils")

package net.maizegenetics.taxa.tree

import net.maizegenetics.taxa.TaxaList
import org.apache.commons.lang3.StringUtils
import org.apache.logging.log4j.LogManager
import java.io.BufferedWriter
import java.io.File
import java.lang.Double
import java.util.*
import kotlin.collections.HashMap

/**
 * These utilities are related to the Newick Tree Format.
 * http://evolution.genetics.washington.edu/phylip/newicktree.html
 */

private val myLogger = LogManager.getLogger("net.maizegenetics.taxa.tree.NewickUtils")

/**
 * Creates a Tree from the given newick formatted file
 */
fun read(filename: String): Tree {

    val newick = File(filename).readLines().joinToString(separator = "") { it }

    if (newick.count { it == '(' } != newick.count { it == ')' }) {
        throw IllegalArgumentException("NewickUtils: read: $filename: number of open parenthesis doesn't match number of close parenthesis.")
    }

    if (newick.last() != ';') {
        throw IllegalArgumentException("NewickUtils: read: $filename: doesn't end with semicolon.")
    }

    return SimpleTree(makeNode(newick.substring(0, newick.length - 1)))

}

/**
 * Creates a node in the tree recursively creates the children for that node.
 */
private fun makeNode(newick: String): Node {

    val nameBranchLength = newick.substringAfterLast(')')

    // The node name comes before the :. If no : present, this is the name
    val name = nameBranchLength.substringBefore(':').replace("'", "").let {
        try {
            // Names that are numbers are changed to empty string
            // These are support values from 0-100, that can be ignored for our purposes
            Double.parseDouble(it)
            ""
        } catch (ne: NumberFormatException) {
            // Newick format considers underscores to be spaces
            it.replace("_", " ")
        }
    }
    val branchLength = if (nameBranchLength.contains(':')) nameBranchLength.substringAfter(':').toDouble() else 0.0

    val result = SimpleNode(name, branchLength)

    if (nameBranchLength != newick) {

        var parenthesisCount = 0
        var currentStr = StringBuilder()
        newick.substringAfter('(').substringBeforeLast(')').forEach {
            when (it) {
                '(' -> {
                    parenthesisCount++
                    currentStr.append(it)
                }
                ')' -> {
                    parenthesisCount--
                    currentStr.append(it)
                }
                ',' -> {
                    if (parenthesisCount == 0) {
                        result.addChild(makeNode(currentStr.toString()))
                        currentStr = StringBuilder()
                    } else {
                        currentStr.append(it)
                    }
                }
                else -> currentStr.append(it)
            }
        }
        if (currentStr.isNotEmpty()) result.addChild(makeNode(currentStr.toString()))

    }

    return result

}

/**
 * Writes give tree to Newick formatted file.
 * http://evolution.genetics.washington.edu/phylip/newicktree.html
 */
fun write(filename: String, tree: Tree, includeBranchLengths: Boolean = true) {

    try {
        File(filename).bufferedWriter().use { writer ->
            write(tree.root, writer, includeBranchLengths)
            writer.append(";\n")
        }
    } catch (e: Exception) {
        myLogger.debug(e.message, e)
        throw IllegalStateException("NewickUtils: write: problem writing: $filename.\n${e.message}")
    }

}

/**
 * Recursively writes nodes and its children to a file
 */
private fun write(node: Node, writer: BufferedWriter, includeBranchLengths: Boolean) {

    if (!node.isLeaf) {
        writer.append("(")
        for (i in 0 until node.childCount) {
            if (i != 0) writer.append(",")
            write(node.getChild(i), writer, includeBranchLengths)
        }
        writer.append(")")
    }

    node.identifier?.name?.let {
        if (it.isNotEmpty()) {
            val name = it.replace(' ', '_')
            writer.append("'")
            writer.append(name)
            writer.append("'")
        }
    }

    if (includeBranchLengths && node.branchLength != 0.0) {
        writer.append(":")
        var lengthStr = "%.7f".format(node.branchLength)
        lengthStr = StringUtils.stripEnd(lengthStr, "0")
        lengthStr = StringUtils.stripEnd(lengthStr, ".")
        writer.append(lengthStr)
    }

}

private const val MERGE_ROOT_NODE = "MERGE_ROOT_NODE"

/**
 * This merges two or more trees into one tree.
 * Nodes are merged when the names are the same.
 *
 * Example:
 * (AA:2,BB:2,CC:2);
 * merged with
 * (XX:2,YY:2,ZZ:2)AA;
 * produces
 * (('XX':2.0000000,'YY':2.0000000,'ZZ':2.0000000)'AA':2.0000000,'BB':2.0000000,'CC':2.0000000);
 */
fun mergeTrees(trees: List<Tree>): Tree {

    if (trees.size < 2) {
        throw IllegalArgumentException("NewickUtils: mergeTrees: must supply at least 2 trees.")
    }

    val nameToNode = HashMap<String, Node>()
    val nodeToNode = IdentityHashMap<Node, Node>()

    lateinit var rootNode: Node

    trees.forEach { tree ->

        val nodes = tree.nodes()

        nodes.forEach { node ->

            val nodeName = node.identifier.name
            val existingNode = nameToNode[nodeName]
            val existingParentNode = nodeToNode[node.parent]
            val existingParentName = existingParentNode?.identifier?.name
            val parentName = node.parent?.identifier?.name

            if (nodeToNode.isEmpty()) {
                rootNode = SimpleNode(nodeName, node.branchLength)
                nodeToNode[node] = rootNode
                if (nodeName.isNotEmpty()) nameToNode[nodeName] = rootNode!!
            } else if (existingNode != null) {
                if (existingParentName != null && parentName != null && existingParentName != parentName) {
                    throw IllegalArgumentException("NewickUtils: mergeTrees: node: $nodeName has different parents ($parentName, ${existingParentName}) in different trees.")
                }
                if (existingNode.branchLength != node.branchLength) {
                    myLogger.warn("mergeTrees: nodes named: $nodeName in different trees have different branch lengths.")
                }
                nodeToNode[node] = existingNode
            } else if (existingParentNode == null) {
                val mergeRoot = if (nameToNode[MERGE_ROOT_NODE] != null) {
                    nameToNode[MERGE_ROOT_NODE]!!
                } else {
                    val temp = SimpleNode("", 0.0)
                    nameToNode[MERGE_ROOT_NODE] = temp
                    temp
                }
                if (rootNode != mergeRoot) mergeRoot.addChild(rootNode)
                val newChild = SimpleNode(nodeName, node.branchLength)
                nodeToNode[node] = newChild
                if (nodeName.isNotEmpty()) nameToNode[nodeName] = newChild
                mergeRoot.addChild(newChild)

                rootNode = mergeRoot
            } else {
                val newChild = SimpleNode(nodeName, node.branchLength)
                nodeToNode[node] = newChild
                if (nodeName.isNotEmpty()) nameToNode[nodeName] = newChild
                existingParentNode.addChild(newChild)
            }

        }

    }

    return SimpleTree(rootNode)

}

/**
 * Extension function that returns all nodes of the tree
 */
fun Tree.nodes(): List<Node> {

    val result = mutableListOf<Node>()
    result.add(this.root)
    addChildren(result, this.root)
    return result

}

/**
 * Recursively adds children to the list
 */
private fun addChildren(nodes: MutableList<Node>, node: Node) {

    for (i in 0 until node.childCount) {
        nodes.add(node.getChild(i))
        addChildren(nodes, node.getChild(i))
    }

}

/**
 * Returns a tree that is a subset of the
 * original containing only the taxa specified.
 */
fun subsetTree(tree: Tree, taxaList: TaxaList): Tree {
    val nameList = taxaList.map { it.name }
    return subsetTree(tree, nameList)
}

/**
 * Returns a tree that is a subset of the
 * original containing only the names specified.
 */
fun subsetTree(tree: Tree, namesToKeep: List<String>): Tree {
    return SimpleTree(keepNode(tree.root, namesToKeep))
}

/**
 * Recursively iterates over the tree looking for nodes
 * that are specified by the list of names.
 */
private fun keepNode(node: Node, namesToKeep: List<String>): Node? {

    if (node.isLeaf) {
        node.identifier?.name?.let {
            if (namesToKeep.contains(it)) {
                return SimpleNode(node.identifier?.name, node.branchLength)
            }
        }
    } else {

        val childNodesToKeep = mutableListOf<Node>()
        for (i in 0 until node.childCount) {
            keepNode(node.getChild(i), namesToKeep)?.let {
                childNodesToKeep.add(it)
            }
        }

        when (childNodesToKeep.size) {
            0 -> return null
            1 -> {
                childNodesToKeep[0].branchLength += node.branchLength
                return childNodesToKeep[0]
            }
            else -> {
                val parent = SimpleNode(node.identifier?.name, node.branchLength)
                childNodesToKeep.forEach {
                    parent.addChild(it)
                }
                return parent
            }
        }

    }

    return null

}

/**
 * This removes the branch lengths from
 * the given tree
 */
fun removeBranchLengths(tree: Tree): Tree {
    return SimpleTree(removeBranchLengths(tree.root))
}

/**
 * Recursively iterates over the tree removing
 * branch lengths.
 */
private fun removeBranchLengths(node: Node): Node {

    val result = SimpleNode(node.identifier?.name, 0.0)

    for (i in 0 until node.childCount) {
        result.addChild(removeBranchLengths(node.getChild(i)))
    }

    return result

}

/**
 * This converts tree node names based on the
 * conversion file specified.
 */
fun convertNames(tree: Tree, filename: String): Tree {

    val conversions = mutableMapOf<String, String>()

    File(filename).bufferedReader().lines().forEach {
        val temp = it.split("\t")
        if (temp.size != 2) {
            throw IllegalArgumentException("NewickUtils: convertNames: each line of file: $filename should have two names separated by a tab.\nOffending line: $it")
        }
        conversions[temp[0].trim().replace('_', ' ')] = temp[1].trim().replace('_', ' ')
    }

    return SimpleTree(convertNames(tree.root, conversions))

}

/**
 * Recursively iterates over the tree changing node names
 * that are specified in the map
 */
private fun convertNames(node: Node, conversions: Map<String, String>): Node {

    val name = conversions[node.identifier?.name] ?: node.identifier?.name

    val result = SimpleNode(name, 0.0)

    for (i in 0 until node.childCount) {
        result.addChild(convertNames(node.getChild(i), conversions))
    }

    return result

}

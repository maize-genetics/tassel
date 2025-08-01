@file:JvmName("TasselVersions")

package net.maizegenetics.tassel

import net.maizegenetics.plugindef.AbstractPlugin
import net.maizegenetics.util.Utils
import org.apache.logging.log4j.LogManager
import org.w3c.dom.Node
import javax.xml.parsers.DocumentBuilderFactory

/**
 * This class is used to get the version information for Tassel and third party libraries.
 * The version information is stored in a file called tassel_info.xml in the resources directory.
 *
 * The file is in the following format:
 *
 * <?xml version="1.0" encoding="UTF-8" standalone="no"?>
 * <Tassel>
 *     <TasselLibrary>
 *         <key>PHG</key>
 *         <name>Practical Haplotype Graph (PHG)</name>
 *         <version>1.3</version>
 *         <date>January 10, 2023</date>
 *         <citation>Bradbury PJ, Casstevens T, Jensen SE, Johnson LC, Miller ZR, Monier B, Romay MC, Song B, Buckler ES. The Practical Haplotype Graph, a platform for storing and using pangenomes for imputation. Bioinformatics. 2022 Aug 2;38(15):3698-3702. doi: 10.1093/bioinformatics/btac410. PMID: 35748708; PMCID: PMC9344836.</citation>
 *     </TasselLibrary>
 * </Tassel>
 *
 */
object TasselVersions {

    private val myLogger = LogManager.getLogger(LibraryInfo::class.java)

    private val infoMap = mutableMapOf<String, LibraryInfo>()

    data class LibraryInfo(val name: String, val version: String, val date: String, val citation: String?)

    init {
        for (current in Utils.getFullyQualifiedResourceNames("tassel_info.xml")) {
            getThirdPartyLibraryInfo(current)
        }
    }

    fun tasselName() = ""

    fun tasselVersion() = TASSELMainFrame.version

    fun tasselVersionDate() = TASSELMainFrame.versionDate

    fun tasselCitation() = AbstractPlugin.DEFAULT_CITATION

    fun phgName() = infoMap["PHG"]!!.name

    fun phgVersion() = infoMap["PHG"]!!.version

    fun phgVersionDate() = infoMap["PHG"]!!.date

    fun phgCitation() = infoMap["PHG"]!!.citation

    fun libraryInfo(library: String): LibraryInfo? = infoMap[library]

    fun libraryInfos() = infoMap.entries

    private fun getThirdPartyLibraryInfo(filename: String) {
        try {
            LibraryInfo::class.java.getResourceAsStream(filename).use { input ->
                val dbFactory = DocumentBuilderFactory.newInstance()
                val dBuilder = dbFactory.newDocumentBuilder()
                val doc = dBuilder.parse(input)
                doc.documentElement.normalize()
                val rootElement = doc.documentElement
                require(
                    rootElement.nodeName.equals(
                        "Tassel",
                        ignoreCase = true
                    )
                ) { """LibraryInfo: getThirdPartyLibraryInfo: Root Node must be Tassel: ${rootElement.nodeName}""" }
                val children = rootElement.childNodes
                for (i in 0 until children.length) {
                    val current = children.item(i)
                    addLibInfo(current)
                }
            }
        } catch (e: Exception) {
            myLogger.debug(e.message, e)
            throw IllegalStateException(
                """
                LibraryInfo: getThirdPartyLibraryInfo: Problem reading XML file: $filename
                ${e.message}
                """.trimIndent()
            )
        }
    }

    private fun addLibInfo(rootElement: Node) {

        try {
            if (!rootElement.nodeName.trim { it <= ' ' }.equals("TasselLibrary", ignoreCase = true)) {
                return
            }
            val children = rootElement.childNodes
            val libraryInfo = mutableMapOf<String, String>()
            for (i in 0 until children.length) {
                val current = children.item(i)
                val elementName = current.nodeName.trim()
                if (current.nodeType == Node.ELEMENT_NODE) {
                    libraryInfo[elementName] = current.textContent.trim()
                }
            }
            infoMap[libraryInfo["key"]!!] =
                LibraryInfo(
                    libraryInfo["name"]!!,
                    libraryInfo["version"]!!,
                    libraryInfo["date"]!!,
                    libraryInfo["citation"]!!
                )
        } catch (e: java.lang.Exception) {
            myLogger.debug(e.message, e)
        }

    }

}
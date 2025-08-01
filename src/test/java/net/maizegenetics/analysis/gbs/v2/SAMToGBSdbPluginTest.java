package net.maizegenetics.analysis.gbs.v2;

import static org.junit.Assert.*;

import java.lang.reflect.Method;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Optional;

import junit.framework.Assert;
import net.maizegenetics.dna.map.Chromosome;
import net.maizegenetics.dna.map.GeneralPosition;
import net.maizegenetics.dna.map.Position;
import net.maizegenetics.dna.tag.SAMUtils;
import net.maizegenetics.dna.tag.Tag;
import net.maizegenetics.dna.tag.TagBuilder;
import net.maizegenetics.util.Tuple;
import org.junit.Ignore;

import org.junit.Test;

/**
 * Defines xxxx
 *
 * @author Ed Buckler
 */
public class SAMToGBSdbPluginTest {

    // bowtie2 default cases (with or without the -M option)
    String[] bowtie2 = {
        // forward strand perfect match:
        "length=64count=309	0	chr10	9373084	11	64M	*	0	0	CAGCAAAATAAATTACAAAACAAAGTGTGATCTACAGGTACCTGTAAGCCTGCTATCTGCTACG	ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff	AS:i:128	XS:i:120	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:64	YT:Z:UU",
        // reverse strand, 51 bp, perfect match (13433871+51M-1 = 13433921):
        "length=51count=974	16	chr9	13433871	32	51M	*	0	0	CAGCAGACATTAGGAGGCTATGGCCTTAGCACTGGGAGTGCTTGTTTGCTG	fffffffffffffffffffffffffffffffffffffffffffffffffff	AS:i:102	XS:i:86	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:51	YT:Z:UU",
        // reverse strand, soft clipping by 2 bases on both ends, one mismatch (19773189+60M+2S-1 = 19773250):
        "length=64count=14	16	chr10	19773189	18	2S60M2S	*	0	0	GCTGTCTCCCAGGACGGCAAAAAGAAGAAGAAGAACCACGACTATGAGAGGACACAGTTTGCTG	ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff	AS:i:112	XS:i:80	XN:i:0	XM:i:1	XO:i:0	XG:i:0	NM:i:1	MD:Z:12T47	YT:Z:UU",
        // multiple alighments of equal score:
        "length=62count=29	0	2	219523831	1	62M	*	0	0	CAGCAAAAAAATCTGAGAGCGGCAAGCGGCAAGGAGATGAAAGATCTCCAGTTAAGCGGCAG	ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff	AS:i:116	XS:i:116	XN:i:0	XM:i:1	XO:i:0	XG:i:0	NM:i:1	MD:Z:18A43	YT:Z:UU",
        // no alignment:
        "length=64count=3137	4	*	0	0	*	*	0	0	CAGCAAAAAAATCCTTGCCCATGTAGGAGAACCGATTTGGCTCTCCCACATGAGCAAGGTCTTG	ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff	YT:Z:UU",
        // no alignment, but reverse complemented anyway:
        "length=64count=3137	20	*	0	0	*	*	0	0	CAGCAAAAAAATCCTTGCCCATGTAGGAGAACCGATTTGGCTCTCCCACATGAGCAAGGTCTTG	ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff	YT:Z:UU",
        // only one good alignment (no "XS:i:#" entry):
        "length=64count=109	0	1	224071175	44	64M	*	0	0	CAGCAAAAAAATCTAGATCATCTACCGACAGCCTTTCATGGAATACAATAATTTGTCTGTCTTC	ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff	AS:i:128	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:64	YT:Z:UU",
//        // ya got any more?:
//        "",
    };
    Tag[] bowtie2ExpTags = {
        TagBuilder.instance("CAGCAAAATAAATTACAAAACAAAGTGTGATCTACAGGTACCTGTAAGCCTGCTATCTGCTACG").build(),
        TagBuilder.instance("CAGCAAACAAGCACTCCCAGTGCTAAGGCCATAGCCTCCTAATGTCTGCTG").build(),
        TagBuilder.instance("CAGCAAACTGTGTCCTCTCATAGTCGTGGTTCTTCTTCTTCTTTTTGCCGTCCTGGGAGACAGC").build(),
        TagBuilder.instance("CAGCAAAAAAATCTGAGAGCGGCAAGCGGCAAGGAGATGAAAGATCTCCAGTTAAGCGGCAG").build(),
        TagBuilder.instance("CAGCAAAAAAATCCTTGCCCATGTAGGAGAACCGATTTGGCTCTCCCACATGAGCAAGGTCTTG").build(),
        TagBuilder.instance("CAAGACCTTGCTCATGTGGGAGAGCCAAATCGGTTCTCCTACATGGGCAAGGATTTTTTTGCTG").build(), //reverse complement of previous
        TagBuilder.instance("CAGCAAAAAAATCTAGATCATCTACCGACAGCCTTTCATGGAATACAATAATTTGTCTGTCTTC").build(),
//        TagBuilder.instance("").build(),
    };
    Position[] bowtie2ExpPositions = {
        new GeneralPosition.Builder(new Chromosome("10"),9373084).addAnno("forward","true").addAnno("cigar","64M").addAnno("supportvalue","128").addAnno("mappingapproach","Bowtie2").strand((byte)1).build(),
        new GeneralPosition.Builder(new Chromosome("9"),13433871+51-1).addAnno("forward","false").addAnno("cigar","51M").addAnno("supportvalue","102").addAnno("mappingapproach","Bowtie2").strand((byte)1).build(),
        new GeneralPosition.Builder(new Chromosome("10"),19773189+60+2-1).addAnno("forward","false").addAnno("cigar","2S60M2S").addAnno("supportvalue","112").addAnno("mappingapproach","Bowtie2").strand((byte)1).build(),
        // this next one (multiple alighments of equal score) probably should be null (so it is commented out and replaced by null):
//        new GeneralPosition.Builder(new Chromosome("2"),219523831).addAnno("forward","true").addAnno("cigar","62M").addAnno("supportvalue","116").addAnno("mappingapproach","Bowtie2").strand((byte)1).build(),
        null,
        null,
        null,
        new GeneralPosition.Builder(new Chromosome("1"),-1).addAnno("forward","true").addAnno("cigar","64M").addAnno("supportvalue","128").addAnno("mappingapproach","Bowtie2").strand((byte)1).build(),
//        new GeneralPosition.Builder(new Chromosome(""),-1).addAnno("forward","").addAnno("cigar","").addAnno("supportvalue","").addAnno("mappingapproach","Bowtie2").strand((byte)1).build(),
    };

    // bowtie2 -k 5 cases
    String[] bowtie2k5 = {
        // five alternative alignments (with unique best primary alignment on forward strand):
        "length=64count=1731	0	chr4	34593346	11	64M	*	0	0	CAGCCGGACTAGCAAAACAAGTGCGTAAAATCGAAGTCTTCCTCAAAGGACCATCCATGCAAAA	ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff	AS:i:96	XS:i:88	XN:i:0	XM:i:4	XO:i:0	XG:i:0	NM:i:4	MD:Z:25C21A3T4G7	YT:Z:UU",
        "length=64count=1731	256	chr4	186122318	11	44M1D20M	*	0	0	CAGCCGGACTAGCAAAACAAGTGCGTAAAATCGAAGTCTTCCTCAAAGGACCATCCATGCAAAA	ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff	AS:i:88	XS:i:88	XN:i:0	XM:i:4	XO:i:1	XG:i:1	NM:i:5	MD:Z:25A18^A3A3T4G7	YT:Z:UU",
        "length=64count=1731	272	chr6	16023721	11	16M1D48M	*	0	0	TTTTGCATGGATGGTCCTTTGAGGAAGACTTCGATTTTACGCACTTGTTTTGCTAGTCCGGCTG	ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff	AS:i:80	XS:i:88	XN:i:0	XM:i:5	XO:i:1	XG:i:1	NM:i:6	MD:Z:7C4A3^T0T6A14T25	YT:Z:UU",
        "length=64count=1731	272	chr6	9678484	11	15M1I48M	*	0	0	TTTTGCATGGATGGTCCTTTGAGGAAGACTTCGATTTTACGCACTTGTTTTGCTAGTCCGGCTG	ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff	AS:i:78	XS:i:88	XN:i:0	XM:i:5	XO:i:1	XG:i:1	NM:i:6	MD:Z:7C4A5C1C16T25	YT:Z:UU",
        "length=64count=1731	272	chr8	115376743	11	25S7M1D32M	*	0	0	TTTTGCATGGATGGTCCTTTGAGGAAGACTTCGATTTTACGCACTTGTTTTGCTAGTCCGGCTG	ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff	AS:i:54	XS:i:88	XN:i:0	XM:i:2	XO:i:1	XG:i:1	NM:i:3	MD:Z:7^G1T4T25	YT:Z:UU",
    };
    Tag[] bowtie2k5ExpTags = {
        TagBuilder.instance("CAGCCGGACTAGCAAAACAAGTGCGTAAAATCGAAGTCTTCCTCAAAGGACCATCCATGCAAAA").build(),
        TagBuilder.instance("CAGCCGGACTAGCAAAACAAGTGCGTAAAATCGAAGTCTTCCTCAAAGGACCATCCATGCAAAA").build(),
        TagBuilder.instance("CAGCCGGACTAGCAAAACAAGTGCGTAAAATCGAAGTCTTCCTCAAAGGACCATCCATGCAAAA").build(),
        TagBuilder.instance("CAGCCGGACTAGCAAAACAAGTGCGTAAAATCGAAGTCTTCCTCAAAGGACCATCCATGCAAAA").build(),
        TagBuilder.instance("CAGCCGGACTAGCAAAACAAGTGCGTAAAATCGAAGTCTTCCTCAAAGGACCATCCATGCAAAA").build(),
    };
    Position[] bowtie2k5ExpPositions = {
        new GeneralPosition.Builder(new Chromosome("4"),34593346).addAnno("forward","true").addAnno("cigar","64M").addAnno("supportvalue","96").addAnno("mappingapproach","Bowtie2").strand((byte)1).build(),
        new GeneralPosition.Builder(new Chromosome("4"),186122318).addAnno("forward","true").addAnno("cigar","44M1D20M").addAnno("supportvalue","88").addAnno("mappingapproach","Bowtie2").strand((byte)1).build(),
        new GeneralPosition.Builder(new Chromosome("6"),16023721).addAnno("forward","false").addAnno("cigar","16M1D48M").addAnno("supportvalue","80").addAnno("mappingapproach","Bowtie2").strand((byte)1).build(),
        new GeneralPosition.Builder(new Chromosome("6"),9678484).addAnno("forward","false").addAnno("cigar","15M1I48M").addAnno("supportvalue","78").addAnno("mappingapproach","Bowtie2").strand((byte)1).build(),
        new GeneralPosition.Builder(new Chromosome("8"),115376743).addAnno("forward","false").addAnno("cigar","25S7M1D32M").addAnno("supportvalue","54").addAnno("mappingapproach","Bowtie2").strand((byte)1).build(),
//        new GeneralPosition.Builder(new Chromosome(""),-1).addAnno("forward","").addAnno("cigar","").addAnno("supportvalue","").addAnno("mappingapproach","Bowtie2").strand((byte)1).build(),
    };

    // bwa-mem -a cases
    String[] bwaMemA = {
        // same tag as the first bowtie2 -k 5 case above:
        "length=64count=1731	0	chr4	34593346	0	64M	*	0	0	CAGCCGGACTAGCAAAACAAGTGCGTAAAATCGAAGTCTTCCTCAAAGGACCATCCATGCAAAA	ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff	NM:i:4	MD:Z:25C21A3T4G7	AS:i:44	XS:i:44",
        "length=64count=1731	256	chr8	15458003	0	64M	*	0	0	*	*	NM:i:4	MD:Z:14C32A3T4G7	AS:i:44",
        "length=64count=1731	256	chr1	264859139	0	64M	*	0	0	*	*	NM:i:4	MD:Z:14C32A3T4G7	AS:i:44",
        "length=64count=1731	256	chr10	91473731	0	64M	*	0	0	*	*	NM:i:4	MD:Z:14C32A3T4G7	AS:i:44",
        "length=64count=1731	256	chr7	62589280	0	64M	*	0	0	*	*	NM:i:4	MD:Z:14C32A3T4G7	AS:i:44",
        "length=64count=1731	256	chr7	62589358	0	64M	*	0	0	*	*	NM:i:4	MD:Z:14C32A3T4G7	AS:i:44",
        "length=64count=1731	272	chr7	96208663	0	64M	*	0	0	*	*	NM:i:4	MD:Z:7C4A3T32G14	AS:i:44",
        "length=64count=1731	272	chr1	252153050	0	17H47M	*	0	0	*	*	NM:i:1	MD:Z:21T25	AS:i:42",
        "length=64count=1731	256	chr4	186122318	0	47M17H	*	0	0	*	*	NM:i:1	MD:Z:25A21	AS:i:42",
        "length=64count=1731	256	chr10	91473809	0	64M	*	0	0	*	*	NM:i:5	MD:Z:3T10C32A3T4G7	AS:i:40",
        "length=64count=1731	256	chr5	150943085	0	48M1D16M	*	0	0	*	*	NM:i:5	MD:Z:14C10A22^A3T4G7	AS:i:38",
        "length=64count=1731	272	chr4	59796234	0	16H48M	*	0	0	*	*	NM:i:2	MD:Z:22T10G14	AS:i:38",
        "length=64count=1731	256	chr10	7672287	0	48M16H	*	0	0	*	*	NM:i:2	MD:Z:14C10A22	AS:i:38",
        "length=64count=1731	272	chr2	136058427	0	15M1I48M	*	0	0	*	*	NM:i:5	MD:Z:7C4A24T10G14	AS:i:38",
        "length=64count=1731	256	chr1	155894848	0	48M2D16M	*	0	0	*	*	NM:i:6	MD:Z:14C10A22^AA3T4G7	AS:i:38",
        "length=64count=1731	272	chr7	96208741	0	60M4H	*	0	0	*	*	NM:i:5	MD:Z:2C4C4A3T32G10	AS:i:38",
        "length=64count=1731	256	chr10	41488892	0	42M22H	*	0	0	*	*	NM:i:1	MD:Z:14C27	AS:i:37",
        "length=64count=1731	272	chr6	9678505	0	22H42M	*	0	0	*	*	NM:i:1	MD:Z:16T25	AS:i:37",
        "length=64count=1731	272	chr6	16023739	0	17H47M	*	0	0	*	*	NM:i:2	MD:Z:6A14T25	AS:i:37",
        "length=64count=1731	256	chr1	264859076	0	15H49M	*	0	0	*	*	NM:i:3	MD:Z:32A3T4G7	AS:i:34",
        "length=64count=1731	272	chr5	75185284	0	16M1D41M1D7M	*	0	0	*	*	NM:i:6	MD:Z:7C4A3^T22T10G7^C7	AS:i:33",
        "length=64count=1731	256	chr4	48349189	0	42M22H	*	0	0	*	*	NM:i:2	MD:Z:8T5C27	AS:i:32",
        "length=64count=1731	272	chr3	191721418	0	22H42M	*	0	0	*	*	NM:i:2	MD:Z:27G9A4	AS:i:32",
        "length=64count=1731	272	chr3	17704874	0	64M	*	0	0	*	*	NM:i:7	MD:Z:7C4A3T30A0A0G0G13	AS:i:32",
    };
    
    // ALignment without MD value (for testing SAMToGBSdbPlugin:calculateAlignedNumber when only have CIGAR
    String[] alignNoMD = {
            // five alternative alignments (with unique best primary alignment on forward strand):
            "length=64count=1731	0	chr4	34593346	11	64M	*	0	0	CAGCCGGACTAGCAAAACAAGTGCGTAAAATCGAAGTCTTCCTCAAAGGACCATCCATGCAAAA	ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff	AS:i:96	XS:i:88	XN:i:0	XM:i:4	XO:i:0	XG:i:0	NM:i:4	YT:Z:UU",
            "length=64count=1731	256	chr4	186122318	11	44M1D20M	*	0	0	CAGCCGGACTAGCAAAACAAGTGCGTAAAATCGAAGTCTTCCTCAAAGGACCATCCATGCAAAA	ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff	AS:i:88	XS:i:88	XN:i:0	XM:i:4	XO:i:1	XG:i:1	NM:i:5	YT:Z:UU",
            "length=64count=1731	272	chr6	16023721	11	16M1D48M	*	0	0	TTTTGCATGGATGGTCCTTTGAGGAAGACTTCGATTTTACGCACTTGTTTTGCTAGTCCGGCTG	ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff	AS:i:80	XS:i:88	XN:i:0	XM:i:5	XO:i:1	XG:i:1	NM:i:6	YT:Z:UU",
            "length=64count=1731	272	chr6	9678484	11	15M1I48M	*	0	0	TTTTGCATGGATGGTCCTTTGAGGAAGACTTCGATTTTACGCACTTGTTTTGCTAGTCCGGCTG	ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff	AS:i:78	XS:i:88	XN:i:0	XM:i:5	XO:i:1	XG:i:1	NM:i:6	YT:Z:UU",
            "length=64count=1731	272	chr8	115376743	11	25S7M1D32M	*	0	0	TTTTGCATGGATGGTCCTTTGAGGAAGACTTCGATTTTACGCACTTGTTTTGCTAGTCCGGCTG	ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff	AS:i:54	XS:i:88	XN:i:0	XM:i:2	XO:i:1	XG:i:1	NM:i:3	YT:Z:UU",
        };
    
    Position[] bwaMemAExpPositions = {   // Nb: some of these are hard clipped on the cut site end, so probably shouldn't be given positions
        new GeneralPosition.Builder(new Chromosome("4"),34593346).addAnno("forward","true").addAnno("cigar","64M").addAnno("supportvalue","44").addAnno("mappingapproach","bwaMem").strand((byte)1).build(),
        new GeneralPosition.Builder(new Chromosome("8"),15458003).addAnno("forward","true").addAnno("cigar","64M").addAnno("supportvalue","44").addAnno("mappingapproach","bwaMem").strand((byte)1).build(),
        new GeneralPosition.Builder(new Chromosome("1"),264859139).addAnno("forward","true").addAnno("cigar","64M").addAnno("supportvalue","44").addAnno("mappingapproach","bwaMem").strand((byte)1).build(),
        new GeneralPosition.Builder(new Chromosome("10"),91473731).addAnno("forward","true").addAnno("cigar","64M").addAnno("supportvalue","44").addAnno("mappingapproach","bwaMem").strand((byte)1).build(),
        new GeneralPosition.Builder(new Chromosome("7"),62589280).addAnno("forward","true").addAnno("cigar","64M").addAnno("supportvalue","44").addAnno("mappingapproach","bwaMem").strand((byte)1).build(),
        new GeneralPosition.Builder(new Chromosome("7"),62589358).addAnno("forward","true").addAnno("cigar","64M").addAnno("supportvalue","44").addAnno("mappingapproach","bwaMem").strand((byte)1).build(),
        new GeneralPosition.Builder(new Chromosome("7"),96208663+64-1).addAnno("forward","false").addAnno("cigar","64M").addAnno("supportvalue","44").addAnno("mappingapproach","bwaMem").strand((byte)1).build(),
        new GeneralPosition.Builder(new Chromosome("1"),252153050+47-1).addAnno("forward","false").addAnno("cigar","17H47M").addAnno("supportvalue","42").addAnno("mappingapproach","bwaMem").strand((byte)1).build(),
        new GeneralPosition.Builder(new Chromosome("4"),186122318).addAnno("forward","true").addAnno("cigar","47M17H").addAnno("supportvalue","42").addAnno("mappingapproach","bwaMem").strand((byte)1).build(),
        new GeneralPosition.Builder(new Chromosome("10"),91473809).addAnno("forward","true").addAnno("cigar","64M").addAnno("supportvalue","40").addAnno("mappingapproach","bwaMem").strand((byte)1).build(),
        new GeneralPosition.Builder(new Chromosome("5"),150943085).addAnno("forward","true").addAnno("cigar","48M1D16M").addAnno("supportvalue","38").addAnno("mappingapproach","bwaMem").strand((byte)1).build(),
        new GeneralPosition.Builder(new Chromosome("4"),59796234+48-1).addAnno("forward","false").addAnno("cigar","16H48M").addAnno("supportvalue","38").addAnno("mappingapproach","bwaMem").strand((byte)1).build(),
        new GeneralPosition.Builder(new Chromosome("10"),7672287).addAnno("forward","true").addAnno("cigar","48M16H").addAnno("supportvalue","38").addAnno("mappingapproach","bwaMem").strand((byte)1).build(),
        new GeneralPosition.Builder(new Chromosome("2"),136058427+15+48-1).addAnno("forward","false").addAnno("cigar","15M1I48M").addAnno("supportvalue","38").addAnno("mappingapproach","bwaMem").strand((byte)1).build(),
        new GeneralPosition.Builder(new Chromosome("1"),155894848).addAnno("forward","true").addAnno("cigar","48M2D16M").addAnno("supportvalue","38").addAnno("mappingapproach","bwaMem").strand((byte)1).build(),
        new GeneralPosition.Builder(new Chromosome("7"),96208741+60+4-1).addAnno("forward","false").addAnno("cigar","60M4H").addAnno("supportvalue","38").addAnno("mappingapproach","bwaMem").strand((byte)1).build(),
        new GeneralPosition.Builder(new Chromosome("10"),41488892).addAnno("forward","true").addAnno("cigar","42M22H").addAnno("supportvalue","37").addAnno("mappingapproach","bwaMem").strand((byte)1).build(),
        new GeneralPosition.Builder(new Chromosome("6"),9678505+42-1).addAnno("forward","false").addAnno("cigar","22H42M").addAnno("supportvalue","37").addAnno("mappingapproach","bwaMem").strand((byte)1).build(),
        new GeneralPosition.Builder(new Chromosome("6"),16023739+47-1).addAnno("forward","false").addAnno("cigar","17H47M").addAnno("supportvalue","37").addAnno("mappingapproach","bwaMem").strand((byte)1).build(),
        null, // the next one has too much hard clipping (15H) on the cut site end, so I am making it null
//        new GeneralPosition.Builder(new Chromosome("1"),264859076-15).addAnno("forward","true").addAnno("cigar","15H49M").addAnno("supportvalue","34").addAnno("mappingapproach","bwaMem").strand((byte)1).build(),
        new GeneralPosition.Builder(new Chromosome("5"),75185284+16+1+41+1+7-1).addAnno("forward","false").addAnno("cigar","16M1D41M1D7M").addAnno("supportvalue","33").addAnno("mappingapproach","bwaMem").strand((byte)1).build(),
        new GeneralPosition.Builder(new Chromosome("4"),48349189).addAnno("forward","true").addAnno("cigar","42M22H").addAnno("supportvalue","32").addAnno("mappingapproach","bwaMem").strand((byte)1).build(),
        new GeneralPosition.Builder(new Chromosome("3"),191721418+42-1).addAnno("forward","false").addAnno("cigar","22H42M").addAnno("supportvalue","32").addAnno("mappingapproach","bwaMem").strand((byte)1).build(),
        new GeneralPosition.Builder(new Chromosome("3"),17704874+64-1).addAnno("forward","false").addAnno("cigar","64M").addAnno("supportvalue","32").addAnno("mappingapproach","bwaMem").strand((byte)1).build(),
//        new GeneralPosition.Builder(new Chromosome(""),-1).addAnno("forward","").addAnno("cigar","").addAnno("supportvalue","").addAnno("mappingapproach","bwaMem").strand((byte)1).build(),
    };
    Tag[] bwaMemAExpTags = {
        TagBuilder.instance("CAGCCGGACTAGCAAAACAAGTGCGTAAAATCGAAGTCTTCCTCAAAGGACCATCCATGCAAAA").build(),
        TagBuilder.instance("CAGCCGGACTAGCAAAACAAGTGCGTAAAATCGAAGTCTTCCTCAAAGGACCATCCATGCAAAA").build(),
        TagBuilder.instance("CAGCCGGACTAGCAAAACAAGTGCGTAAAATCGAAGTCTTCCTCAAAGGACCATCCATGCAAAA").build(),
        TagBuilder.instance("CAGCCGGACTAGCAAAACAAGTGCGTAAAATCGAAGTCTTCCTCAAAGGACCATCCATGCAAAA").build(),
        TagBuilder.instance("CAGCCGGACTAGCAAAACAAGTGCGTAAAATCGAAGTCTTCCTCAAAGGACCATCCATGCAAAA").build(),
        TagBuilder.instance("CAGCCGGACTAGCAAAACAAGTGCGTAAAATCGAAGTCTTCCTCAAAGGACCATCCATGCAAAA").build(),
        TagBuilder.instance("CAGCCGGACTAGCAAAACAAGTGCGTAAAATCGAAGTCTTCCTCAAAGGACCATCCATGCAAAA").build(),
        TagBuilder.instance("CAGCCGGACTAGCAAAACAAGTGCGTAAAATCGAAGTCTTCCTCAAAGGACCATCCATGCAAAA").build(),
        TagBuilder.instance("CAGCCGGACTAGCAAAACAAGTGCGTAAAATCGAAGTCTTCCTCAAAGGACCATCCATGCAAAA").build(),
        TagBuilder.instance("CAGCCGGACTAGCAAAACAAGTGCGTAAAATCGAAGTCTTCCTCAAAGGACCATCCATGCAAAA").build(),
        TagBuilder.instance("CAGCCGGACTAGCAAAACAAGTGCGTAAAATCGAAGTCTTCCTCAAAGGACCATCCATGCAAAA").build(),
        TagBuilder.instance("CAGCCGGACTAGCAAAACAAGTGCGTAAAATCGAAGTCTTCCTCAAAGGACCATCCATGCAAAA").build(),
        TagBuilder.instance("CAGCCGGACTAGCAAAACAAGTGCGTAAAATCGAAGTCTTCCTCAAAGGACCATCCATGCAAAA").build(),
        TagBuilder.instance("CAGCCGGACTAGCAAAACAAGTGCGTAAAATCGAAGTCTTCCTCAAAGGACCATCCATGCAAAA").build(),
        TagBuilder.instance("CAGCCGGACTAGCAAAACAAGTGCGTAAAATCGAAGTCTTCCTCAAAGGACCATCCATGCAAAA").build(),
        TagBuilder.instance("CAGCCGGACTAGCAAAACAAGTGCGTAAAATCGAAGTCTTCCTCAAAGGACCATCCATGCAAAA").build(),
        TagBuilder.instance("CAGCCGGACTAGCAAAACAAGTGCGTAAAATCGAAGTCTTCCTCAAAGGACCATCCATGCAAAA").build(),
        TagBuilder.instance("CAGCCGGACTAGCAAAACAAGTGCGTAAAATCGAAGTCTTCCTCAAAGGACCATCCATGCAAAA").build(),
        TagBuilder.instance("CAGCCGGACTAGCAAAACAAGTGCGTAAAATCGAAGTCTTCCTCAAAGGACCATCCATGCAAAA").build(),
        TagBuilder.instance("CAGCCGGACTAGCAAAACAAGTGCGTAAAATCGAAGTCTTCCTCAAAGGACCATCCATGCAAAA").build(),
        TagBuilder.instance("CAGCCGGACTAGCAAAACAAGTGCGTAAAATCGAAGTCTTCCTCAAAGGACCATCCATGCAAAA").build(),
        TagBuilder.instance("CAGCCGGACTAGCAAAACAAGTGCGTAAAATCGAAGTCTTCCTCAAAGGACCATCCATGCAAAA").build(),
        TagBuilder.instance("CAGCCGGACTAGCAAAACAAGTGCGTAAAATCGAAGTCTTCCTCAAAGGACCATCCATGCAAAA").build(),
        TagBuilder.instance("CAGCCGGACTAGCAAAACAAGTGCGTAAAATCGAAGTCTTCCTCAAAGGACCATCCATGCAAAA").build(),
//        TagBuilder.instance("").build(),
    };
    
    // bwa cases
    String[] bwa = {
        // same tag as the first bowtie2 -k 5 case above (X0:i:10 = 10 best hits, but only 1 reported):
        "length=64count=1731	0	chr4	34593346	0	64M	*	0	0	CAGCCGGACTAGCAAAACAAGTGCGTAAAATCGAAGTCTTCCTCAAAGGACCATCCATGCAAAA	ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff	XT:A:R	NM:i:4	X0:i:10	X1:i:0	XM:i:4	XO:i:0	XG:i:0	MD:Z:25C21A3T4G7",
        // another tag with multiple positions reported (requires parseResult to return Tuple<Tag,Optional<Position>[]>):
        "length=64count=17	16	chr1	149645339	18	64M	*	0	0	CTGCACCCTTGGGCGACTTTTGCACGGAGAGTGGTTTTGTGCACTTATTTTGCTAGTCCGGCTG	ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff	XT:A:U	NM:i:1	X0:i:1	X1:i:3	XM:i:1	XO:i:0	XG:i:0	MD:Z:49G14	XA:Z:chr4,+98993172,64M,2;chr4,-33326770,64M,2;chr2,-51515187,64M,2;",
//        "",
    };
    Tag[] bwaExpTags = {
        TagBuilder.instance("CAGCCGGACTAGCAAAACAAGTGCGTAAAATCGAAGTCTTCCTCAAAGGACCATCCATGCAAAA").build(),
        TagBuilder.instance("CAGCCGGACTAGCAAAATAAGTGCACAAAACCACTCTCCGTGCAAAAGTCGCCCAAGGGTGCAG").build(),
//        TagBuilder.instance("").build(),
    };
    Position[] bwaExpPositions1 = {
        new GeneralPosition.Builder(new Chromosome("4"),34593346).addAnno("forward","true").addAnno("cigar","64M").addAnno("supportvalue","-4").addAnno("mappingapproach","bwa").strand((byte)1).build(),
//        new GeneralPosition.Builder(new Chromosome(""),0).addAnno("forward","").addAnno("cigar","").addAnno("supportvalue","").addAnno("mappingapproach","bwa").strand((byte)1).build(),
    };
    Position[] bwaExpPositions2 = {
        new GeneralPosition.Builder(new Chromosome("1"),149645339+64-1).addAnno("forward","false").addAnno("cigar","64M").addAnno("supportvalue","-1").addAnno("mappingapproach","bwa").strand((byte)1).build(),
        new GeneralPosition.Builder(new Chromosome("4"),98993172).addAnno("forward","true").addAnno("cigar","64M").addAnno("supportvalue","-2").addAnno("mappingapproach","bwa").strand((byte)1).build(),
        new GeneralPosition.Builder(new Chromosome("4"),33326770+64-1).addAnno("forward","false").addAnno("cigar","64M").addAnno("supportvalue","-2").addAnno("mappingapproach","bwa").strand((byte)1).build(),
        new GeneralPosition.Builder(new Chromosome("2"),51515187+64-1).addAnno("forward","false").addAnno("cigar","64M").addAnno("supportvalue","-2").addAnno("mappingapproach","bwa").strand((byte)1).build(),
//        new GeneralPosition.Builder(new Chromosome(""),0).addAnno("forward","").addAnno("cigar","").addAnno("supportvalue","").addAnno("mappingapproach","bwa").strand((byte)1).build(),
    };
    Position[][] bwaExpPositions = {bwaExpPositions1, bwaExpPositions2};

    @Test
    public void testAdjustForClipping() throws Exception {
        System.out.println("Begin test testAdjustForClipping ... ");
        String cigar = "3S8M1D6M4S";
        int cutPosition = 54;
 
        int[] pos =  SAMUtils.adjustCoordinates(cigar,cutPosition);
        assertEquals(pos[0],51); // answer is 51 54 - 3
        
        cigar = "48M";
        pos = SAMUtils.adjustCoordinates(cigar,cutPosition);
        assertEquals(pos[0],54); // nothing changes - no clipping
        
        cigar = "37H52M1D6M";
        pos = SAMUtils.adjustCoordinates(cigar,cutPosition);
        assertEquals(pos[0],17); // should be 54-37 = 17
        
        cigar = "52M1D6M26H";
        pos = SAMUtils.adjustCoordinates(cigar,cutPosition); // strand was true
        assertEquals(pos[0],54); // no clipping at beginning 
        assertEquals(pos[1],138); // add up to the end: 54 (start) + 52 + 6 + 26
        
        System.out.println("SAMToGBSdbPlugin: end testAdjustForClipping");
    }

    @Test
    public void testParseRowWithNegativePosition() throws Exception {
        // use reflection to access the private method SAMToGBSdbPlugin.parseRow(String inputLine)
        SAMToGBSdbPlugin sAMToGBSdb = new SAMToGBSdbPlugin();
        Class SAMToGBSdbClass = sAMToGBSdb.getClass();
        Method theMethod = SAMToGBSdbClass.getDeclaredMethod("parseRow", new Class[] {String.class});
        theMethod.setAccessible(true);
        
        // the sequence below has a "POS" of 1, but has 4 bp's soft-clipped from the beginning.
        // This should translate to a -3 as the starting position, which should result in 
        // parseRow tossing the read.
        String samLine = 
                "tagSeq=CAGCACGTGATGCATAACAACGCAATGCCTCCTCATACTCATGCACAACAGAATCACGGTTTTT        16      790     1       0      4S60M    *       0       0        AAAAACCGTGATTCTGTTGTGCATGAGTATGAGGAGGCATTGCGTTGTTATGCATCACGTGCTG      ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff AS:i:60 XS:i:60 XF:i:1  XE:i:1  NM:i:1";
        String[] tokens=samLine.split("\\s+");
        String origSeq = tokens[0].split("=")[1]; //s[0] is tagSeg=<original sequence here
        System.out.println("LCJ - in testParseRowWithNegavitePosition: number of tokens in string: " + tokens.length);
        
        Tag tag = TagBuilder.instance(origSeq).build();
        Tuple<Tag, Optional<Position>> expectedReturn = new Tuple<>(tag,Optional.<Position>empty());
        Tuple<Tag,Optional<Position>> actualReturn = (Tuple<Tag, Optional<Position>>) theMethod.invoke(sAMToGBSdb,samLine);
        assertEquals(expectedReturn,actualReturn);
    }
    
    // Working on this one
    @Ignore
    @Test
    public void testParseRow() throws Exception {
        // use reflection to access the private method SAMToGBSdbPlugin.parseRow(String inputLine)
        SAMToGBSdbPlugin SAMToGBSdb = new SAMToGBSdbPlugin();
        Class SAMToGBSdbClass = SAMToGBSdb.getClass();
        Method parseRow = SAMToGBSdbClass.getDeclaredMethod("parseRow", new Class[] {String.class});
        parseRow.setAccessible(true);
        
        SAMToGBSdb.mappingApproach("Bowtie2") ; // SAMToGBSdbPlugin.processData() sets this if "bowtie2" shows up in the file
        
        // loop through various bowtie2 default test cases
        for (int i = 0; i < bowtie2.length; i++) {
            System.out.println("bowtie2 line: "+i);
            Tuple<Tag,Optional<Position>> parsedResult = (Tuple<Tag,Optional<Position>>) parseRow.invoke(SAMToGBSdb,(Object) bowtie2[i]);
            Assert.assertEquals("Tags at testCase"+i+" differ:", bowtie2ExpTags[i], parsedResult.x);
        //    Assert.assertEquals("Positions at testCase"+i+" differ:", Optional.of(bowtie2ExpPositions[i]), parsedResult.y);
        //    Map.Entry<String, String>[] expAnnos = bowtie2ExpPositions[i].getAnnotation().getAllAnnotationEntries();
        //    for (Entry<String, String> expEntry : expAnnos) {
        //        Assert.assertTrue(
        //            "Position annotation "+expEntry.getKey()+"="+expEntry.getValue()+ " not present in testCase"+i,
        //            parsedResult.y.get().getAnnotation().isAnnotatedWithValue(expEntry.getKey(),expEntry.getValue())
        //        );
        //    }
        }
        
        // // bowtie2 -k 5 cases
        // for (int i = 0; i < bowtie2k5.length; i++) {
        //     Tuple<Tag,Optional<Position>> parsedResult = (Tuple<Tag,Optional<Position>>) parseRow.invoke(SAMToGBSdb,(Object) bowtie2k5[i]);
        //     Assert.assertEquals("Tags at testCase"+i+" differ:", bowtie2k5ExpTags[i], parsedResult.x);
        //     Assert.assertEquals("Positions at testCase"+i+" differ:", Optional.of(bowtie2k5ExpPositions[i]), parsedResult.y);
        //     Map.Entry<String, String>[] expAnnos = bowtie2k5ExpPositions[i].getAnnotation().getAllAnnotationEntries();
        //     for (Entry<String, String> expEntry : expAnnos) {
        //         Assert.assertTrue(
        //             "Position annotation "+expEntry.getKey()+"="+expEntry.getValue()+ " not present in testCase"+i,
        //             parsedResult.y.get().getAnnotation().isAnnotatedWithValue(expEntry.getKey(),expEntry.getValue())
        //         );
        //     }
        // }
        // 
        // SAMToGBSdb.isBowtie = false;
        // // bwa mem -a cases
        // for (int i = 0; i < bwaMemA.length; i++) {
        //     Tuple<Tag,Optional<Position>> parsedResult = (Tuple<Tag,Optional<Position>>) parseRow.invoke(SAMToGBSdb,(Object) bwaMemA[i]);
        //     Assert.assertEquals("Tags at testCase"+i+" differ:", bwaMemAExpTags[i], parsedResult.x);
        //     Assert.assertEquals("Positions at testCase"+i+" differ:", Optional.of(bwaMemAExpPositions[i]), parsedResult.y);
        //     Map.Entry<String, String>[] expAnnos = bwaMemAExpPositions[i].getAnnotation().getAllAnnotationEntries();
        //     for (Entry<String, String> expEntry : expAnnos) {
        //         Assert.assertTrue(
        //             "Position annotation "+expEntry.getKey()+"="+expEntry.getValue()+ " not present in testCase"+i,
        //             parsedResult.y.get().getAnnotation().isAnnotatedWithValue(expEntry.getKey(),expEntry.getValue())
        //         );
        //     }
        // }
    }
    
    //Skipping this test for now. There isn't a parseBWARow method any more, but is there some functionality here that should be incorporated into other tests?
    @Ignore
    @Test
    public void testParseBWARow() throws Exception {
        // use reflection to access the private method SAMToGBSdbPlugin.parseBWARow(String inputLine)
        SAMToGBSdbPlugin SAMToGBSdb = new SAMToGBSdbPlugin();
        Class SAMToGBSdbClass = SAMToGBSdb.getClass();
        Method parseBWARow = SAMToGBSdbClass.getDeclaredMethod("parseBWARow", new Class[] {String.class});
        parseBWARow.setAccessible(true);
        
        // bwa test cases (returns an array of Optional<Position>'s for each SAM row)
        for (int t = 0; t < bwa.length; t++) {
            Tuple<Tag,Optional<Position>[]> parsedResult = (Tuple<Tag,Optional<Position>[]>) parseBWARow.invoke(SAMToGBSdb,(Object) bwa[t]);
            Assert.assertEquals("Tags at testCase"+t+" differ:", bwaExpTags[t], parsedResult.x);
            for (int p = 0; p < parsedResult.y.length; p++) {
                Assert.assertEquals("Positions at testCase"+t+", position"+p+" differ:", Optional.of(bwaExpPositions[t][p]), parsedResult.y[p]);
                Map.Entry<String, String>[] expAnnos = bwaExpPositions[t][p].getAnnotation().getAllAnnotationEntries();
                for (Entry<String, String> expEntry : expAnnos) {
                    Assert.assertTrue(
                        "Position annotation "+expEntry.getKey()+"="+expEntry.getValue()+ " not present in testCase"+t,
                        parsedResult.y[p].get().getAnnotation().isAnnotatedWithValue(expEntry.getKey(),expEntry.getValue())
                    );
                }
            }
        }
    }

    @Test
    public void testHasMinAlignLength() throws Exception {
        // use reflection to access the private method SAMToGBSdbPlugin.hasMinAlignLength(String inputLine)
       //SAMToGBSdbPlugin SAMToGBSdb = new SAMToGBSdbPlugin();
        
        SAMToGBSdbPlugin SAMToGBSdb = new SAMToGBSdbPlugin()
        							.minAlignLength(52)
        							.minAlignProportion(0.0);
        						
        Class SAMToGBSdbClass = SAMToGBSdb.getClass();
        Method hasMinAlignLength = SAMToGBSdbClass.getDeclaredMethod("hasMinAlignLength", new Class[] {String[].class});
        hasMinAlignLength.setAccessible(true);
        

        String[] alignString=bowtie2[0].split("\\s");
        boolean hasMinAlign = (boolean)hasMinAlignLength.invoke(SAMToGBSdb,(Object)alignString);
        Assert.assertTrue("String has min length of 52",hasMinAlign); // 64 of 64 matched
        
        alignString = bowtie2[1].split("\\s");
        hasMinAlign = (boolean)hasMinAlignLength.invoke(SAMToGBSdb,(Object)alignString);
        Assert.assertFalse("String has min length of 52",hasMinAlign); // 51 of 64 matched - less than 52
        
        alignString = bowtie2[2].split("\\s");
        hasMinAlign = (boolean)hasMinAlignLength.invoke(SAMToGBSdb,(Object)alignString); // more to parse
        Assert.assertTrue("String has min length of 52",hasMinAlign); // 59 of 64 matched
        
        alignString = bowtie2[4].split("\\s");
        hasMinAlign = (boolean)hasMinAlignLength.invoke(SAMToGBSdb,(Object)alignString); // no CIGAR or MD
        Assert.assertFalse("String has min length of 52",hasMinAlign);
    }
    
    @Test
    public void testHasMinAlignProportion() throws Exception {
        // use reflection to access the private method SAMToGBSdbPlugin.hasMinAlignProportion(String[] inputLine)
       //SAMToGBSdbPlugin SAMToGBSdb = new SAMToGBSdbPlugin();
        
        SAMToGBSdbPlugin SAMToGBSdb = new SAMToGBSdbPlugin()
        							.minAlignLength(0)
        							.minAlignProportion(0.9);
        						
        Class SAMToGBSdbClass = SAMToGBSdb.getClass();
        Method hasMinAlignProportion = SAMToGBSdbClass.getDeclaredMethod("hasMinAlignProportion", new Class[] {String[].class});
        hasMinAlignProportion.setAccessible(true);
        
        // get alignment proportion test.  Returns proportion of string that has minimum sequence

        String[] alignString=bowtie2[0].split("\\s");
        boolean hasMinProportion = (boolean)hasMinAlignProportion.invoke(SAMToGBSdb,(Object)alignString);
        Assert.assertTrue("String matches 64 of 64",hasMinProportion); // 64 of 64 matched
        
        alignString=bowtie2[1].split("\\s");
        hasMinProportion = (boolean)hasMinAlignProportion.invoke(SAMToGBSdb,(Object)alignString);
        Assert.assertTrue("String matches 51 of 51",hasMinProportion); // 51 of 51 matched - string is shorter
        
        alignString=bowtie2[2].split("\\s");
        hasMinProportion = (boolean)hasMinAlignProportion.invoke(SAMToGBSdb,(Object)alignString); // more to parse
        Assert.assertTrue("String matches 59 of 64",hasMinProportion); // 59 of 64 matched
        
        alignString=bowtie2[4].split("\\s");
        hasMinProportion = (boolean)hasMinAlignProportion.invoke(SAMToGBSdb,(Object)alignString); // no MD, CIGAR=* (no matches)
        Assert.assertFalse("String has no matches",hasMinProportion);
    }
    
    @Test
    public void testCalculateNumberAligned() throws Exception {
        // use reflection to access the private method SAMToGBSdbPlugin.calculateNumberAligned(String[] samRead)
        SAMToGBSdbPlugin SAMToGBSdb = new SAMToGBSdbPlugin()
        							.minAlignLength(52)
        							.minAlignProportion(0.0);
        						
        Class SAMToGBSdbClass = SAMToGBSdb.getClass();
        Method calculateNumberAligned = SAMToGBSdbClass.getDeclaredMethod("calculateNumberAligned", new Class[] {String[].class});
        calculateNumberAligned.setAccessible(true);
        
        // Test using MD values
        String[] alignString=bowtie2[0].split("\\s");
        // From:http://yourmitra.wordpress.com/2008/09/26/using-java-reflection-to-invoke-a-method-with-array-parameters/
       //int numAligned = (int)calculateNumberAligned.invoke(SAMToGBSdb,new Object[]{alignString});
        int numAligned = (int)calculateNumberAligned.invoke(SAMToGBSdb, (Object)alignString);
        Assert.assertEquals(64,numAligned); // 64 of 64 matched: MD 64
        
        alignString=bowtie2[1].split("\\s");
        numAligned = (int)calculateNumberAligned.invoke(SAMToGBSdb,(Object)alignString);
        Assert.assertEquals(51,numAligned); // 51 of 51 matched: MD 51
        
        alignString=bowtie2[2].split("\\s");
        numAligned = (int)calculateNumberAligned.invoke(SAMToGBSdb,(Object)alignString); 
        Assert.assertEquals(59,numAligned); // 59 of 64 matched: MD 12T47
        
        alignString=bowtie2[4].split("\\s");
        numAligned = (int)calculateNumberAligned.invoke(SAMToGBSdb,(Object)alignString); 
        Assert.assertEquals(0,numAligned); // no MD, CIGAR=*
        
        alignString=bowtie2k5[2].split("\\s");
        numAligned = (int)calculateNumberAligned.invoke(SAMToGBSdb,(Object)alignString); 
        Assert.assertEquals(59,numAligned); // 59 of 64 matched: MD 7C4A3^T0T6A14T25
        
        // Using CIGAR values to find matches
        alignString=alignNoMD[0].split("\\s");       
        numAligned = (int)calculateNumberAligned.invoke(SAMToGBSdb,new Object[]{alignString});
        Assert.assertEquals(64,numAligned); // 64 of 64 matched: CIGAR 64M
        
        alignString=alignNoMD[2].split("\\s");
        numAligned = (int)calculateNumberAligned.invoke(SAMToGBSdb,(Object)alignString);
        Assert.assertEquals(64,numAligned); // 64 of 64 matched:  CIGAR 16M1D48M
        
        alignString=alignNoMD[3].split("\\s");
        numAligned = (int)calculateNumberAligned.invoke(SAMToGBSdb,(Object)alignString); 
        Assert.assertEquals(63,numAligned); // 59 of 64 matched: CIGAR 15M1I48M
        
        alignString=alignNoMD[4].split("\\s");
        numAligned = (int)calculateNumberAligned.invoke(SAMToGBSdb,(Object)alignString); // no CIGAR or MD
        Assert.assertEquals(39,numAligned);  // 39 of 64 matched: CIGAR 25S7M1D32M
    }

//   @Test
//    public void testRealignTags() throws Exception {
//        Tag aTag;  //assign
//        Position aPosition;
//        Genome aGenome;
//        Tuple<Integer,String> result=SAMToGBSdbPlugin.realignToGenome(aTag, aPostion, aGenome);
//
//
//    }
}

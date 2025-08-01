package net.maizegenetics.dna.tag;

import junit.framework.Assert;
import org.junit.Test;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import static org.junit.Assert.*;

public class TagBuilderTest {

    @Test
    public void testVariableLengthTags() throws Exception {
        Set<Tag> tagMap=new HashSet<>();

        String tagS="GACGGCAGTACTGGCTAGCTGAACTGTGGCCGTAGCTTTCAGCAGGAATTTTTTTTTTTTGCTG";
        Tag tag=TagBuilder.instance(tagS).build();
        Assert.assertEquals(tagS.length(), tag.seqLength());
        Assert.assertEquals(tagS, tag.sequence());
        tagMap.add(tag);
        Assert.assertTrue(tagMap.contains(tag));

        //63bp in length
        tagS="GACGGCAGTACTGGCTAGCTGAACTGTGGCCGTAGCTTTCAGCAGGAATTTTTTTTTTTTGCT";
        tag=TagBuilder.instance(tagS).reference().name("Bob").build();
        Assert.assertEquals(tagS.length(),tag.seqLength());
        Assert.assertEquals(tagS,tag.sequence());
        tagMap.add(tag);
        Assert.assertTrue(tagMap.contains(tag));
        Assert.assertTrue(tag.isReference());
        Assert.assertEquals("Bob", tag.name());

        tagS="CTGCCCCCCCTGGAATTCTCCATGGCGGCTG";
        tag=TagBuilder.instance(tagS).build();
        Assert.assertEquals(tagS.length(),tag.seqLength());
        Assert.assertEquals(tagS,tag.sequence());
        Assert.assertEquals("",tag.name());
        tagMap.add(tag);
        Assert.assertTrue(tagMap.contains(tag));
        Assert.assertFalse(tag.isReference());
    }

    @Test
    public void testReverseComplement() throws Exception {
        String tagS=      "GACGGCT";
        String tagReverse="AGCCGTC";
        Tag tag=TagBuilder.instance(tagS).build();
        Tag reverseTag=TagBuilder.reverseComplement(tag).build();
        Assert.assertEquals(tagReverse,reverseTag.sequence());
        Assert.assertFalse(reverseTag.equals(tag));
        Assert.assertEquals(tag.seqLength(),reverseTag.seqLength());
    }
}
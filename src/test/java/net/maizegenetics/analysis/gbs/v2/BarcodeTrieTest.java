package net.maizegenetics.analysis.gbs.v2;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableMap;
import junit.framework.Assert;
import net.maizegenetics.analysis.gbs.Barcode;
import org.junit.Test;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.OptionalInt;

import static org.junit.Assert.*;

public class BarcodeTrieTest {

    @Test
    public void testBarcodeParsing() throws Exception {
        GBSEnzyme enzyme=new GBSEnzyme("ApeKI");
        List<Barcode> barcodes=new ArrayList<>();
        barcodes.add(new Barcode("ACGT",enzyme.initialCutSiteRemnant(),"T1",1,"",""));
        barcodes.add(new Barcode("AGGG",enzyme.initialCutSiteRemnant(),"T2",2,"",""));
        barcodes.add(new Barcode("AGGGGT",enzyme.initialCutSiteRemnant(),"T3",3,"",""));
        barcodes.add(new Barcode("AGGGAAA",enzyme.initialCutSiteRemnant(),"T4",4,"",""));
        BarcodeTrie bt=new BarcodeTrie();
        bt.addAllBarcodes(barcodes);

        Map<String, Integer> tests= ImmutableMap.of(
                "ACGTCAGCTTTTTTTTTTTTTT", 1,
                "AGGGAAACAGCTTTTTTTTTTTTTT", 4,
                "AGGGAAACTGCTTTTTTTTTTTTTT", 4,
                "AGGGAAACTGCACGTACAGT", 4,
                "AGgGCAGCTTTTTTTTTTTTTT", 2
                );

        tests.entrySet().stream()
                .forEach(e -> {
                    Barcode b=bt.longestPrefix(e.getKey());
                    Assert.assertTrue(e.getValue() == b.getTaxaIndex());
                    System.out.println(e.getKey() +"="+b);
                });

        List<String> testList= ImmutableList.of(
                "ACNTCAGCTTTTTTTTTTTTTT",
                "A.GGAAACAGCTTTTTTTTTTTTTT",
                "A*GGAAACTGCTTTTTTTTTTTTTT",
                "A_GGAAACTGCACGTACAGT",
                "CGGGCAGCTTTTTTTTTTTTTT"
        );

        testList.stream()
                .forEach(s -> {
                    Barcode b=bt.longestPrefix(s);
                    Assert.assertNull(b);
                    System.out.println(s +"="+b);
                });

//        Barcode b1=bt.longestPrefix("ACGTCAGCTTTTTTTTTTTTTT");
//        System.out.println(b1);
//        Assert.assertEquals(1,b1.getTaxaIndex());
//        b=bt.longestPrefix("AAGTCAGCTTTTTTTTTTTTTT");
//        System.out.println(b);
//        Assert.assertNull(b);
//
//        b=bt.longestPrefix("AGGGCAGCTTTTTTTTTTTTTT");
//        System.out.println(b2);
//        Assert.assertNull(b2);
    }
}
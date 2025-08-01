/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.analysis.gbs;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import net.maizegenetics.dna.map.AbstractTagsOnPhysicalMap;
import net.maizegenetics.dna.map.TagMappingInfoV3;
import net.maizegenetics.dna.map.TagsOnPhysMapHDF5;
import net.maizegenetics.dna.map.TagsOnPhysicalMap;
import net.maizegenetics.dna.map.TagsOnPhysicalMapV3;
import net.maizegenetics.dna.tag.TagsByTaxa;
import net.maizegenetics.dna.tag.TagsByTaxaByteHDF5TagGroups;

/**
 * Stores physical position of tags. The physical position of tag is used to block the corresponding marker in genetic mapping if the tag is mapping to the marker coming from itself
 * Positions come from TOPM alignment hypothesis or the best position from machine learning prediction
 * @author Fei Lu
 */
public class TagBlockPosition {
    int[] blockChr;
    int[] blockPos;
    
    public TagBlockPosition (String tbtH5FileS, String topmFileS, int TOPMVersionValue) {
        long lastTimePoint = System.nanoTime();
        AbstractTagsOnPhysicalMap topm = null;
        if (TOPMVersionValue == 0) {
            topm = new TagsOnPhysicalMap (topmFileS, true);
        }
        else if (TOPMVersionValue == 1){
            topm = new TagsOnPhysMapHDF5(topmFileS);
        }
        else if (TOPMVersionValue == 2) {
            topm = new TagsOnPhysicalMapV3(topmFileS);
        }
        else {
            System.out.println("Input TOPM version value is not supported");
            System.exit(0);
        }
        TagsByTaxaByteHDF5TagGroups tbt = new TagsByTaxaByteHDF5TagGroups (tbtH5FileS);
        blockChr = new int[tbt.getTagCount()];
        blockPos = new int[tbt.getTagCount()];
        long[] t;
        int index;
        for (int i = 0; i < tbt.getTagCount(); i++) {
            t = tbt.getTag(i);
            index = topm.getTagIndex(t);
            if (index < 0) {
                blockChr[i] = Integer.MIN_VALUE;
                blockPos[i] = Integer.MIN_VALUE;
            }
            else {
                blockChr[i] = topm.getChromosome(index);
                blockPos[i] = topm.getStartPosition(index);
            }
        }
        System.out.println("Generating TagBlockPosition from TOPM took " + String.valueOf(this.getTimeSpanSecond(lastTimePoint)) + " seconds\n");
    }
    
    public TagBlockPosition (String tbtH5FileS, String topmH5FileS, String software) {
        long lastTimePoint = System.nanoTime();
        TagsByTaxaByteHDF5TagGroups tbt = new TagsByTaxaByteHDF5TagGroups (tbtH5FileS);
        System.out.println("Loading TBT HDF5 took " + String.valueOf(this.getTimeSpanSecond(lastTimePoint)) + " seconds");
        System.out.println("TBT has " + tbt.getTagCount() + " tags and " + tbt.getTaxaCount() + " taxa\n");
        lastTimePoint = System.nanoTime();
        TagMappingInfoV3.Aligner alignerName = TagMappingInfoV3.Aligner.getAlignerFromName(software);
        if (alignerName == null) {
            System.out.println("Input software is not Bowtie2, BWA or Blast, not supporting other aligner for now.");
            System.out.println("Program stops.");
            System.exit(0);
        }
        blockChr = new int[tbt.getTagCount()];
        blockPos = new int[tbt.getTagCount()];
        TagsOnPhysicalMapV3 topm = new TagsOnPhysicalMapV3(topmH5FileS);
        long[] t;
        int index;
        for (int i = 0; i < tbt.getTagCount(); i++) {
            t = tbt.getTag(i);
            index = topm.getTagIndex(t);
            if (index < 0) {
                blockChr[i] = Integer.MIN_VALUE;
                blockPos[i] = Integer.MIN_VALUE;
            }
            else {
                int[] chrPos = topm.getUniqueMappingOfAligner(index, alignerName);
                if (chrPos == null) {
                    blockChr[i] = Integer.MIN_VALUE;
                    blockPos[i] = Integer.MIN_VALUE;
                }
                else {
                    blockChr[i] = chrPos[0];
                    blockPos[i] = chrPos[1];
                }
            }
        }
        System.out.println("Generating TagBlockPosition from TOPM took " + String.valueOf(this.getTimeSpanSecond(lastTimePoint)) + " seconds\n");
    }
    
    public TagBlockPosition (String blockFileS) {
        this.readTagBlockPostition(blockFileS);
    }
    
    public void readTagBlockPostition (String blockFileS) {
        try {
            DataInputStream dis = new DataInputStream(new BufferedInputStream(new FileInputStream(blockFileS), 65536));
            int tagCount = dis.readInt();
            blockChr = new int[tagCount];
            blockPos = new int[tagCount];
            for (int i = 0; i < tagCount; i++) {
                blockChr[i] = dis.readInt();
                blockPos[i] = dis.readInt();
            }
            dis.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void writeTagBlockPosition (String blockFileS, int[] selectIndex) {
        try {
            DataOutputStream dos = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(blockFileS), 65536));
            dos.writeInt(selectIndex.length);
            for (int i = 0; i < selectIndex.length; i++) {
                dos.writeInt(blockChr[selectIndex[i]]);
                dos.writeInt(blockPos[selectIndex[i]]);
            }
            dos.flush();
            dos.close();
        }
        catch (Exception e) {
            System.out.println(e.toString());
            System.exit(1);
        }
    }
    
    public void writeTagBlockPosition (String blockFileS) {
        try {
            DataOutputStream dos = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(blockFileS), 65536));
            dos.writeInt(blockChr.length);
            for (int i = 0; i < blockChr.length; i++) {
                dos.writeInt(blockChr[i]);
                dos.writeInt(blockPos[i]);
            }
            dos.flush();
            dos.close();
        }
        catch (Exception e) {
            System.out.println(e.toString());
            System.exit(1);
        }
    }
    
    public int[] getBlockChr () {
        return blockChr;
    }
    
    public int[] getBlockPos () {
        return blockPos;
    }
    
    public int getBlockChr (int index) {
        return blockChr[index];
    }
    
    public int getBlockPos (int index) {
        return blockPos[index];
    }
    
    private double getTimeSpanSecond (long lastTimePoint) {
        return (double)this.getTimeSpanNano(lastTimePoint)/1000000000;
    }
    
    private long getTimeSpanNano (long lastTimePoint) {
        return System.nanoTime()- lastTimePoint;
    }
}

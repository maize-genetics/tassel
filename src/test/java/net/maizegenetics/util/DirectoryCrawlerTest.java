package net.maizegenetics.util;

import junit.framework.Assert;
import net.maizegenetics.constants.TutorialConstants;
import org.junit.Ignore;
import org.junit.Test;

import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;

public class DirectoryCrawlerTest {
    @Ignore
    @Test
    public void testGlobRegex() throws Exception {

        Path inDir = Paths.get(TutorialConstants.TUTORIAL_DIR);
        //test glob
        List<Path> txtList = DirectoryCrawler.listPaths("glob:*.txt", inDir);
        Assert.assertEquals("Number of text files not correct", 18, txtList.size());
        txtList = DirectoryCrawler.listPaths("glob:*geno*", inDir);
        Assert.assertEquals("Number of text files not correct", 20, txtList.size());
        //test regex
        txtList = DirectoryCrawler.listPaths("regex:.*\\.txt", inDir);
        Assert.assertEquals("Number of text files not correct", 18, txtList.size());
        //test default regex
        txtList = DirectoryCrawler.listPaths(".*\\.txt", inDir);
        Assert.assertEquals("Number of text files not correct", 18, txtList.size());

    }
}

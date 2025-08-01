package net.maizegenetics.analysis.data;

import net.maizegenetics.constants.GeneralConstants;
import org.junit.Test;

import java.util.Random;

/**
 * Defines xxxx
 *
 * @author Ed Buckler
 */
public class MigrateHDF5FromT4T5Test {
    static final String t4HDF5File=GeneralConstants.DATA_DIR+"CandidateTests/regHmpInT4.hmp.h5";

    @Test
    public void testMigrateT4WithGenotypes() throws Exception {
        Random r=new Random();
        MigrateHDF5FromT4T5.copyGenotypes(t4HDF5File,GeneralConstants.TEMP_DIR+"r"+r.nextInt(1000)+".t5.h5");

    }
}

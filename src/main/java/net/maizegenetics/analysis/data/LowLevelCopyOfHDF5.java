package net.maizegenetics.analysis.data;

import ch.systemsx.cisd.hdf5.HDF5Factory;
import ch.systemsx.cisd.hdf5.HDF5LinkInformation;
import ch.systemsx.cisd.hdf5.IHDF5Reader;
import ch.systemsx.cisd.hdf5.IHDF5Writer;
import net.maizegenetics.dna.snp.GenotypeTableBuilder;
import net.maizegenetics.dna.snp.NucleotideAlignmentConstants;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.util.HDF5Utils;
import net.maizegenetics.util.Tassel5HDF5Constants;

import java.util.List;

/**
 * Provides low level copy and migration tool.
 *
 * @author Ed Buckler
 */
public class LowLevelCopyOfHDF5 {
    public static void subsetGenotypesToNewFile(String t5File, String subT5File, TaxaList subsetTaxa) {
        IHDF5Reader reader=HDF5Factory.openForReading(t5File);
        IHDF5Writer writer=HDF5Factory.open(subT5File);

        writer.object().createGroup(Tassel5HDF5Constants.GENOTYPES_MODULE);
        HDF5Utils.unlockHDF5GenotypeModule(writer);
        writer.object().createGroup(Tassel5HDF5Constants.TAXA_MODULE);
        HDF5Utils.unlockHDF5TaxaModule(writer);

        int numTaxa = 0;
        HDF5Utils.writeHDF5GenotypesAlleleStates(writer,NucleotideAlignmentConstants.NUCLEOTIDE_ALLELES);
        HDF5Utils.writeHDF5GenotypesMaxNumAlleles(writer,NucleotideAlignmentConstants.NUMBER_NUCLEOTIDE_ALLELES);

        HDF5Utils.writeHDF5GenotypesRetainRareAlleles(writer,false);
        List<HDF5LinkInformation> fields = reader.object().getAllGroupMemberInformation(Tassel5HDF5Constants.GENOTYPES_MODULE, true);
        for (HDF5LinkInformation is : fields) {
            if (is.isGroup() == false) continue;
            String taxonName=is.getName();
            if(subsetTaxa.indexOf(taxonName)<0) continue;
           // System.out.println(taxonName);
            reader.object().copy(Tassel5HDF5Constants.GENOTYPES_MODULE+"/"+taxonName, writer,
                    Tassel5HDF5Constants.GENOTYPES_MODULE+"/"+taxonName);
            HDF5Utils.addTaxon(writer, new Taxon(taxonName));
            numTaxa++;
        }
        HDF5Utils.writeHDF5GenotypesNumTaxa(writer,numTaxa);
        HDF5Utils.writeHDF5TaxaNumTaxa(writer,numTaxa);

        //Position module

        reader.object().copy(Tassel5HDF5Constants.POSITION_MODULE, writer, Tassel5HDF5Constants.POSITION_MODULE);

        //Precalculated Stats
        GenotypeTableBuilder.annotateHDF5File(writer);

        
        HDF5Utils.lockHDF5GenotypeModule(writer);
        HDF5Utils.lockHDF5TaxaModule(writer);
        reader.close();
        writer.close();
    }

}

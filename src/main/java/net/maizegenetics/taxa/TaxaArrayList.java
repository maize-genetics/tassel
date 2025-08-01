package net.maizegenetics.taxa;

import it.unimi.dsi.fastutil.objects.Object2IntOpenHashMap;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.util.List;
import java.util.Map;

/**
 * In memory immutable instance of {@link TaxaList}. Basic list of taxa
 * (samples) that are used in Alignments and other purposes.
 *
 * Use {@link TaxaListBuilder} to instantiate.
 *
 * @author Ed Buckler
 * @author Terry Casstevens
 */
class TaxaArrayList extends TaxaNoIndexList {

    private static final Logger myLogger = LogManager.getLogger(TaxaArrayList.class);
    private final Map<String, Integer> myNameToIndex;

    TaxaArrayList(List<Taxon> srcList) {
        super(srcList);
        myNameToIndex = new Object2IntOpenHashMap(numberOfTaxa());
        for (int i = 0; i < numberOfTaxa(); i++) {
            myNameToIndex.put(srcList.get(i).getName(), i);
        }
    }

    @Override
    public int indexOf(String name) {
        Integer index = myNameToIndex.get(name);
        if (index == null) return -1;
        return index;
    }

    @Override
    public boolean contains(Object o) {
        if (o instanceof String) {
            Integer index = myNameToIndex.get((String) o);
            if (index == null) return false;
            else return true;
        } else if (o instanceof Taxon) {
            Integer index = myNameToIndex.get(((Taxon) o).getName());
            if (index == null) return false;
            else return true;
        } else {
            return false;
        }
    }

}

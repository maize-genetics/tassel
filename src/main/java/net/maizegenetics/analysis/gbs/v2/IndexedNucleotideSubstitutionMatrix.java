package net.maizegenetics.analysis.gbs.v2;

import java.util.IdentityHashMap;
import java.util.List;
import java.util.Map;

import org.biojava.nbio.core.alignment.template.SubstitutionMatrix;
import org.biojava.nbio.core.sequence.compound.NucleotideCompound;
import org.biojava.nbio.core.sequence.template.CompoundSet;

/**
 * A drop-in wrapper around a BioJava nucleotide {@link SubstitutionMatrix} that makes
 * {@link #getValue(NucleotideCompound, NucleotideCompound)} O(1).
 * <p>
 * BioJava's {@code SimpleSubstitutionMatrix.getValue(from,to)} maps each compound to its matrix
 * row/column with a linear {@code List.indexOf()} on <em>every</em> dynamic-programming cell.
 * Profiling the GBSv2 alignment phase showed
 * {@code SimpleSubstitutionMatrix.getValue -> getIndexOfCompound -> ArrayList.indexOf} was ~28% of
 * the whole pipeline. This wrapper precomputes a compound&rarr;index map and a dense score table,
 * populated once by calling the delegate's {@code getValue}, so the scores are byte-identical to the
 * delegate; each subsequent lookup is two {@link HashMap} gets and an array read.
 * <p>
 * {@link #getCompoundSet()} returns the caller-supplied compound set (the aligned sequences' set) so
 * it passes the reference-equality check in {@code Alignments.getMultipleSequenceAlignment}.
 * {@code AbstractProfileProfileAligner} scores using {@code query.getCompoundSet().getAllCompounds()}
 * and {@code matrix.getValue(...)} (not {@code matrix.getCompoundSet()}), so building the table over
 * that same compound set means every lookup the aligner makes hits the fast path.
 */
final class IndexedNucleotideSubstitutionMatrix implements SubstitutionMatrix<NucleotideCompound> {

    private final SubstitutionMatrix<NucleotideCompound> delegate;
    private final CompoundSet<NucleotideCompound> compoundSet;
    private final Map<NucleotideCompound, Integer> index;
    private final short[][] scores;

    IndexedNucleotideSubstitutionMatrix(SubstitutionMatrix<NucleotideCompound> delegate,
                                        CompoundSet<NucleotideCompound> compoundSet) {
        this.delegate = delegate;
        this.compoundSet = compoundSet;
        List<NucleotideCompound> all = compoundSet.getAllCompounds();
        int n = all.size();
        // IdentityHashMap, not HashMap: the aligner passes the compound-set's own singleton
        // NucleotideCompound instances, so identity lookup works and we avoid NucleotideCompound
        // .hashCode()/equals() (which are toString().hashCode()/equalsIgnoreCase — String garbage
        // per lookup). Any non-identity compound simply falls back to delegate.getValue below.
        this.index = new IdentityHashMap<>(n * 2);
        for (int i = 0; i < n; i++) {
            index.put(all.get(i), i);
        }
        this.scores = new short[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                scores[i][j] = delegate.getValue(all.get(i), all.get(j)); // slow lookup, but only n*n times, once
            }
        }
    }

    @Override
    public short getValue(NucleotideCompound from, NucleotideCompound to) {
        Integer i = index.get(from);
        Integer j = index.get(to);
        if (i == null || j == null) {
            return delegate.getValue(from, to); // compound outside the prebuilt set: exact (slow) fallback
        }
        return scores[i][j];
    }

    @Override
    public CompoundSet<NucleotideCompound> getCompoundSet() {
        return compoundSet;
    }

    // Everything else delegates unchanged (not on the per-cell scoring path).
    @Override public short[][] getMatrix() { return delegate.getMatrix(); }
    @Override public String getMatrixAsString() { return delegate.getMatrixAsString(); }
    @Override public short getMaxValue() { return delegate.getMaxValue(); }
    @Override public short getMinValue() { return delegate.getMinValue(); }
    @Override public String getName() { return delegate.getName(); }
    @Override public String getDescription() { return delegate.getDescription(); }
    @Override public void setName(String name) { delegate.setName(name); }
    @Override public void setDescription(String description) { delegate.setDescription(description); }
    @Override public Map<NucleotideCompound, Short> getRow(NucleotideCompound row) { return delegate.getRow(row); }
    @Override public Map<NucleotideCompound, Short> getColumn(NucleotideCompound column) { return delegate.getColumn(column); }
    @Override public SubstitutionMatrix<NucleotideCompound> normalizeMatrix(short scale) { return delegate.normalizeMatrix(scale); }
}

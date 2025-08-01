/*
 *  GenericArrayTaxa
 * 
 *  Created on Oct 10, 2016
 */
package net.maizegenetics.analysis.avro;

import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.ListIterator;
import net.maizegenetics.taxa.TaxaList;
import org.apache.avro.Schema;
import org.apache.avro.generic.GenericArray;
import org.apache.avro.generic.GenericRecord;

/**
 *
 * @author Terry Casstevens
 */
public class GenericArrayTaxa implements GenericArray<GenericRecord> {
    
    private final TaxaList myTaxa;
    private final int myNumTaxa;
    
    public GenericArrayTaxa(TaxaList taxa) {
        myTaxa = taxa;
        myNumTaxa = taxa.numberOfTaxa();
    }
    
    @Override
    public GenericRecord peek() {
        throw new UnsupportedOperationException("Not supported yet.");
    }
    
    @Override
    public void reverse() {
        throw new UnsupportedOperationException("Not supported yet.");
    }
    
    @Override
    public int size() {
        return myNumTaxa;
    }
    
    @Override
    public boolean isEmpty() {
        throw new UnsupportedOperationException("Not supported yet.");
    }
    
    @Override
    public boolean contains(Object o) {
        throw new UnsupportedOperationException("Not supported yet.");
    }
    
    @Override
    public Iterator<GenericRecord> iterator() {
        
        return new Iterator<GenericRecord>() {
            
            private int myIndex = 0;
            
            @Override
            public boolean hasNext() {
                return myIndex < myNumTaxa;
            }
            
            @Override
            public GenericRecord next() {
                return new GenericRecordTaxon(myTaxa.get(myIndex++));
            }
        };
        
    }
    
    @Override
    public Object[] toArray() {
        throw new UnsupportedOperationException("Not supported yet.");
    }
    
    @Override
    public <T> T[] toArray(T[] a) {
        throw new UnsupportedOperationException("Not supported yet.");
    }
    
    @Override
    public boolean add(GenericRecord e) {
        throw new UnsupportedOperationException("Not supported yet.");
    }
    
    @Override
    public boolean remove(Object o) {
        throw new UnsupportedOperationException("Not supported yet.");
    }
    
    @Override
    public boolean containsAll(Collection<?> c) {
        throw new UnsupportedOperationException("Not supported yet.");
    }
    
    @Override
    public boolean addAll(Collection<? extends GenericRecord> c) {
        throw new UnsupportedOperationException("Not supported yet.");
    }
    
    @Override
    public boolean addAll(int index, Collection<? extends GenericRecord> c) {
        throw new UnsupportedOperationException("Not supported yet.");
    }
    
    @Override
    public boolean removeAll(Collection<?> c) {
        throw new UnsupportedOperationException("Not supported yet.");
    }
    
    @Override
    public boolean retainAll(Collection<?> c) {
        throw new UnsupportedOperationException("Not supported yet.");
    }
    
    @Override
    public void clear() {
        throw new UnsupportedOperationException("Not supported yet.");
    }
    
    @Override
    public GenericRecord get(int index) {
        throw new UnsupportedOperationException("Not supported yet.");
    }
    
    @Override
    public GenericRecord set(int index, GenericRecord element) {
        throw new UnsupportedOperationException("Not supported yet.");
    }
    
    @Override
    public void add(int index, GenericRecord element) {
        throw new UnsupportedOperationException("Not supported yet.");
    }
    
    @Override
    public GenericRecord remove(int index) {
        throw new UnsupportedOperationException("Not supported yet.");
    }
    
    @Override
    public int indexOf(Object o) {
        throw new UnsupportedOperationException("Not supported yet.");
    }
    
    @Override
    public int lastIndexOf(Object o) {
        throw new UnsupportedOperationException("Not supported yet.");
    }
    
    @Override
    public ListIterator<GenericRecord> listIterator() {
        throw new UnsupportedOperationException("Not supported yet.");
    }
    
    @Override
    public ListIterator<GenericRecord> listIterator(int index) {
        throw new UnsupportedOperationException("Not supported yet.");
    }
    
    @Override
    public List<GenericRecord> subList(int fromIndex, int toIndex) {
        throw new UnsupportedOperationException("Not supported yet.");
    }
    
    @Override
    public Schema getSchema() {
        return AvroConstants.TAXA_SCHEMA;
    }
    
}

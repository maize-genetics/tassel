/*
 *  GenericMapAnnotations
 * 
 *  Created on Oct 11, 2016
 */
package net.maizegenetics.analysis.avro;

import java.util.Arrays;
import java.util.Collection;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import net.maizegenetics.util.GeneralAnnotation;

/**
 *
 * @author Terry Casstevens
 */
public class GenericMapAnnotations implements Map<String, String> {

    private final GeneralAnnotation myAnnotations;

    public GenericMapAnnotations(GeneralAnnotation annotations) {
        myAnnotations = annotations;
    }

    @Override
    public int size() {
        if (myAnnotations == null) {
            return 0;
        } else {
            return myAnnotations.numAnnotations();
        }
    }

    @Override
    public boolean isEmpty() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public boolean containsKey(Object key) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public boolean containsValue(Object value) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public String get(Object key) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public String put(String key, String value) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public String remove(Object key) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public void putAll(Map<? extends String, ? extends String> m) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public void clear() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public Set<String> keySet() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public Collection<String> values() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public Set<Entry<String, String>> entrySet() {
        if (myAnnotations == null) {
            return null;
        } else {
            return new HashSet<>(Arrays.asList(myAnnotations.getAllAnnotationEntries()));
        }
    }

}

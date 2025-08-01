/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.util;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;
import java.util.Stack;

/**
 * Functions used for Graph traversal and analysis
 * @author Eli Rodgers-Melnick
 */
public final class GraphUtils<T> {
    /**
     * Produces edges in a depth-first search starting at source and labeled
     * -1,0,1 for forward, reverse, and nontree (direction type). Based on
     * http://www.ics.uci.edu/~eppstein/PADS/DFS.py by D. Eppstein, July, 2004
     * @param G A graph
     * @param source The source node of the graph
     * @return A tuple of an arraylist containing the edges and an ArrayList containing the
     * direction types
     */
    public static <T> Tuple<ArrayList<Tuple<T,T>>, ArrayList<Byte>> dfsLabeledEdges(Graph<T> G, T source) {
        Set<T> visited = new HashSet();
        ArrayList<Tuple<T,T>> edges = new ArrayList();
        ArrayList<Byte> directions = new ArrayList();
        // Add source reference
        edges.add(new Tuple(source,source));
        directions.add((byte)1);
        visited.add(source);
        // Instantiate stack
        Stack<Tuple<T,Iterator<T>>> stack = new Stack();
        stack.push(new Tuple(source, G.neighbors(source).iterator()));
        while (!stack.empty()) {
            Tuple<T, Iterator<T>> toVisit = stack.peek();
            // Check for child
            if(toVisit.y.hasNext()) {
                T child = toVisit.y.next();
                // Check if child already visited
                if (visited.contains(child)) {
                    // This is a nontree type direction. Do not add child to stack
                    edges.add(new Tuple(toVisit.x, child));
                    directions.add((byte)0);
                } else {
                    // This is a forward type direction. Add child to stack after 
                    // adding edge to return edges
                    edges.add(new Tuple(toVisit.x, child));
                    directions.add((byte)1);
                    stack.push(new Tuple(child, G.neighbors(child).iterator()));
                    // Add child to visited
                    visited.add(child);
                }
            } else {
                // Remove from top of stack
                stack.pop();
                // If stack not empty, put in a direction reversal
                if (!stack.empty()) {
                    Tuple<T, Iterator<T>> nextVisit = stack.peek();
                    edges.add(new Tuple(nextVisit.x, toVisit.x));
                    directions.add((byte)-1);
                }
            }
        }
        // Put in reversal for source node
        edges.add(new Tuple(source, source));
        directions.add((byte)-1);
        return new Tuple(edges, directions);
    }
    /**
     * Produces nodes in a depth-first search pre-ordering starting at source (i.e. listing node
     * starting with the source node)
     * @param <T> The class of the node
     * @param G A graph
     * @param source The source node of the graph
     * @return An ArrayList with a depth-first search pre-ordering starting from source
     */
    public static <T> ArrayList<T> dfsPreorderNodes(Graph<T> G, T source) {
        Tuple<ArrayList<Tuple<T,T>>, ArrayList<Byte>> dfs = dfsLabeledEdges(G, source);
        ArrayList<T> pre = new ArrayList();
        // Go through labeled edges, adding all forward traversals target nodes
        for (int i = 0; i < dfs.x.size(); i++) {
            if (dfs.y.get(i) == 1) {
                pre.add(dfs.x.get(i).y);
            }
        }
        return pre;
    }
    /**
     * Produces nodes in a depth-first search post-ordering starting at source (i.e. listing node
     * starting with the last node and headed toward source)
     * @param <T> The class of the node
     * @param G A graph
     * @param source The source node of the graph
     * @return An ArrayList with a depth-first search post-ordering starting from source
     */
    public static <T> ArrayList<T> dfsPostorderNodes(Graph<T> G, T source) {
        Tuple<ArrayList<Tuple<T,T>>, ArrayList<Byte>> dfs = dfsLabeledEdges(G, source);
        ArrayList<T> pre = new ArrayList();
        // Go through labeled edges, adding all reverse traversals target nodes
        for (int i = 0; i < dfs.x.size(); i++) {
            if (dfs.y.get(i) == -1) {
                pre.add(dfs.x.get(i).y);
            }
        }
        return pre;
    }
}

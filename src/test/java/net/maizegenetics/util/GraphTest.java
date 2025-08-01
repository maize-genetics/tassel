/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.util;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author Eli Rodgers-Melnick
 */
public class GraphTest {
    UndirectedGraph<Integer> gU;
    DirectedGraph<Integer> gD;
    public GraphTest() {
        /** Undirected graph **/
        GraphBuilder<Integer> undirectedGraphBuilder = new GraphBuilder(Graph.GraphType.UNDIRECTED);
        /** Directed graph **/
        GraphBuilder<Integer> directedGraphBuilder = new GraphBuilder(Graph.GraphType.DIRECTED);
        // Create graph of 3 nodes (1,2,3), as a path
        undirectedGraphBuilder.addEdge(1, 2).addEdge(2, 3);
        directedGraphBuilder.addEdge(1, 2).addEdge(2, 3);
        gU = (UndirectedGraph)undirectedGraphBuilder.build();
        gD = (DirectedGraph)directedGraphBuilder.build();
    }
    
    @BeforeClass
    public static void setUpClass() {
    }
    
    @AfterClass
    public static void tearDownClass() {
    }
    
    @Before
    public void setUp() {
    }
    
    @After
    public void tearDown() {
    }
    @Test
    public void testNodesIterUndirected() {
        System.out.println("Testing nodesIter() on undirected graph");
        HashSet<Integer> expected = new HashSet();
        HashSet<Integer> observed = new HashSet();
        for (int i=1; i<= 3; i++) {
            expected.add(i);
        }
        Iterator<Integer> it = gU.nodesIter();
        while (it.hasNext()) {
            observed.add(it.next());
        }
        assertEquals("Set mismatch: ",expected, observed);
    }
    @Test
    public void testNodesIterDirected() {
        System.out.println("Testing nodesIter() on directed graph");
        HashSet<Integer> expected = new HashSet();
        HashSet<Integer> observed = new HashSet();
        for (int i=1; i<= 3; i++) {
            expected.add(i);
        }
        Iterator<Integer> it = gD.nodesIter();
        while (it.hasNext()) {
            observed.add(it.next());
        }
        assertEquals("Set mismatch: ",expected, observed);
    }
    @Test
    public void testNodesUndirected() {
        System.out.println("Testing nodes() on undirected graph");
        HashSet<Integer> expected = new HashSet();
        for (int i=1; i<= 3; i++) {
            expected.add(i);
        }
        assertEquals("Set mismatch: ",expected, gU.nodes());
    }
    @Test
    public void testNodesDirected() {
        System.out.println("Testing nodes() on directed graph");
        HashSet<Integer> expected = new HashSet();
        for (int i=1; i<= 3; i++) {
            expected.add(i);
        }
        assertEquals("Set mismatch: ",expected, gD.nodes());
    }
    @Test
    public void testHasNodeUndirected() {
        System.out.println("Testing hasNode() on undirected graph");
        assertTrue("hasNode mismatch", gU.hasNode(1));
        assertFalse("hasNode mismatch", gU.hasNode(4));
    }
    @Test
    public void testHasNodeDirected() {
        System.out.println("Testing hasNode() on directed graph");
        assertTrue("hasNode mismatch", gD.hasNode(1));
        assertFalse("hasNode mismatch", gD.hasNode(4));
    }
    @Test
    public void testNumberOfNodesUndirected() {
        System.out.println("Testing numberOfNodes() on undirected graph");
        assertEquals("Size mismatch: ", gU.numberOfNodes(), 3);
    }
    @Test
    public void testNumberOfNodesDirected() {
        System.out.println("Testing numberOfNodes() on directed graph");
        assertEquals("Size mismatch: ", gD.numberOfNodes(), 3);
    }
    @Test
    public void testHasEdgeUndirected() {
        System.out.println("Testing hasEdge() on undirected graph");
        assertTrue("Edge mismatch: ", gU.hasEdge(1, 2));
        assertTrue("Edge mismatch: ", gU.hasEdge(2, 1));
        assertFalse("Edge mismatch: ", gU.hasEdge(1, 3));
        assertFalse("Edge mismatch: ", gU.hasEdge(3, 1));
    }
    @Test
    public void testHasEdgeDirected() {
        System.out.println("Testing hasEdge() on directed graph");
        assertTrue("Edge mismatch: ", gD.hasEdge(1, 2));
        assertTrue("Edge mismatch: ", gD.hasEdge(2, 3));
        assertFalse("Edge mismatch: ", gD.hasEdge(3, 2));
        assertFalse("Edge mismatch: ", gD.hasEdge(2, 1));
        assertFalse("Edge mismatch: ", gD.hasEdge(1, 3));
        assertFalse("Edge mismatch: ", gD.hasEdge(3, 1));
    }
    @Test
    public void testNeighborsUndirected() {
        System.out.println("Testing neighbors() on undirected graph");
        Collection<Integer> neighbors1 = gU.neighbors(1);
        Collection<Integer> neighbors2 = gU.neighbors(2);
        assertEquals("Neighborhood size error: ",neighbors1.size(), 1);
        assertTrue("Neighborhood error: ", neighbors1.contains(2));
        assertEquals("Neighborhood size error: ", neighbors2.size(),2);
        assertTrue("Neighborhood error: ", neighbors2.contains(1));
        assertTrue("Neighborhood error: ", neighbors2.contains(3));
    }
    @Test
    public void testNeighborsDirected() {
        System.out.println("Testing neighbors() on directed graph");
        Collection<Integer> neighbors1 = gD.neighbors(1);
        Collection<Integer> neighbors2 = gD.neighbors(2);
        assertEquals("Neighborhood size error: ",neighbors1.size(), 1);
        assertTrue("Neighborhood error: ", neighbors1.contains(2));
        assertEquals("Neighborhood size error: ", neighbors2.size(),1);
        assertFalse("Neighborhood error: ", neighbors2.contains(1));
        assertTrue("Neighborhood error: ", neighbors2.contains(3));
    }
    @Test
    public void testEdgesUndirected() {
        System.out.println("Testing edges() on undirected graph");
        Collection<Map.Entry<Integer, Integer>> edges = gU.edges();
        assertEquals("Size mismatch: ", 2, edges.size());
    }
    @Test
    public void testEdgesDirected() {
        System.out.println("Testing edges() on directed graph");
        Collection<Map.Entry<Integer, Integer>> edges = gD.edges();
        assertEquals("Size mismatch: ", 2, edges.size());
    }
    @Test
    public void testEdgesIterUndirected() {
        System.out.println("Testing edgesIter() on undirected graph");
        Iterator<Map.Entry<Integer, Integer>> it = gU.edgesIter();
        HashSet<Tuple<Integer, Integer>> added = new HashSet();
        while (it.hasNext()) {
            Map.Entry<Integer, Integer> entry = it.next();
            assertFalse("Itererror: ",added.contains(new Tuple(entry.getValue(), 
                    entry.getKey())));
            added.add(new Tuple(entry.getKey(),entry.getValue()));
        }
        assertEquals("Size mismatch: ", 2, added.size());
    }
    @Test
    public void testEdgesIterDirected() {
        System.out.println("Testing edgesIter() on directed graph");
        Iterator<Map.Entry<Integer, Integer>> it = gD.edgesIter();
        HashSet<Tuple<Integer, Integer>> added = new HashSet();
        while (it.hasNext()) {
            Map.Entry<Integer, Integer> entry = it.next();
            assertFalse("Itererror: ",added.contains(new Tuple(entry.getValue(), 
                    entry.getKey())));
            added.add(new Tuple(entry.getKey(),entry.getValue()));
        }
        assertEquals("Size mismatch: ", 2, added.size());
    }
    @Test
    public void testInDegree() {
        System.out.println("Testing inDegree for directed graph");
        assertEquals("In degree mismatch: ", 1, gD.inDegree(2));
        assertEquals("In degree mismatch: ", 0, gD.inDegree(1));
    }
    @Test
    public void testOutDegree() {
        System.out.println("Testing outDegree for directed graph");
        assertEquals("Out degree mismatch: ", 1, gD.outDegree(2));
        assertEquals("Out degree mismatch: ", 1, gD.outDegree(1));
        assertEquals("Out degree mismatch: ", 0, gD.outDegree(3));
    }
    @Test
    public void testSuccessors() {
        System.out.println("Testing successors for directed graph");
        Collection<Integer> successors2 = gD.successors(2);
        assertEquals("Successors size mismatch: ", 1, successors2.size());
        assertTrue("Successors mismatch: ", successors2.contains(3));
    }
    @Test
    public void testPredecessors() {
        System.out.println("Testing predecessors for directed graph");
        Collection<Integer> predecessors2 = gD.predecessors(2);
        assertEquals("Predecessors size mismatch: ", 1, predecessors2.size());
        assertTrue("Predecessors mismatch: ", predecessors2.contains(1));
    }
    @Test
    public void testSizeDirected() {
        System.out.println("Testing size of undirected graph");
        assertEquals("Graph size mismatch: ", 2, gD.size());
        assertEquals("Graph size mismatch: ", (Double)2., (Double)gD.size(true));
    }
    @Test
    public void testSizeUndirected() {
        System.out.println("Testing size of directed graph");
        assertEquals("Graph size mismatch: ", 2, gU.size());
        assertEquals("Graph size mismatch: ", (Double)2., (Double)gU.size(true));
    }
}
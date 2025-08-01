/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.util;
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
public class GraphBuilderTest {
    private GraphBuilder<Integer> undirectedGraphBuilder;
    private GraphBuilder<Integer> directedGraphBuilder;
    public GraphBuilderTest() {}
    @BeforeClass
    public static void setUpClass() {
    }

    @AfterClass
    public static void tearDownClass() {
    }
    @Before
    public void setUp() {
        directedGraphBuilder = new GraphBuilder(Graph.GraphType.DIRECTED);
        undirectedGraphBuilder = new GraphBuilder(Graph.GraphType.UNDIRECTED);
    }
    @After
    public void tearDown() {
    }
    @Test
    public void testAddNodeUndirected() {
        System.out.println("Testing addNode on Undirected Graph");
        undirectedGraphBuilder.addNode(1);
        UndirectedGraph g = (UndirectedGraph)undirectedGraphBuilder.build();
        assertEquals("Size Mismatch: ", 1, g.numberOfNodes());
    }
    @Test
    public void testAddNodeDirected() {
        System.out.println("Testing addNode on Directed Graph");
        directedGraphBuilder.addNode(1);
        DirectedGraph g = (DirectedGraph)directedGraphBuilder.build();
        assertEquals("Size Mismatch: ", 1, g.numberOfNodes());
    }
    @Test
    public void testAddEdgeUndirected() {
        System.out.println("Testing addEdge on Undirected Graph");
        undirectedGraphBuilder.addEdge(1, 2);
        UndirectedGraph g = (UndirectedGraph)undirectedGraphBuilder.build();
        assertEquals("Size Mismatch: ", 2, g.numberOfNodes());
        assertEquals("Size Mismatch: ", 1, g.size());
    }
    @Test
    public void testAddEdgeDirected() {
        System.out.println("Testing addEdge on Directed Graph");
        directedGraphBuilder.addEdge(1, 2);
        DirectedGraph g = (DirectedGraph)directedGraphBuilder.build();
        assertEquals("Size Mismatch: ", 2, g.numberOfNodes());
        assertEquals("Size Mismatch: ", 1, g.size());
    }
}

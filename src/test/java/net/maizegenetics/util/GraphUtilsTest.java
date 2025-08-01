/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.util;

import java.util.ArrayList;
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
public class GraphUtilsTest {
    UndirectedGraph<Integer> gU;
    DirectedGraph<Integer> gD;
    public GraphUtilsTest() {
        /** Undirected graph **/
        GraphBuilder<Integer> undirectedGraphBuilder = new GraphBuilder(Graph.GraphType.UNDIRECTED);
        /** Directed graph **/
        GraphBuilder<Integer> directedGraphBuilder = new GraphBuilder(Graph.GraphType.DIRECTED);
        // Create graph of 3 nodes (1,2,3), as a path
        undirectedGraphBuilder.addEdge(1, 2).addEdge(2, 3).addEdge(2, 4).addEdge(3, 6).addEdge(4, 5);
        directedGraphBuilder.addEdge(1, 2).addEdge(2, 3).addEdge(2, 4).addEdge(3, 6).addEdge(4, 5);
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
    public void testDfsPreorderNodesDirected() {
        System.out.println("Testing DFS preorder nodes on directed graph");
        int[] exp_order = {1,2,3,6,4,5};
        ArrayList<Integer> pre_dfs = GraphUtils.dfsPreorderNodes(gD, 1);
        for (int i = 0; i < pre_dfs.size(); i++) {
            assertEquals("Order error: ", exp_order[i], (int)pre_dfs.get(i));
        }
    }
    @Test
    public void testDfsPreorderNodesUnirected() {
        System.out.println("Testing DFS preorder nodes on undirected graph");
        int[] exp_order = {1,2,3,6,4,5};
        ArrayList<Integer> pre_dfs = GraphUtils.dfsPreorderNodes(gU, 1);
        for (int i = 0; i < pre_dfs.size(); i++) {
            assertEquals("Order error: ", exp_order[i], (int)pre_dfs.get(i));
        }
    }
    @Test
    public void testDfsPostorderNodesDirected() {
        System.out.println("Testing DFS postorder nodes on directed graph");
        int[] exp_order = {6,3,5,4,2,1};
        ArrayList<Integer> post_dfs = GraphUtils.dfsPostorderNodes(gD, 1);
        for (int i = 0; i < post_dfs.size(); i++) {
            assertEquals("Order error: ", exp_order[i], (int)post_dfs.get(i));
        }
    }
    @Test
    public void testDfsPostorderNodesUndirected() {
        System.out.println("Testing DFS postorder nodes on undirected graph");
        int[] exp_order = {6,3,5,4,2,1};
        ArrayList<Integer> post_dfs = GraphUtils.dfsPostorderNodes(gU, 1);
        for (int i = 0; i < post_dfs.size(); i++) {
            assertEquals("Order error: ", exp_order[i], (int)post_dfs.get(i));
        }
    }
}
#include <gtest/gtest.h>
#include "graph.hpp"
#include "exceptions.hpp"

TEST(Graph, InitGraph)
{
    Matrix<double> graphMat = Matrix<double>::random(4,4, 0.0, 10.0 );
    Graph g = Graph(graphMat);
    ASSERT_TRUE( graphMat.compare(g.getGraph()));

}

TEST(Graph, InitGraphFail)
{
    Matrix<double> graphMat = Matrix<double>::random(4,3, 0.0, 10.0 );
    ASSERT_THROW(Graph g = Graph(graphMat), SquareMatrixException);
}


/*

 v0 ----> v1 ----> v2 ----> v0

 */

TEST(Graph, DirectConnected)
{
    double gData[] = {0,1,0,  0,0,1,  1,0,0};
    Matrix<double> gMat = Matrix<double>(3,3, gData);

    Graph g = Graph(gMat);
    ASSERT_TRUE( g.isDirectedGraph() );

    ASSERT_TRUE(g.isAdjacent(0, 1) );
    ASSERT_TRUE(g.isAdjacent(1, 2) );
    ASSERT_TRUE(g.isAdjacent(2, 0) );

    ASSERT_FALSE(g.isAdjacent(1, 0) );
    ASSERT_FALSE(g.isAdjacent(2, 1) );
    ASSERT_FALSE(g.isAdjacent(0, 2) );

    g.toUndirectedGraph();
    ASSERT_FALSE( g.isDirectedGraph() );

    ASSERT_TRUE(g.isAdjacent(0, 1) );
    ASSERT_TRUE(g.isAdjacent(1, 2) );
    ASSERT_TRUE(g.isAdjacent(2, 0) );

    ASSERT_TRUE(g.isAdjacent(1, 0) );
    ASSERT_TRUE(g.isAdjacent(2, 1) );
    ASSERT_TRUE(g.isAdjacent(0, 2) );
}

/*

 v0 ----> v1 ----> v2 ----> v1

 */

TEST(Graph, Connected)
{
    double gData[] = {0,1,0,  0,0,1,  0,1,0};
    Matrix<double> gMat = Matrix<double>(3,3, gData);

    Graph g = Graph(gMat);

    ASSERT_TRUE(g.isConnected(0, 1) );
    ASSERT_TRUE(g.isConnected(0, 2) );
    ASSERT_TRUE(g.isConnected(1, 2) );
    ASSERT_TRUE(g.isConnected(2, 1) );

    ASSERT_FALSE(g.isConnected(2, 0) );
    ASSERT_FALSE(g.isConnected(1, 0) );
}
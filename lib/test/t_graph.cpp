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

    ASSERT_TRUE( g.isDirectlyConnected(0,1) );
    ASSERT_TRUE( g.isDirectlyConnected(1,2) );
    ASSERT_TRUE( g.isDirectlyConnected(2,0) );

    ASSERT_FALSE( g.isDirectlyConnected(1,0) );
    ASSERT_FALSE( g.isDirectlyConnected(2,1) );
    ASSERT_FALSE( g.isDirectlyConnected(0,2) );
}
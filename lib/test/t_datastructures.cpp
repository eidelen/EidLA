#include <gtest/gtest.h>
#include <inc/datastructures.hpp>

TEST(Heap, NodeConstructor)
{
    HeapNode<int,4> node(5);

    ASSERT_EQ( node.m_children.size(), 4 );
    ASSERT_EQ( node.m_value, 5 );

    for( HeapNode<int,4>* n: node.m_children )
        ASSERT_EQ(n, nullptr);



}
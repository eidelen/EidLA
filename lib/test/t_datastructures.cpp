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

TEST(Heap, NodeInsert)
{
    std::vector<int> data = {3,2,0,1};
    HeapNode<int,3> root(4);

    for( int i : data )
        ASSERT_TRUE(root.insert(i));

    ASSERT_FALSE(root.insert(5));

    ASSERT_EQ( root.m_value, 4 );
    ASSERT_EQ( root.m_children.at(0)->m_value, 3 );
    ASSERT_EQ( root.m_children.at(1)->m_value, 2 );
    ASSERT_EQ( root.m_children.at(2)->m_value, 0 );
    ASSERT_EQ( root.m_children.at(0)->m_children.at(0)->m_value, 1 );

    std::cout << &root << std::endl;
}

TEST(Heap, HeapInsert)
{
    std::vector<int> data = {4,3,2,0,1};
    Heap<int,3>* heap = new Heap<int,3>();

    for( int i : data )
        heap->insert(i);

    ASSERT_EQ( heap->m_root->m_value, 4 );
    ASSERT_EQ( heap->m_root->m_children.at(0)->m_value, 3 );
    ASSERT_EQ( heap->m_root->m_children.at(1)->m_value, 2 );
    ASSERT_EQ( heap->m_root->m_children.at(2)->m_value, 0 );
    ASSERT_EQ( heap->m_root->m_children.at(0)->m_children.at(0)->m_value, 1 );

    std::cout << heap;

    // crack top node

    heap->insert(5);

    std::cout << heap;

    ASSERT_EQ( heap->m_root->m_value, 5 );
}

TEST(Heap, Find)
{
    std::vector<int> data = {6,5,3,2,1};
    Heap<int,3>* heap = new Heap<int,3>();

    for( int i : data )
        heap->insert(i);

    std::cout << heap;

    for( int i : data )
        ASSERT_TRUE(heap->find(i) != nullptr);

    ASSERT_TRUE(heap->find(-1) == nullptr);


}
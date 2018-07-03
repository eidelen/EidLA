#include <gtest/gtest.h>
#include <inc/datastructures.hpp>
#include "matrix.hpp"

TEST(Heap, NodeConstructor)
{
    HeapNodeMax<int,4> node(5);

    ASSERT_EQ( node.m_children.size(), 4 );
    ASSERT_EQ( node.m_value, 5 );
    ASSERT_EQ( node.m_nbrOfElements, 1 );

    for( HeapNode<int,4>* n: node.m_children )
        ASSERT_EQ(n, nullptr);
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

    delete heap;

}

bool maxHeapCheck( HeapNode<int, 2>* node )
{
    int node_val = node->m_value;

    for( auto e : node->m_children )
    {
        if( e != nullptr )
        {
            if( maxHeapCheck(e) )
            {
                // children value must be smaller or equal. if bigger, invalid heap
                if( node_val < e->m_value )
                    return false;
            }
            else
            {
                return false;
            }
        }
    }

    return true;
}

TEST(Heap, BigHeap)
{
    size_t nbrOfElementsToInsert = 20000;
    Matrix<int> data = Matrix<int>::random(nbrOfElementsToInsert, 1, 0, 10000);
    Heap<int, 2>* binaryHeap = new Heap<int, 2>();

    for(size_t i = 0; i < data.getNbrOfElements(); i++)
        binaryHeap->insert(data.data()[i]);

    ASSERT_EQ( binaryHeap->m_root->m_nbrOfElements, nbrOfElementsToInsert );

    // assure max is in root
    auto max = data.max();
    ASSERT_EQ( std::get<2>(max), binaryHeap->m_root->m_value);

    // assure all values in heap
    for(size_t i = 0; i < data.getNbrOfElements(); i++)
    {
        int val = data.data()[i];
        ASSERT_TRUE( binaryHeap->find(val) != nullptr );
    }

    // value -1 not in heap
    ASSERT_TRUE( binaryHeap->find(-1) == nullptr );

    // check heap is valid
    ASSERT_TRUE(maxHeapCheck( binaryHeap->m_root ));

    delete binaryHeap;
}

TEST(Heap, MinNodeConstructor)
{
    HeapNodeMin<int,4>* node = new HeapNodeMin<int,4>(5);

    ASSERT_EQ( node->m_children.size(), 4 );
    ASSERT_EQ( node->m_value, 5 );

    for( HeapNode<int,4>* n: node->m_children )
        ASSERT_EQ(n, nullptr);

    delete node;
}

TEST(Heap, MinNodeInsert)
{
    std::vector<int> data = {3,2,0,1};
    HeapNodeMin<int,3> root(0);

    for( int i : data )
        ASSERT_TRUE(root.insert(i));

    std::cout << &root << std::endl;

    ASSERT_EQ( root.m_value, 0 );

    ASSERT_EQ( root.m_children.at(0)->m_value, 3);
    ASSERT_EQ( root.m_children.at(1)->m_value, 2 );
    ASSERT_EQ( root.m_children.at(2)->m_value, 0 );

    ASSERT_EQ( root.m_children.at(2)->m_children.at(0)->m_value, 1);
}

TEST(Heap, BigMinHeap)
{
    Matrix<int> data = Matrix<int>::random(10000, 1, 0, 10000);
    Heap<int, 2>* binaryHeap = new Heap<int, 2>(HeapType::Min);

    for(size_t i = 0; i < data.getNbrOfElements(); i++)
        binaryHeap->insert(data.data()[i]);

    // assure min is in root
    auto min = data.min();
    ASSERT_EQ( std::get<2>(min), binaryHeap->m_root->m_value);

    // assure all values in heap
    for(size_t i = 0; i < data.getNbrOfElements(); i++)
    {
        int val = data.data()[i];
        ASSERT_TRUE( binaryHeap->find(val) != nullptr );
    }

    // value -1 not in heap
    ASSERT_TRUE( binaryHeap->find(-1) == nullptr );

    delete binaryHeap;
}

TEST(Heap, NodeDepth)
{
    HeapNode<int,2>* root = new HeapNodeMax<int,2>(5);
    ASSERT_EQ( root->depth(), 1);

    ASSERT_TRUE(root->insert(0));
    ASSERT_EQ( root->depth(), 2);

    ASSERT_TRUE(root->insert(4));
    ASSERT_EQ( root->depth(), 2);
    ASSERT_TRUE(root->insert(3));
    ASSERT_TRUE(root->insert(2));
    ASSERT_TRUE(root->insert(1));

    ASSERT_EQ( root->depth(), 4 );
    ASSERT_EQ( root->m_children.at(0)->depth(), 1);
    ASSERT_EQ( root->m_children.at(1)->depth(), 3);

    delete root;


    HeapNode<int,1>* listNode = new HeapNodeMax<int,1>(10);
    listNode->insert(9);
    listNode->insert(8);
    listNode->insert(7);
    listNode->insert(6);

    ASSERT_EQ( listNode->depth(), 5);

    delete listNode;
}



TEST(Heap, Balance)
{
    Matrix<int> data = Matrix<int>::random(300, 1, 0, 1000);
    Heap<int, 3>* binaryHeap = new Heap<int, 3>(HeapType::Min);

    for(size_t i = 0; i < data.getNbrOfElements(); i++)
        binaryHeap->insert(data.data()[i]);

    for( const HeapNode<int,3>* n: binaryHeap->m_root->m_children )
        std::cout << "Child depth = " << n->depth() << ", nbr of elements = " << n->m_nbrOfElements << std::endl;

    delete binaryHeap;
}

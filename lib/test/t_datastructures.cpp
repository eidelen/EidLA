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

    ASSERT_EQ(binaryHeap->size(), 0);

    for(size_t i = 0; i < data.getNbrOfElements(); i++)
        binaryHeap->insert(data.data()[i]);

    ASSERT_EQ(binaryHeap->size(), nbrOfElementsToInsert);

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
    Matrix<int> data = Matrix<int>::random(16, 1, 0, 100);
    Heap<int, 2>* binaryHeap = new Heap<int, 2>(HeapType::Max);

    for(size_t i = 0; i < data.getNbrOfElements(); i++)
        binaryHeap->insert(data.data()[i]);

    for( const HeapNode<int,2>* n: binaryHeap->m_root->m_children )
        std::cout << "Child depth = " << n->depth() << ", nbr of elements = " << n->m_nbrOfElements << std::endl;

    delete binaryHeap;
}


TEST(BST, NodeConstruct)
{
    BSTNode<int>* node = new BSTNode<int>(5);

    ASSERT_EQ(node->m_val, 5);
    ASSERT_TRUE(node->m_right == nullptr);
    ASSERT_TRUE(node->m_left == nullptr);

    node->insert(5); // again 5 -> nothing happens
    ASSERT_TRUE(node->m_right == nullptr);
    ASSERT_TRUE(node->m_left == nullptr);

    node->insert(6);
    ASSERT_EQ(node->m_right->m_val,6);

    node->insert(4);
    ASSERT_EQ(node->m_right->m_val, 6);

    delete node;
}

BST<int>* getTestBST()
{
    BST<int>* bst = new BST<int>();
    int data[] = {8,3,10,1,6,14,13,4,9,12,7};
    for( const int& k : data )
    {
        bst->insert(k);
    }

    return bst;
}

TEST(BST, Insert)
{
    auto bst = getTestBST();

    ASSERT_EQ( bst->m_root->m_val, 8 );

    ASSERT_EQ( bst->m_root->m_left->m_val, 3 );
    ASSERT_EQ( bst->m_root->m_right->m_val, 10 );

    ASSERT_EQ( bst->m_root->m_left->m_left->m_val, 1);
    ASSERT_EQ( bst->m_root->m_left->m_left->m_left, nullptr);
    ASSERT_EQ( bst->m_root->m_left->m_left->m_right, nullptr);
    ASSERT_EQ( bst->m_root->m_left->m_right->m_val, 6);
    ASSERT_EQ( bst->m_root->m_left->m_right->m_left->m_val, 4);
    ASSERT_EQ( bst->m_root->m_left->m_right->m_right->m_val, 7);

    ASSERT_EQ( bst->m_root->m_right->m_right->m_val, 14 );
    ASSERT_EQ( bst->m_root->m_right->m_left->m_val, 9 );
    ASSERT_EQ( bst->m_root->m_right->m_right->m_left->m_val, 13 );
    ASSERT_EQ( bst->m_root->m_right->m_right->m_right, nullptr );
    ASSERT_EQ( bst->m_root->m_right->m_right->m_left->m_left->m_val, 12 );
    ASSERT_EQ( bst->m_root->m_right->m_right->m_left->m_right, nullptr );

    delete bst;
}

TEST(BST, Search)
{
    auto bst = getTestBST();
    int bstVals[] = {8,3,10,1,6,14,13,4,9,12,7};

    for( const int& k : bstVals )
    {
        ASSERT_TRUE( bst->search(k) );
    }

    ASSERT_FALSE( bst->search(20) );
    ASSERT_FALSE( bst->search(0) );

    delete bst;
}

TEST(BST, NodeType)
{
    auto bst = getTestBST();

    ASSERT_EQ( BSTNode<int>::getNodeType(nullptr), BSTNode<int>::NodeType::NoNode);
    ASSERT_EQ( BSTNode<int>::getNodeType( bst->m_root ), BSTNode<int>::NodeType::TwoChildren);
    ASSERT_EQ( BSTNode<int>::getNodeType( bst->m_root->m_left->m_left ), BSTNode<int>::NodeType::NoChild);
    ASSERT_EQ( BSTNode<int>::getNodeType( bst->m_root->m_right->m_right ), BSTNode<int>::NodeType::OnlyLeftChild);

    bst->insert(2);
    ASSERT_EQ( BSTNode<int>::getNodeType( bst->m_root->m_left->m_left ), BSTNode<int>::NodeType::OnlyRightChild);
}

TEST(BST, LeftMostChild)
{
    auto bst = getTestBST();
    ASSERT_EQ(bst->m_root->getLeftMostChild()->m_val, 1);
    ASSERT_EQ(bst->m_root->m_right->getLeftMostChild()->m_val, 9);
    ASSERT_EQ(bst->m_root->m_right->m_right->getLeftMostChild()->m_val, 12);
    ASSERT_EQ(bst->m_root->m_right->m_right->getLeftMostChild()->getLeftMostChild()->m_val, 12);
}

TEST(BST, FindRemoveProblem)
{
    // linked list
    BSTNode<int>* root = new BSTNode<int>(0);
    root->insert(1);
    root->insert(2);
    root->insert(3);

    ASSERT_EQ(root->m_val, 0);
    ASSERT_EQ(root->m_left, nullptr);

    ASSERT_EQ(root->m_right->m_val, 1);
    ASSERT_EQ(root->m_right->m_left, nullptr);

    ASSERT_EQ(root->m_right->m_right->m_val, 2);
    ASSERT_EQ(root->m_right->m_right->m_left, nullptr);

    ASSERT_EQ(root->m_right->m_right->m_right->m_val, 3);
    ASSERT_EQ(root->m_right->m_right->m_right->m_left, nullptr);
    ASSERT_EQ(root->m_right->m_right->m_right->m_right, nullptr);

    root->remove(1);

    delete root;
}

TEST(BST, Remove)
{
    auto bst = getTestBST();

    // remove when no child -> delete node
    ASSERT_EQ( bst->m_root->m_left->m_right->m_left->m_val, 4);
    bst->remove(4);
    ASSERT_EQ( bst->m_root->m_left->m_right->m_left, nullptr);

    // remove when 1 child -> 6 replaced by 7
    ASSERT_EQ( bst->m_root->m_left->m_right->m_val, 6);
    bst->remove(6);
    ASSERT_EQ( bst->m_root->m_left->m_right->m_val, 7);

    // remove when 1 child -> replace node with child node
    ASSERT_EQ( bst->m_root->m_right->m_right->m_val, 14 );
    bst->remove(14);
    ASSERT_EQ( bst->m_root->m_right->m_right->m_val, 13 ); // 13 was only child of 14, therefore it was replaced

    // remove when 2 children -> replace node left most element of right subtree
    // 10 is replaced with 12
    ASSERT_EQ( bst->m_root->m_right->m_val, 10 );
    bst->remove(10);
    ASSERT_EQ( bst->m_root->m_right->m_val, 12 );
    ASSERT_EQ( bst->m_root->m_right->m_left->m_val, 9 );
    ASSERT_EQ( bst->m_root->m_right->m_right->m_val, 13 );
    ASSERT_EQ( bst->m_root->m_right->m_right->m_left, nullptr );
    ASSERT_EQ( bst->m_root->m_right->m_right->m_right, nullptr );

    bst->remove(12);
    ASSERT_EQ( bst->m_root->m_right->m_left->m_val, 9 );
    ASSERT_EQ( bst->m_root->m_right->m_val, 13 );
    ASSERT_EQ( bst->m_root->m_right->m_right, nullptr );

    delete bst;
}
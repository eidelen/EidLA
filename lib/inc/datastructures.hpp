/****************************************************************************
** Copyright (c) 2017 Adrian Schneider
**
** Permission is hereby granted, free of charge, to any person obtaining a
** copy of this software and associated documentation files (the "Software"),
** to deal in the Software without restriction, including without limitation
** the rights to use, copy, modify, merge, publish, distribute, sublicense,
** and/or sell copies of the Software, and to permit persons to whom the
** Software is furnished to do so, subject to the following conditions:
**
** The above copyright notice and this permission notice shall be included in
** all copies or substantial portions of the Software.
**
** THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
** IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
** FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
** AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
** LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
** FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
** DEALINGS IN THE SOFTWARE.
**
*****************************************************************************/

#ifndef MY_DATASTRUCTURES_H
#define MY_DATASTRUCTURES_H

#include <vector>
#include <cmath>
#include <iostream>

template <typename T, size_t C>
class HeapNode
{
public:
    HeapNode( T value ) : m_children(C), m_value(value), m_nbrOfElements(1)
    {
        std::fill(m_children.begin(), m_children.end(), nullptr );
    }

    virtual ~HeapNode()
    {
        for( HeapNode<T,C>* e : m_children )
            delete e;
    }

    bool insert(T val)
    {
        // val must be smaller/bigger or equal to this node's value
        if( compareFunction(val, m_value) )
            return false;

        // if empty child place, put it there
        ssize_t newChildIdx = getEmptyChildIdx();
        if( newChildIdx >= 0 )
        {
            m_children[newChildIdx] = createNode(val);
            m_nbrOfElements++;
            return true;
        }

        // all children spots are occupied
        // try to insert into every child
        for( HeapNode<T,C>* e : m_children )
        {
            if( e->insert(val) )
            {
                reorderChildren();
                m_nbrOfElements++;
                return true;
            }
        }

        // the value is bigger/smaller than any of the children ones, but smaller/bigger or equal
        // of the this node. replace the smallest/biggest child with the new value node.
        size_t replace_idx = childToReplace();
        HeapNode<T,C>* newNode = createNode(val);
        newNode->m_nbrOfElements = 1 + m_children.at(replace_idx)->m_nbrOfElements; // overtake number of elements from former node
        newNode->m_children.at(0) = m_children.at(replace_idx);
        m_children.at(replace_idx) = newNode;
        reorderChildren();

        m_nbrOfElements++;
        return true;
    }

    ssize_t getEmptyChildIdx() const
    {
        ssize_t idx = -1;
        for( ssize_t i = 0; i < m_children.size(); i++ )
            if( m_children.at(i) == nullptr )
                return i;

        return -1;
    }

    HeapNode<T,C>* find(T val)
    {
        // if val is bigger/smaller than the value in this node,
        // the value val cannot be found further down this node.
        if( compareFunction(val, m_value) )
            return nullptr;

        if( m_value == val )
            return this;

        for( HeapNode<T,C>* e : m_children )
        {
            if (e != nullptr)
            {
                HeapNode<T, C> *found = e->find(val);
                if (found != nullptr)
                    return found;
            }
        }

        return nullptr;
    }

    size_t depth() const
    {
        // Find max depth among children
        size_t maxDepth = 0;
        for( HeapNode<T,C>* e : m_children )
        {
            if( e != nullptr )
            {
                size_t cDepth = e->depth();
                if( cDepth > maxDepth )
                    maxDepth = cDepth;
            }
        }

        return maxDepth + 1;
    }

    virtual bool compareFunction(T a, T b) const = 0;
    virtual size_t childToReplace() const = 0;
    virtual HeapNode<T,C>* createNode(T val) const = 0;

public:
    std::vector<HeapNode<T,C>*> m_children;
    size_t m_nbrOfElements;
    T m_value;

private:
    static bool sortChildren(const HeapNode<T,C>* a, const HeapNode<T,C>* b)
    {
        return a->m_nbrOfElements < b->m_nbrOfElements;
    }

    void reorderChildren()
    {
        // reorder children according to their number of containing elements, so that
        // the next insertion is first tried at the child with the least amount of elements.
        // this keeps the heap balanced
        std::sort( m_children.begin(), m_children.end(), sortChildren );
    }

};

template <typename T, size_t C>
class HeapNodeMax: public HeapNode<T,C>
{
public:
    HeapNodeMax( T value ) : HeapNode<T,C>(value)
    {
    }

    bool compareFunction(T a, T b) const override
    {
        return a > b;
    }

    size_t childToReplace() const override
    {
        // find the index of the child to replace. in max heap,
        // it is the child with the smallest value.
        size_t smallest_idx = 0;
        T smallestVal = std::numeric_limits<T>::max();

        for( size_t i = 0; i < this->m_children.size(); i++ )
        {
            T cVal = this->m_children.at(i)->m_value;
            if( cVal < smallestVal )
            {
                smallestVal = cVal;
                smallest_idx = i;
            }
        }

        return smallest_idx;
    }

    virtual HeapNode<T,C>* createNode(T val) const override
    {
        return new HeapNodeMax<T,C>(val);
    }
};

template <typename T, size_t C>
class HeapNodeMin: public HeapNode<T,C>
{
public:
    HeapNodeMin( T value ) : HeapNode<T,C>(value)
    {
    }

    bool compareFunction(T a, T b) const override
    {
        return a < b;
    }

    size_t childToReplace() const override
    {
        // find the index of the child to replace. in min heap,
        // it is the child with the biggest value.
        size_t biggest_idx = 0;
        T biggest_val = std::numeric_limits<T>::min();

        for( size_t i = 0; i < this->m_children.size(); i++ )
        {
            T cVal = this->m_children.at(i)->m_value;
            if( cVal > biggest_val )
            {
                biggest_val = cVal;
                biggest_idx = i;
            }
        }

        return biggest_idx;
    }

    virtual HeapNode<T,C>* createNode(T val) const override
    {
        return new HeapNodeMin<T,C>(val);
    }
};

enum class HeapType
{
    Max,
    Min
};

template <typename T, size_t C>
class Heap
{

public:

    Heap(HeapType heapType = HeapType::Max) : m_root(nullptr), m_type(heapType)
    {
    }

    ~Heap()
    {
        if(m_root)
            delete m_root;
    }

    void insert(T val)
    {
        if( !m_root )
        {
            // first element -> root does not exist yet
            m_root = createNode(val);
        }
        else
        {
            // try to insert
            if(!(m_root->insert(val)))
            {
                // insert does not work ->
                // val will become the new root element.
                HeapNode<T,C>* tmp = m_root;
                m_root = createNode(val);
                m_root->m_nbrOfElements = 1 + tmp->m_nbrOfElements; // overtake number of elements
                m_root->m_children.at(0) = tmp;
            }
        }
    }

    HeapNode<T,C>* find( T val )
    {
        if( !m_root )
        {
            return nullptr;
        }
        else
        {
            return m_root->find(val);
        }
    };

    size_t size() const
    {
        if( m_root == nullptr )
            return 0;

        return m_root->m_nbrOfElements;
    }

    HeapNode<T,C>* m_root;
    HeapType m_type;

private:
    HeapNode<T,C>* createNode(T val)
    {
        if( m_type == HeapType::Max )
            return new HeapNodeMax<T,C>(val);
        else
            return new HeapNodeMin<T,C>(val);
    }
};


template <typename T, size_t C>
std::ostream& operator<<(std::ostream& os, const HeapNode<T,C>* node)
{
    os << "(" << node->m_value;

    for( const HeapNode<T,C>* e : node->m_children )
        if( e!= nullptr )
            os << e;
        else
            os << "(-)";

    os << ")";

    return os;
}

template <typename T, size_t C>
std::ostream& operator<<(std::ostream& os, const Heap<T,C>* heap)
{
    os << heap->m_root << std::endl;
    return os;
}


template <typename T>
class BSTNode
{

public:

    BSTNode<T>* m_left;
    BSTNode<T>* m_right;
    T m_val;
    size_t m_balanceFactor;

    BSTNode( const T& val) : m_left(nullptr), m_right(nullptr), m_val(val), m_balanceFactor(0)
    {
    }

    ~BSTNode()
    {
        if(m_left)
            delete m_left;
        if(m_right)
            delete m_right;
    }

    bool insert(const T& val)
    {
        bool insertResult;
        if(val == m_val)
        {
            // no two same values in the bst
            insertResult = false;
        }
        else
        {
            if (val < m_val)
                m_left = insertToNode(m_left, val, insertResult);
            else if (val > m_val)
                m_right = insertToNode(m_right, val, insertResult);
        }

        return insertResult;
    }

    bool search( const T& val ) const
    {
        if( val == m_val )
            return true;

        // run down left branch
        if( val < m_val && m_left != nullptr )
            return m_left->search(val);

        // run down right branch
        if( val > m_val && m_right != nullptr )
            return m_right->search(val);

        return false;
    }

    BSTNode<T>* remove( const T& val, bool& removeResult )
    {
        if( val == m_val )
        {
            removeResult = true;
            return removeMySelf();
        }
        else if( m_right != nullptr && val > m_val )
        {
            BSTNode<T>* newRight = m_right->remove(val, removeResult);
            if( newRight == nullptr )
                delete m_right;

            m_right = newRight;
        }
        else if( m_left != nullptr && val < m_val )
        {
            BSTNode<T>* newLeft = m_left->remove(val, removeResult);
            if( newLeft == nullptr )
                delete m_left;

            m_left = newLeft;
        }

        return this;
    }

    enum class NodeType
    {
        NoNode,
        NoChild,
        OnlyLeftChild,
        OnlyRightChild,
        TwoChildren
    };

    static NodeType getNodeType( const BSTNode<T>* node )
    {
        if( node == nullptr )
            return NodeType::NoNode;

        if( node->m_left == nullptr && node->m_right == nullptr )
            return NodeType::NoChild;

        if( node->m_left != nullptr && node->m_right == nullptr )
            return NodeType::OnlyLeftChild;

        if( node->m_left == nullptr && node->m_right != nullptr )
            return NodeType::OnlyRightChild;

        if( node->m_left != nullptr && node->m_right != nullptr )
            return NodeType::TwoChildren;

        return NodeType::NoNode;
    }

    BSTNode<T>* getLeftMostChild()
    {
        BSTNode<T>* lChild = this;
        while(lChild->m_left != nullptr)
            lChild = lChild->m_left;

        return lChild;
    }

    void deleteButNotChildren()
    {
        this->m_right = nullptr;
        this->m_left = nullptr;
        delete this;
    }

    size_t height( ) const
    {
        size_t leftHeight = 0;
        size_t rightHeight = 0;

        if( m_left != nullptr )
            leftHeight = m_left->height();

        if( m_right != nullptr )
            rightHeight = m_right->height();

        return 1 + std::max( leftHeight, rightHeight );
    }

    size_t updateBalanceFactors( )
    {
        size_t leftHeight = 0;
        size_t rightHeight = 0;

        if( m_left != nullptr )
            leftHeight = m_left->updateBalanceFactors();

        if( m_right != nullptr )
            rightHeight = m_right->updateBalanceFactors();

        m_balanceFactor = rightHeight - leftHeight;

        return 1 + std::max( leftHeight, rightHeight );
    }

    bool rotateRight()
    {
        // check if rotation possible
        if( m_left == nullptr )
            return false;

        // create a new node insert in the right branch
        BSTNode<T>* downNode = new BSTNode<T>(m_val);
        downNode->m_right = m_right;
        downNode->m_left = m_left->m_right;
        m_right = downNode;

        // copy the value of the left child
        m_val = m_left->m_val;

        // connect the left branch to left grand child
        BSTNode<T>* keepForDeletion = m_left;
        m_left = m_left->m_left;

        keepForDeletion->deleteButNotChildren();

        return true;
    };

    bool rotateLeft()
    {
        // check if rotation possible
        if( m_right == nullptr )
            return false;

        // create a new node insert in the left branch
        BSTNode<T>* downNode = new BSTNode<T>(m_val);
        downNode->m_left = m_left;
        downNode->m_right = m_right->m_left;
        m_left = downNode;

        // copy the value of the right child
        m_val = m_right->m_val;

        // connect the right branch to right grand child
        BSTNode<T>* keepForDeletion = m_right;
        m_right = m_right->m_right;

        keepForDeletion->deleteButNotChildren();

        return true;
    };

    /**
     * Degenerating the tree to a ascending linked list
     */
    void flatten()
    {
        while( this->rotateRight() ) {}

        if( m_right != nullptr )
            m_right->flatten();
    }

    /**
     * Compute the maximum number of nodes in a tree
     * with a certain height.
     * @param height
     * @return Max number of elements.
     */
    static size_t getMaxNumberOfNodes(size_t height)
    {
        if( height == 0 )
            return 0;

        return -(1-std::pow(2,height));
    }

    /**
     * Compute the height of a compact bst for
     * a certain number of nodes.
     * @param nbrNodes Number of nodes.
     * @return Min height.
     */
    static size_t getCompactHeight(size_t nbrNodes)
    {
        double height = std::log2(nbrNodes+1.0); // -1 removed, as we start at height 1
        return (size_t) std::ceil(height);
    }

    /**
     * Compute the number of nodes in the lowest layer when
     * dealing with a compact tree. If the lowest layer is full,
     * zero is returned.
     * @param nbrNodes
     * @return
     */
    static size_t getCompactNbrNodesLowestLayer(size_t nbrNodes)
    {
        if( nbrNodes == 0)
            return 0;

        size_t expectedHeight = getCompactHeight(nbrNodes);
        size_t possibleNodes = getMaxNumberOfNodes(expectedHeight);

        if( possibleNodes == nbrNodes )
            return 0;

        size_t freeNodesInLastLayer = possibleNodes - nbrNodes;
        size_t usedNodesInLastLayer = std::pow(2, expectedHeight - 1) - freeNodesInLastLayer;

        return usedNodesInLastLayer;
    }

    bool compare(const BSTNode<T>* node)
    {
        // "this" cannot be null
        if( node == nullptr )
            return false;

        if( node->m_val != m_val )
            return false;

        if( getNodeType(this) != getNodeType(node) )
            return false;

        if( m_left != nullptr )
            if( !m_left->compare(node->m_left) )
                return false;

        if( m_right != nullptr )
            if( !m_right->compare(node->m_right) )
                return false;

        return true;
    }

private:

    static BSTNode<T>* insertToNode(BSTNode<T>* node, const T& val, bool& insertResult)
    {
        if( node )
        {
            insertResult = node->insert(val);
        }
        else
        {
            // free spot -> add a new node
            node = new BSTNode<T>(val);
            insertResult = true;
        }

        return node;
    }

    BSTNode<T>* removeMySelf()
    {
        NodeType type = getNodeType(this);
        if( type == NodeType::NoChild )
        {
            return nullptr;
        }
        else if( type == NodeType::OnlyRightChild )
        {
            // keep for deletion
            BSTNode<T>* tmp = m_right;

            // overtake content and children of right child
            copyNode(m_right, this);

            tmp->deleteButNotChildren();
        }
        else if( type == NodeType::OnlyLeftChild )
        {
            // keep for deletion
            BSTNode<T>* tmp = m_left;

            // overtake content and children of left child
            copyNode(m_left, this);

            tmp->deleteButNotChildren();
        }
        else if( type == NodeType::TwoChildren )
        {
            // copy left most child of right brunch but keep this children
            m_val = m_right->getLeftMostChild()->m_val;

            // delete copied node beneath
            bool notUsed;
            m_right = m_right->remove(m_val,notUsed);
        }

        return this;
    }

    static void copyNode( const BSTNode<T>* src, BSTNode<T>* dst )
    {
        dst->m_val = src->m_val;
        dst->m_left = src->m_left;
        dst->m_right = src->m_right;
    }
};

template <typename T>
class BST
{
public:
    BST() : m_root(nullptr), m_size(0)
    {
    }

    ~BST()
    {
        if( m_root )
            delete m_root;
    }

    bool insert( const T& val )
    {
        bool insertOk;

        if( m_root == nullptr )
        {
            m_root = new BSTNode<T>(val);
            m_size = 1;
            insertOk = true;
        }
        else
        {
            insertOk = m_root->insert(val);
            if(insertOk)
                m_size++;
        }

        return insertOk;
    }

    bool find( const T& val ) const
    {
        if( m_root == nullptr )
            return false;

        return m_root->search(val);
    }

    bool remove( const T& val )
    {
        bool removeRes = false;

        if( m_root != nullptr )
        {
            m_root = m_root->remove(val, removeRes);
            if( removeRes )
                m_size--;
        }

        return removeRes;
    }

    size_t size() const
    {
        return m_size;
    }

    size_t height() const
    {
        size_t s = 0;

        if( m_root != nullptr )
            s = m_root->height();

        return s;
    }

    void computeBalanceFactors()
    {
        if( m_root != nullptr )
            m_root->updateBalanceFactors();
    }

    void print() const
    {
        std::cout << "Number of elements: " << this->size() << std::endl;
        print("", m_root, false);
    }

    /**
     * Form a balanced tree
     * Info: DAY-STOUT-WARREN algorithm http://www.smunlisted.com/day-stout-warren-dsw-algorithm.html
     */
    void balance()
    {
        if( m_root == nullptr )
            return;

        // step 1 -> make a linked list
        m_root->flatten();

        // step 2 -> form balanced tree
        size_t nNodes = size();
        size_t nNodesLowestLayer = BSTNode<T>::getCompactNbrNodesLowestLayer(nNodes);

        // left rotations on first nNodesLowestLayer on every odd node from root
        rotateBackbone(m_root, nNodesLowestLayer);

        size_t times = nNodes - nNodesLowestLayer;
        while( times > 1 )
        {
            times = times / 2;
            rotateBackbone(m_root, times);
        }
    }

    bool compare(const BST<T>* bst)
    {
        if( m_root == nullptr && bst->m_root == nullptr )
            return true;

        if( m_root == nullptr && bst->m_root != nullptr )
            return false;

        return m_root->compare(bst->m_root);
    }

    bool isBalanced()
    {
        size_t actualHeight = height();
        size_t balancedHeight = BSTNode<T>::getCompactHeight(size());
        return actualHeight == balancedHeight;
    }

    BSTNode<T>* m_root;
    size_t m_size;

private:

    // from java implementation: https://stackoverflow.com/a/42449385/2631225
    void print(const std::string& prefix, const BSTNode<T>* node, bool isLeft) const
    {
        if (node != nullptr)
        {
            std::cout << prefix;

            std::cout << (isLeft ? "├──" : "└──");

            // print the value of the node
            std::cout << node->m_val << std::endl;

            // enter the next tree level
            print(prefix + (isLeft ? "│   " : "    "), node->m_left, true);
            print(prefix + (isLeft ? "│   " : "    "), node->m_right, false);
        }
    }

    void rotateBackbone(BSTNode<T>* root, size_t count)
    {
        BSTNode<T>* rotNode = root;
        for( size_t k = 0; k < count; k++ )
        {
            BSTNode<T>* nextRot = rotNode->m_right->m_right;
            rotNode->rotateLeft();
            rotNode = nextRot;
        }
    }

};


#endif //MY_DATASTRUCTURES_H

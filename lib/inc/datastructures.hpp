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

template <typename T, size_t C>
class HeapNode
{
public:
    HeapNode( T value ) : m_children(C), m_value(value)
    {
        std::fill(m_children.begin(), m_children.end(), nullptr );
    }

    ~HeapNode()
    {
        for( HeapNode<T,C>* e : m_children )
            delete e;
    }

    virtual bool compareFunction(T a, T b)
    {
        return a > b;
    }

    virtual size_t childToReplace()
    {
        // find the index of the child to replace. in max heap,
        // it is the child with the smallest value.

        size_t smallest_idx = 0;
        T smallestVal = std::numeric_limits<T>::max();

        for( size_t i = 0; i < m_children.size(); i++ )
        {
            T cVal = m_children.at(i)->m_value;
            if( cVal < smallestVal )
            {
                smallestVal = cVal;
                smallest_idx = i;
            }
        }

        return smallest_idx;
    }

    bool insert(T val)
    {
        // val must be smaller or equal to this node's value
        if( compareFunction(val, m_value) )
            return false;

        // if empty child place, put it there
        ssize_t newChildIdx = getEmptyChildIdx();
        if( newChildIdx >= 0 )
        {
            m_children[newChildIdx] = new HeapNode<T,C>(val);
            return true;
        }

        // all children spots are occupied
        // try to insert into every child
        for( HeapNode<T,C>* e : m_children )
        {
            if( e->insert(val) )
            {
                return true;
            }
        }


        // the value is bigger than any of the children ones, but smaller or equal
        // of the this node. replace the smallest child with the new value node.
        size_t smallest_idx = childToReplace();
        HeapNode<T,C>* newNode = new HeapNode<T,C>(val);
        newNode->m_children.at(0) = m_children.at(smallest_idx);
        m_children.at(smallest_idx) = newNode;

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
        // if val is bigger than the value in this node,
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

public:
    std::vector<HeapNode<T,C>*> m_children;
    T m_value;
};

template <typename T, size_t C>
class HeapNodeMin: public HeapNode<T,C>
{
public:
    HeapNodeMin( T value ) : HeapNode<T,C>(value)
    {
    }

    ~HeapNodeMin()
    {
    }

    bool compareFunction(T a, T b) override
    {
        return a < b;
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


enum class HeapType
{
    Max,
    Min
};

template <typename T, size_t C>
class Heap
{

public:
    Heap(HeapType heapType = HeapType::Max) : m_root(nullptr)
    {

    }

    ~Heap()
    {
        delete m_root;
    }

    HeapNode<T,C>* m_root;
    HeapType m_type;

    void insert(T val)
    {
        if( !m_root )
        {
            // first element -> root does not exist yet
            m_root = new HeapNode<T,C>(val);
        }
        else
        {
            // try to insert
            if(!(m_root->insert(val)))
            {
                // insert does not work ->
                // val will become the new root element.
                HeapNode<T,C>* tmp = m_root;
                m_root = new HeapNode<T,C>(val);
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
};

template <typename T, size_t C>
std::ostream& operator<<(std::ostream& os, const Heap<T,C>* heap)
{
    os << heap->m_root << std::endl;
    return os;
}

#endif //MY_DATASTRUCTURES_H

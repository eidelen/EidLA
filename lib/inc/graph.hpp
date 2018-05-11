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

#ifndef MY_GRAPH_H
#define MY_GRAPH_H

#include "matrix.hpp"
#include "exceptions.hpp"

#include <iostream>
#include <limits>

/**
 * This class encapsulates functions related to graphs.
 * The graph is represented by a matrix.
 */
class Graph
{

public:

    /**
     * Graph constructor.
     * @param graphMat The matrix representing the graph.
     */
    template <class R>
    Graph(const Matrix<R> &graphMat) : m_graph( graphMat )
    {
        if(!m_graph.isSquare())
        {
            std::cerr << "Invalid graph matrix. Matrix needs to be of type  m x m" << std::endl;
            throw squareMatException;
        }
    }

    /**
     * Get the graph represented by a matrix.
     * @return Graph as matrix
     */
    const Matrix<double>& getGraph() const
    {
        return m_graph;
    }

    /**
     * Checks if there is a direct edge from vertex v0 to v1
     * @param v0 Start vertex
     * @param v1 End vertex
     * @return True if direct edge exists
     */
    bool isDirectlyConnected( size_t v0, size_t v1 )
    {
        return isEdge( m_graph(v0,v1) );
    }

private:

    bool isEdge( double weight )
    {
        return weight > std::numeric_limits<double>::min();
    }

    Matrix<double> m_graph;
};



#endif //MY_GRAPH_H

# My own C++ linear algebra library

This is a C++ implementation of common linear algebra operations and data types.
The main motivation behind this project is to refresh and deepen my knowledge in numerics and matrix operations.
Further, I am working on several projects which are using OpenCV or Eigen for matrix operations.
At one point, I hope to replace these with my own implementation without loosing too much performance.

Each library function has corresponding unit tests.
Google's gtest was used as a testing library.


Usage
-----

The main type is the templated class Matrix. There exist different constructors, such

.. code:: cpp

    #include "matrix.hpp"

    // creates an empty 4 x 6 double matrix
    Matrix<double> m1 = Matrix<double>(4, 6);

    // creates 2 x 2 int matrix with content [1,2; 3,4]
    int  data[] = {1, 2, 3, 4};
    Matrix<int> m2 = Matrix<int>(2, 2, data);


Matrix elements can be accessed with

.. code:: cpp

    // assign the value 5 to element m = 0, n = 1
    m2(0,1) = 5;

    // read the value on position m = 1, n = 0
    int v10 = m2(1,0);

    // read entire row 1
    Matrix<int> row1 = m2.row(1);

    // write entire row 1 to row 0.
    m2.setRow(0, row1);


Assign special content

.. code:: cpp

    // creates a 3x3 identity matrix
    Matrix<int> ident3x3 = Matrix<int>::identity(3);

    // sets an existing matrix to identity
    m2.setToIdentity();

    // fills the entire matrix with the value "pi"
    auto piMat = Matrix<double>(2,4);
    piMat.fill(3.14);


Arithmetic functions

.. code:: cpp

    int  d1[] = {1, 2, 3, 4};
    Matrix<int> a = Matrix<int>(2, 2, d1);

    int  d2[] = {1, 2, 3, 4};
    Matrix<int> b = Matrix<int>(2, 2, d2);

    // add the elements of two matrices
    Matrix<int> sum = a + b;

    // subtract
    Matrix<int> sub = a - b;

    // matrix multiplication
    Matrix<int> prod = a * b;

    // matrix multiplication with scalar
    Matrix<int> scale = a * 5;


Matrix properties

.. code:: cpp

    auto mat = Matrix<int>(2,2);

    // get number of rows and columns
    mat->rows();
    mat->cols();

    // get matrix rank
    size_t rank = mat.getRank();

    // get matrix inverse
    bool invertable;
    Matrix<double> inv = mat.inverted(&invertable);


Matrix transformations

.. code:: cpp

    #include "decomposition.hpp"

    // Lower-Upper triangle decomposition
    Decomposition::LUResult luRes = Decomposition::luDecomposition(mat);
    Matrix<double> lowerTriangle = luRes.L;
    Matrix<double> upperTriangle = luRes.U;


    #include "transformation.hpp"

    // Echelon transformations (possible with list of basic operations)
    Matrix<double> echelon = Transformation::echelon(mat);
    Matrix<double> reducedEchelon = Transformation::reduced_echelon(mat);



License
-------

MIT license: To my understanding, you can do whatever you wish to do with the code. However, no warranty is given that the written code is correct.



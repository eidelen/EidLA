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

    // creates a 4 x 6 double matrix
    Matrix<double> m1 = Matrix<double>(4, 6);
    
    // creates 2 x 2 int matrix [1,2; 3,4]
    int  data[] = {1, 2, 3, 4};
    Matrix<int> m2 = Matrix<int>(2, 2, data);
    
The matrix elements can be accessed with

    // assign the value 5 to element m = 0, n = 1
    m2(0,1) = 5;
    
    // read the value on position m = 1, n = 0
    int v10 = m2(1,0);


Arithmetic functions

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
    
    
To be extended....

License
-------
MIT license: To my understanding, you can do whatever you wish to do with the code. However, no warranty is given that the written code is correct.



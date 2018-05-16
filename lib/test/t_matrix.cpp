#include <gtest/gtest.h>
#include "matrix.hpp"

TEST(Matrix, InitializationCheckSizes)
{
    size_t nRows     = 5;
    size_t nCols     = 10;
    size_t nElements = nRows * nCols;

    MatrixISP m(new Matrix<int>(nRows, nCols));

    ASSERT_EQ(m->rows(), nRows);
    ASSERT_EQ(m->cols(), nCols);
    ASSERT_EQ(m->getNbrOfElements(), nElements);
}

TEST(Matrix, ConstructorWithData)
{
    int  src[] = {1, 2, 3, 4};
    auto mat   = Matrix<int>(2, 2, src);

    ASSERT_EQ(mat(0, 0), 1);
    ASSERT_EQ(mat(0, 1), 2);
    ASSERT_EQ(mat(1, 0), 3);
    ASSERT_EQ(mat(1, 1), 4);
}

TEST(Matrix, CopyConstructor)
{
    auto orig = Matrix<int>(2, 2);
    orig.fill(2);
    orig(1, 1) = 5;

    auto cpy = orig;
    ASSERT_TRUE(orig.compare(cpy));

    cpy(0, 0) = 3;
    ASSERT_FALSE(orig.compare(cpy));

    orig(0, 0) = 3;
    ASSERT_TRUE(orig.compare(cpy));
}

TEST(Matrix, CopyConstructorInt2Double)
{
    auto in_int = Matrix<int>(2, 2);
    in_int.setToIdentity();

    auto copy_double = Matrix<double>(in_int);
    auto soll        = Matrix<double>::identity(2);

    ASSERT_TRUE(soll.compare(in_int));
}

TEST(Matrix, SetGetSingleValue)
{
    size_t      nRows      = 5;
    size_t      nCols      = 10;
    int         currentVal = 2;
    Matrix<int> mat(nRows, nCols);

    for (size_t m = 0; m < nRows; m++)
    {
        for (size_t n = 0; n < nCols; n++)
        {
            currentVal++;
            mat.setValue(m, n, currentVal);
            ASSERT_EQ(mat.getValue(m, n), currentVal);

            currentVal++; //operators
            mat(m, n) = currentVal;
            ASSERT_EQ(mat(m, n), currentVal);
        }
    }
}

TEST(Matrix, SetElementsToSameValue)
{
    size_t nRows = 5;
    size_t nCols = 10;
    int    val   = 99;

    MatrixISP mat(new Matrix<int>(nRows, nCols));
    mat->fill(val);

    for (size_t m = 0; m < nRows; m++)
        for (size_t n = 0; n < nCols; n++)
            ASSERT_EQ(mat->getValue(m, n), val);
}

TEST(Matrix, Identity)
{
    size_t      s    = 5;
    Matrix<int> eyeM = Matrix<int>::identity(s);

    for (size_t m = 0; m < s; m++)
        for (size_t n = 0; n < s; n++)
            if (n == m)
                ASSERT_EQ(eyeM(m, n), 1);
            else
                ASSERT_EQ(eyeM(m, n), 0);

    auto secondEye = Matrix<int>(s, s);
    secondEye.setToIdentity();
    ASSERT_TRUE(secondEye.compare(eyeM));
}

TEST(Matrix, Row)
{
    auto mat = Matrix<int>::identity(3);

    auto r3 = mat.row(2); // r3 -> 0,0,1

    auto s3 = Matrix<int>(1, 3);
    s3(0, 0) = 0;
    s3(0, 1) = 0;
    s3(0, 2) = 1;

    ASSERT_TRUE(s3.compare(r3));
}

TEST(Matrix, Column)
{
    auto mat = Matrix<int>::identity(3);

    auto r3 = mat.column(2); // r3 -> 0;0;1

    auto s3 = Matrix<int>(3, 1);
    s3(0, 0) = 0;
    s3(1, 0) = 0;
    s3(2, 0) = 1;

    ASSERT_TRUE(s3.compare(r3));
}

TEST(Matrix, Swap)
{
    auto mat = Matrix<int>::identity(3);
    mat.swapRows(1, 2);

    int  sollData[9] = {1, 0, 0, 0, 0, 1, 0, 1, 0};
    auto soll        = Matrix<int>(3, 3, sollData);

    ASSERT_TRUE(soll.compare(mat));
}

TEST(Matrix, SetRow)
{
    auto soll = Matrix<int>::identity(3);
    soll.swapRows(1, 2);

    auto mat = Matrix<int>::identity(3);
    auto r1  = mat.row(1);
    auto r2  = mat.row(2);

    mat.setRow(1, r2);
    mat.setRow(2, r1);

    ASSERT_TRUE(soll.compare(mat));
}

TEST(Matrix, SetColumn)
{
    int  inData[] = {1, 0, 1, 0};
    auto in       = Matrix<int>(2, 2, inData);

    auto soll = Matrix<int>(2, 2);
    soll.fill(1);

    auto cpy = in.column(0);
    in.setColumn(1, cpy);

    ASSERT_TRUE(soll.compare(in));
}

TEST(Matrix, Transpose)
{
    int  inData[6] = {1, 1, 1, 3, 3, 3};
    auto in        = Matrix<int>(2, 3, inData);

    int  sollData[6] = {1, 3, 1, 3, 1, 3};
    auto soll        = Matrix<int>(3, 2, sollData);

    auto inTranspose = in.transpose();

    ASSERT_TRUE(inTranspose.compare(soll));
}

TEST(Matrix, Diagonal)
{
    int  inData[9] = {1, 1, 1, 3, 3, 3, 5, 5, 5};
    auto in        = Matrix<int>(3, 3, inData);

    int  diagData[3] = {1, 3, 5};
    auto diagSoll    = Matrix<int>(3, 1, diagData);

    ASSERT_TRUE(diagSoll.compare(in.diagonal()));
}

TEST(Matrix, Sum)
{
    int  inData[9] = {1, 1, 1, 3, 3, 3, 5, 5, 5};
    auto in        = Matrix<int>(3, 3, inData);

    ASSERT_EQ(in.sum(), 27);
}

TEST(Matrix, AssignmentOperator)
{
    auto mat = Matrix<int>(3, 3);
    mat.setToIdentity();
    int* addressAllocOrigin = mat.data();

    auto newMat = Matrix<int>(2, 2);
    newMat.fill(5);

    mat                  = newMat;
    int* addressAllocNew = mat.data();

    ASSERT_TRUE(mat.compare(newMat));
    ASSERT_TRUE(addressAllocNew != addressAllocOrigin); // different number of elements -> new alloc -> different address
}

TEST(Matrix, AssignmentOperatorSameNbrElement)
{
    auto mat = Matrix<int>(3, 3);
    mat.setToIdentity();
    int* addrOriginal = mat.data();

    auto newMat = Matrix<int>(1, 9);
    newMat.fill(5);

    mat                   = newMat;
    int* addrAfterAssignt = mat.data();

    ASSERT_TRUE(mat.compare(newMat));
    ASSERT_TRUE(addrOriginal == addrAfterAssignt); // should be same, since same number of elements
}

TEST(Matrix, RemoveRow)
{
    // remove row 1

    auto mat = Matrix<int>::random(3,2,0,10);

    auto soll = Matrix<int>(2,2);
    soll.setRow(0, mat.row(0));
    soll.setRow(1, mat.row(2));

    mat.removeRow(1);

    ASSERT_TRUE(mat.compare(soll));
}

TEST(Matrix, RemoveColumn)
{
    // remove col 0

    auto mat = Matrix<int>::random(2,3,0,10);

    auto soll = Matrix<int>(2,2);
    soll.setColumn(0, mat.column(1));
    soll.setColumn(1, mat.column(2));

    mat.removeColumn(0);

    ASSERT_TRUE(mat.compare(soll));
}

TEST(Matrix, SubMatrix)
{
    int matData[] = {0,1,2,   3,4,5};
    auto mat = Matrix<int>(2,3, matData);

    int sollData[] = {1, 4};
    auto soll = Matrix<int>(2,1, sollData);
    auto subMat = mat.subMatrix(0,1,2,1);
    ASSERT_TRUE(soll.compare(subMat));

    int sollData2[] = {1,2, 4,5};
    auto soll2 = Matrix<int>(2,2, sollData2);
    auto subMat2 = mat.subMatrix(0,1,2,2);
    ASSERT_TRUE(soll2.compare(subMat2));

    int sollData3[] = {3};
    auto soll3 = Matrix<int>(1,1, sollData3);
    auto subMat3 = mat.subMatrix(1,0,1,1);
    ASSERT_TRUE(soll3.compare(subMat3));
}

TEST(Matrix, SetSubMatrix)
{
    int matData[] = {0, 1, 2, 3, 4, 5};
    auto mat = Matrix<int>(2, 3, matData);

    int submatData[] = {8, 9, 10, 11};
    auto sub = Matrix<int>(2, 2, submatData);

    int sollData[] = {0, 8, 9, 3, 10, 11};
    auto soll = Matrix<int>(2, 3, sollData);

    mat.setSubMatrix(0, 1, sub);

    ASSERT_TRUE(soll.compare(mat));

    Matrix<int> nextTest = Matrix<int>(2,3);
    nextTest.setSubMatrix(0,0, mat);

    ASSERT_TRUE(nextTest.compare(mat));
}

TEST(Matrix, Symmetric)
{
    auto mat = Matrix<int>(3,3);
    mat.fill(1);
    ASSERT_TRUE(mat.isSymmetric());

    mat(0,0) = 2;
    ASSERT_TRUE(mat.isSymmetric());

    mat(0,1) = 5;
    ASSERT_FALSE(mat.isSymmetric());
}

TEST(Matrix, Square)
{
    auto mat = Matrix<int>(3,3);
    ASSERT_TRUE(mat.isSquare());
    mat.removeColumn(0);
    ASSERT_FALSE(mat.isSquare());
}

TEST(Matrix, SerializeDeserialize)
{
    int data[] = {1,2,3,  4,5,6};
    Matrix<int> serialize = Matrix<int>(2, 3, data);
    std::string serialData = serialize.serialize();

    Matrix<int> deserialize(1,1);
    deserialize.deserialize(serialData);

    ASSERT_TRUE(serialize.compare(deserialize));
}


TEST(Matrix, SaveLoad)
{
    std::string filename = "tmpmat.mat";

    int data[] = {1,2,3,  4,5,6};
    Matrix<int> save = Matrix<int>(2, 3, data);
    ASSERT_TRUE(save.save(filename));

    Matrix<int> load(filename);

    ASSERT_TRUE(save.compare(load));
}

TEST(Matrix, SaveLoadDouble)
{
    std::string filename = "tmpmat.mat";

    double data[] = {1,2,3,  4,5,6};
    Matrix<double> save = Matrix<double>(2, 3, data);
    ASSERT_TRUE(save.save(filename));

    Matrix<double> load(filename);

    ASSERT_TRUE(save.compare(load));
}


#ifdef OPENCVEIDLA
TEST(MATRIX, FromOpenCV)
{
    double matDataD[] = { 1,2,  3,4,  5,6,  7,8};
    cv::Mat opencvMatD(2,4,CV_64F, matDataD); // double
    float matDataF[] = { 1,2,  3,4,  5,6,  7,8};
    cv::Mat opencvMatF(2,4,CV_32F, matDataF); // float
    int matDataI[] = { 1,2,  3,4,  5,6,  7,8};
    cv::Mat opencvMatI(2,4,CV_32S, matDataI); // int

    Matrix<double> shouldBeMatD(2,4, matDataD);
    Matrix<float> shouldBeMatF(2,4, matDataF);
    Matrix<int> shouldBeMatI(2,4, matDataI);

    Matrix<double> fromOpenCV( opencvMatD );
    ASSERT_TRUE( fromOpenCV.compare(shouldBeMatD) );

    Matrix<float> fromOpenCVF( opencvMatF );
    ASSERT_TRUE( fromOpenCVF.compare(shouldBeMatF) );

    Matrix<int> fromOpenCVI( opencvMatI );
    ASSERT_TRUE( fromOpenCVI.compare(shouldBeMatI) );
}

TEST(MATRIX, ToOpenCV)
{
    double matDataD[] = { 1,2,  3,4,  5,6,  7,8};
    Matrix<double> shouldBeMatD(2,4, matDataD);

    float matDataF[] = { 1,2,  3,4,  5,6,  7,8};
    Matrix<float> shouldBeMatF(2,4, matDataF);

    int matDataI[] = { 1,2,  3,4,  5,6,  7,8};
    Matrix<int> shouldBeMatI(2,4, matDataI);

    short matDataShort[] = { 1,2,  3,4,  5,6,  7,8};
    Matrix<short > shouldBeMatShort(2,4, matDataShort);

    cv::Mat cvMatD = shouldBeMatD.toOpenCVMatrix();
    cv::Mat cvMatF = shouldBeMatF.toOpenCVMatrix();
    cv::Mat cvMatI = shouldBeMatI.toOpenCVMatrix();
    cv::Mat cvMatSD = shouldBeMatShort.toOpenCVMatrix();

    // uses template specializations
    ASSERT_TRUE( Matrix<double>(cvMatD).compare(shouldBeMatD) );
    ASSERT_TRUE( Matrix<float>(cvMatF).compare(shouldBeMatF) );
    ASSERT_TRUE( Matrix<int>(cvMatI).compare(shouldBeMatI) );

    // uses general implementation -> double returned
    ASSERT_TRUE( Matrix<double>(cvMatSD).compare(shouldBeMatD) );
}

#endif // OPENCVEIDLA

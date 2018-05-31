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

#include <iostream>
#include <opencv2/highgui/highgui.hpp>

#include "matrix.hpp"

int main(int argc, char* argv[])
{
    // read example image with OpenCV
    cv::Mat rawImg_cv = cv::imread("ins.jpg", CV_LOAD_IMAGE_GRAYSCALE);
    if( !rawImg_cv.data )
    {
        std::cerr << "Cannot load image from file";
        return -1;
    }

    // transform cv-matrix to EidLA matrix
    Matrix<uchar> rawImg_char( rawImg_cv );
    Matrix<double> rawImg(rawImg_char);

    // SVD decomposition
    Decomposition::SVDResult deco = Decomposition::svd(rawImg);

    // check if decomposition is correct
    if( !rawImg.compare(deco.U * deco.S * deco.V.transpose(), true, 0.05 ) )
    {
        std::cerr << "SVD decompsition failed";
        return -2;
    }

    float reducing = 0.2;
    size_t mods = std::round(deco.S.cols() * (1.0 - reducing));

    Matrix<double> compressed = deco.U.subMatrix(0,0,rawImg.rows(),mods) * deco.S.subMatrix(0,0,mods,mods) * deco.V.transpose().subMatrix(0,0,mods,rawImg.cols());
    // correct values
    for( size_t i = 0; i < compressed.getNbrOfElements(); i++ )
        if( compressed.data()[i] > 255.0 )
            compressed.data()[i] = 255;
        else if( compressed.data()[i] < 0.0 )
            compressed.data()[i] = 0;

    Matrix<uchar> compressed_uchar( compressed );

    std::cout << compressed_uchar;

    cv::Mat compressed_cv = compressed_uchar.toOpenCVMatrix();

    cv::imshow("compressed", compressed_cv);
    cv::waitKey(0);
/*

    TEST(Decomposition, SVDCompression)
    {
        // third line is almost linear combination of both upper
        // therefore, using first two singular values should lead to a accurate representation of the matrix.
        double data[] = {1, 0, 0,   0, 1, 0,   1, 1, 0.01}; // third line is almost linear combination of both upper
        Matrix<double> mat(3,3, data);

        Decomposition::SVDResult res = Decomposition::svd(mat);

        ASSERT_TRUE( res.U.isOrthogonal(0.05) );
        ASSERT_TRUE( res.V.isOrthogonal(0.05) );
        ASSERT_TRUE( mat.compare(res.U * res.S * res.V.transpose(), true, 0.05 ) );

        auto compressed = res.U.subMatrix(0,0,3,2) * res.S.subMatrix(0,0,2,2) * res.V.transpose().subMatrix(0,0,2,3);

        ASSERT_TRUE( mat.compare(compressed, true, 0.01));
    }



    if(! image.data )                              // Check for invalid input
    {
        cout <<  "Could not open or find the image" << std::endl ;
        return -1;
    }*/

    return 0;
}

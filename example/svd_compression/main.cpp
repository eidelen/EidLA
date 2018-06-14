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
    // read image with OpenCV
    cv::Mat rawImg_cv = cv::imread("ins.jpg", CV_LOAD_IMAGE_GRAYSCALE);
    if( !rawImg_cv.data )
    {
        std::cerr << "Cannot load image from file";
        return -1;
    }

    // transform cv-matrix to EidLA matrix
    Matrix<uchar> rawImg( rawImg_cv );

    // normalize the image (transpose -> needs to be higher than wide)
    double mean, scale;
    Matrix<double> normalized(rawImg.normalize(mean,scale).transpose());

    // SVD decomposition
    Decomposition::SVDResult deco = Decomposition::svd(normalized);

    // as the SVD decomposition may take several hours, it makes sense
    // to keep the resulting matrices
    deco.S.save("S.mat"); deco.U.save("U.mat"); deco.V.save("V.mat");

    // generate several compressed images, by increasing the number of mods used
    // for the reconstruction.
    for( size_t mods = 1; mods < deco.S.cols(); mods++)
    {
        // create compressed image
        Matrix<double> compressed =
                deco.U.subMatrix(0, 0, normalized.rows(), mods) * deco.S.subMatrix(0, 0, mods, mods) *
                deco.V.transpose().subMatrix(0, 0, mods, normalized.cols());

        // denormalize the compressed image
        compressed = compressed.denormalize(mean, scale);

        // correct pixel values if out of range
        for (size_t i = 0; i < compressed.getNbrOfElements(); i++)
            if (compressed.data()[i] > 255.0)
                compressed.data()[i] = 255;
            else if (compressed.data()[i] < 0.0)
                compressed.data()[i] = 0;

        // compute compression rate
        size_t compressNumberOfPixels = normalized.rows()*mods + mods*mods + mods*normalized.cols();
        double compressionRatio = static_cast<double>(compressNumberOfPixels)/static_cast<double>(normalized.getNbrOfElements());

        // generate file name
        std::string wName("m_");
        wName.append(std::to_string(mods));
        wName.append("_red_");
        wName.append(std::to_string(compressionRatio));
        wName.append(".png");

        // save the image
        cv::imwrite(wName, compressed.transpose().toOpenCVMatrix());
    }

    return 0;
}

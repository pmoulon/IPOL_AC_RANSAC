// Copyright (c) 2007-2011 libmv authors.
// 
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to
// deal in the Software without restriction, including without limitation the
// rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
// sell copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
// 
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
// IN THE SOFTWARE.

#include <vector>
#include "homography_model.hpp"
#include "libNumerics/homography.h"
#include "CppUnitLite/TestHarness.h"

#define EXPECT_MATRIX_NEAR(a, b, tolerance) \
do { \
  bool dims_match = (a.nrow() == b.nrow()) && (a.ncol() == b.ncol()); \
  CHECK_EQUAL(a.nrow(),b.nrow()); \
  CHECK_EQUAL(a.ncol(),b.ncol()); \
  if (dims_match) { \
    for (int r = 0; r < a.nrow(); ++r) { \
      for (int c = 0; c < a.ncol(); ++c) { \
        DOUBLES_EQUAL(a(r, c), b(r, c), tolerance); \
      } \
    } \
  } \
} while(false);

TEST(HomographyKernelTest, Fitting) {
  //------------------------//
  //-- Unit test Objective--//
  //------------------------//
  // Given two points set linked by a known homography :
  // -=> Compute the homography trough sampled points
  // -=> Assert that estimated H is good and residuals are near 0 to a given epsilon.
  //--
  //--

  typedef libNumerics::matrix<double> Mat;

  // Define a few homographies.
  std::vector<Mat> H_gt(3, Mat::zeros(3,3));

  H_gt[0] = Mat::eye(3);
  double matGTd1[] = { 1, 0, -4,  0, 1, 5,  0, 0, 1};
  H_gt[1].read(matGTd1);
  double matGTd2[] = { 1, 0, 3,  0, 1, -6,  0, 0, 1};
  H_gt[2].read(matGTd2);

  // Define a set of points.
  Mat x(2, 9);
  double points[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8,
                      0,.4,.6, 0, 1, 2, 0, 1, 2};
  x.read(points);

  for (int i = 0; i < H_gt.size(); ++i) {

    bool bFound = false;
    // Transform points by the ground truth homography.
    Mat y=x;
    for(int k = 0; k < 9; ++k)
    {
        libNumerics::TransformH(H_gt[i], y(0,k), y(1,k) );
    }
    std::vector<Match> m(9);
    for(int i=0; i<9; i++)
        m[i] = Match(x(0,i), x(1,i), y(0,i), y(1,i));

    orsa::HomographyModel kernel(m);

    //-- Fit a model and check re-projection error.
    int samples_[5]={0,1,2,3,4};
    std::vector<int> samples(samples_,samples_+5);
    for (int j = 4; samples.size() < x.ncol(); samples.push_back(j++)) {
      Mat H(3,3);
      if(kernel.ComputeModel(samples, &H))
      {
        EXPECT_MATRIX_NEAR(H_gt[i], H, 5e-7);
        bFound = true;
      }
    }
    CHECK(bFound);
  }
}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */

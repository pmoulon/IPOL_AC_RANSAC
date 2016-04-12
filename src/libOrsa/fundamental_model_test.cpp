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

#include <iostream>
#include "CppUnitLite/TestHarness.h"
#include "libOrsa/fundamental_model.hpp"

typedef libNumerics::matrix<double> Mat;

// Check the properties of a fundamental matrix:
//
//   1. The determinant is 0 (rank deficient)
//   2. The condition x'T*F*x = 0 is satisfied to precision.
//
bool ExpectFundamentalProperties(const Mat &F,
                                 const Mat &ptsA, int w1, int h1,
                                 const Mat &ptsB, int w2, int h2,
                                 double precision) {
  bool bOk = true;
  bOk &= F.det() < precision;
  std::cout << std::endl << F << std::endl;
  orsa::FundamentalModel model(ptsA, w1, h1, ptsB, w2, h2);
  assert(ptsA.ncol() == ptsB.ncol());
  for (int i = 0; i < ptsA.ncol(); ++i) {
    double residual = model.Error(F,i);
    bOk &= residual < precision;
  }
  return bOk;
}

// Check the fundamental fitting:
//
//   1. Estimate the fundamental matrix
//   2. Check epipolar distance.
//
template <class Kernel>
bool ExpectKernelProperties(const Mat &x1, int w1, int h1,
                            const Mat &x2, int w2, int h2) {
  bool bOk = true;
  orsa::OrsaModel* kernel = new Kernel(x1, w1, h1, x2, w2, h2);
  std::vector<int> samples;
  for (int i = 0; i < x1.ncol(); ++i) {
    samples.push_back(i);
  }
  std::vector<Mat> vec_F;
  kernel->Fit(samples, &vec_F);
  for (size_t i = 0; i  < vec_F.size(); ++i)
  {
    bOk &= ExpectFundamentalProperties(vec_F[i], x1, w1, h1, x2, w2, h2, 1e-8);
  }  
  delete kernel;
  return bOk;
}

TEST(SevenPointTest, EasyCase) {
  Mat x1(2, 7), x2(2, 7);
  double points[] = { 0, 0, 0, 1, 1, 1, 2,
                      0, 1, 2, 0, 1, 2, 0};
  x1.read(points);
  double points2[] = { 0, 0, 0, 1, 1, 1, 2,
                       1, 2, 3, 1, 2, 3, 1};
  x2.read(points2);
  typedef orsa::FundamentalModel Model;
  CHECK(ExpectKernelProperties<Model>(x1, 3, 3, x2, 3, 3));
}

TEST(SevenPointTest, RealCorrespondences) {
  Mat x1(2, 7), x2(2, 7);
  double points[] = { 723, 1091, 1691, 447,  971, 1903, 1483,
                      887,  699,  811, 635,   91,  447, 1555};
  x1.read(points);
  double points2[] = { 1251, 1603, 2067, 787, 1355, 2163, 1875,
                       1243,  923, 1031, 484,  363,  743, 1715};
  x2.read(points2);

  typedef orsa::FundamentalModel Model;
  CHECK(ExpectKernelProperties<Model>(x1, 1500, 1800, x2, 2200, 1800));
}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */

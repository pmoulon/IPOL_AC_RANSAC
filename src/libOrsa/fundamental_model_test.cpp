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
                                 const std::vector<Match>& m,
                                 double precision) {
  bool bOk = true;
  bOk &= F.det() < precision;
  std::cout << std::endl << F << std::endl;
  orsa::FundamentalModel model(m);
  for (int i = 0; i < model.NbData(); ++i) {
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
bool ExpectKernelProperties(const std::vector<Match> &m) {
  bool bOk = true;
  orsa::ModelEstimator* kernel = new Kernel(m);
  std::vector<int> samples;
  for (size_t i = 0; i < m.size(); ++i) {
    samples.push_back(static_cast<int>(i));
  }
  std::vector<Mat> vec_F;
  kernel->Fit(samples, &vec_F);
  for (size_t i = 0; i  < vec_F.size(); ++i)
  {
    bOk &= ExpectFundamentalProperties(vec_F[i], m, 1e-8);
  }  
  delete kernel;
  return bOk;
}

TEST(SevenPointTest, EasyCase) {
  double points1[2*7] = { 0, 0, 0, 1, 1, 1, 2,
                          0, 1, 2, 0, 1, 2, 0};
  double points2[2*7] = { 0, 0, 0, 1, 1, 1, 2,
                          1, 2, 3, 1, 2, 3, 1};
  std::vector<Match> m(7);
  for(int i=0; i<7; i++)
      m[i] = Match(points1[i], points1[i+7], points2[i], points2[i+7]);
  typedef orsa::FundamentalModel Model;
  CHECK(ExpectKernelProperties<Model>(m));
}

TEST(SevenPointTest, RealCorrespondences) {
  double points1[2*7] = { 723, 1091, 1691, 447,  971, 1903, 1483,
                          887,  699,  811, 635,   91,  447, 1555};
  double points2[2*7] = { 1251, 1603, 2067, 787, 1355, 2163, 1875,
                       1243,  923, 1031, 484,  363,  743, 1715};

  std::vector<Match> m(7);
  for(int i=0; i<7; i++)
      m[i] = Match(points1[i], points1[i+7], points2[i], points2[i+7]);

  typedef orsa::FundamentalModel Model;
  CHECK(ExpectKernelProperties<Model>(m));
}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */

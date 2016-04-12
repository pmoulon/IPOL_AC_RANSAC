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
#include <algorithm>

#include "libOrsa/fundamental_model.hpp"
#include "CppUnitLite/TestHarness.h"

typedef libNumerics::matrix<double> Mat;

// Unit test :
//-  It considers 16 points and an outlier.
//-  It asserts that the F matrix and inliers number are correct.
// F should be similar to:
// 0, -a,  -b,
// a,  0,  -c,
// b,  c,   0
//
// Note :
// This test use not normalized coordinates.
TEST(Fundamental_Orsa, Test)
{
  const int n = 16;
  Mat x1(2,n);
  double points[] = {
    0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4,   5,
    0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2,   5};
  x1.read(points);

  Mat x2 = x1;
  for (int i = 0; i < n; ++i) {
    x2(0, i) += i % 2; // Multiple horizontal disparities.
  }
  // The outlier has vertical disparity.
  x2(0, n - 1) = 10;
  x2(1, n - 1) = 10;

  /// Create the Kernel (Model estimator and tester)
  orsa::FundamentalModel model(x1, 5, 5, x2, 5, 5);
  Mat F(3,3);
  std::vector<int> vec_inliers;
  /// Orsa (robust estimation routine)
  model.orsa(vec_inliers, 100, NULL, &F, true);

  const double expectedPrecision = 1e-8;
  const double & ep = expectedPrecision;
  DOUBLES_EQUAL(0.0, F(0,0), ep);
  DOUBLES_EQUAL(0.0, F(1,1), ep);
  DOUBLES_EQUAL(0.0, F(2,2), ep);
  DOUBLES_EQUAL(F(0,1), -F(1,0), ep);
  DOUBLES_EQUAL(F(0,2), -F(2,0), ep);
  DOUBLES_EQUAL(F(1,2), -F(2,1), ep);
  CHECK_EQUAL(n - 1, vec_inliers.size());
  CHECK_EQUAL(true, std::find(vec_inliers.begin(), vec_inliers.end(), 15) == vec_inliers.end());
}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */

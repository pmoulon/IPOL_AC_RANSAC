/**
 * @file cc_math_svd_test.cpp
 * @brief SVD test
 * @author Pascal Monasse
 * 
 * Copyright (c) 2010 Pascal Monasse
 * All rights reserved.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "numerics.h"
#include "CppUnitLite/TestHarness.h"

static const int NTESTS=100;
static const double EPS=1E-5;

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

static double max(const libNumerics::matrix<double>& A) {
    double m = A(0,0);
    for(int i=0; i<A.nrow(); i++)
        for(int j=0; j<A.ncol(); j++) {
            double abs = fabs(A(i,j));
            if(abs>m)
                m = abs;
        }
    return m;
}

TEST(SVDTest, Decomposition) {
  //------------------------//
  //-- Unit test Objective--//
  //------------------------//
  // Test SVD decomposition
  //--
  //--

    for(int k=0; k<NTESTS; k++) {
        int m=1+rand()%50;
        int n=1+rand()%50;
        libNumerics::matrix<double> A(m,n);
        for(int i=0; i<m; i++)
            for(int j=0; j<n; j++)
                A(i,j) = (rand()/(double)RAND_MAX)*50-10;
        libNumerics::SVD svd(A);
        libNumerics::matrix<double> B=svd.compose();
        double err = max(A-B);
        std::cout << k << ": " << m << 'x' << n << "  " << err << std::endl;
        if(err > EPS) {
            std::cout << "A=" << A << std::endl;
            std::cout << "U=" << svd.U << std::endl;
            std::cout << "D=" << svd.D << std::endl;
            std::cout << "V=" << svd.V << std::endl;
        }
        EXPECT_MATRIX_NEAR(A, B, EPS);
    }
}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */

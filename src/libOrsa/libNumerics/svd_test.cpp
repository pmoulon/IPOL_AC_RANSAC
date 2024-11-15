/**
 * @file svd_test.cpp
 * @brief SVD test
 * @author Pascal Monasse
 * 
 * Copyright (c) 2010, 2024 Pascal Monasse
 */

#include "numerics.h"
#include "CppUnitLite/TestHarness.h"

static const int NTESTS=100;
static const double EPS=1E-5;
static const int SIZE_MAX=50;

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
        int m=1+rand()%SIZE_MAX;
        int n=1+rand()%SIZE_MAX;
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

        libNumerics::SVD svd2(A, libNumerics::SVD::compact);
        B = svd2.compose();
        if(! err > EPS) {
            std::cout << "U=" << svd2.U << std::endl;
            std::cout << "D=" << svd2.D << std::endl;
            std::cout << "V=" << svd2.V << std::endl;
            std::cout << "compact -> FAILURE" << std::endl;
        }
        std::cout << k << ": (compact) " << m<<'x'<<n << "  " << err<<std::endl;
        EXPECT_MATRIX_NEAR(A, B, EPS);
    }
}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */

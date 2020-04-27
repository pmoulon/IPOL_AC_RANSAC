/**
 * @file homography_model.cpp
 * @brief Homography matrix model
 * @author Pascal Monasse, Pierre Moulon
 * 
 * Copyright (c) 2011,2020 Pascal Monasse
 * Copyright (c) 2011 Pierre Moulon
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

#include "homography_model.hpp"
#include "conditioning.hpp"
#include "libNumerics/numerics.h"

#include <algorithm>
#include <iostream>
#include <limits>

#define _USE_MATH_DEFINES // For Windows (M_PI)
#include <math.h>

namespace orsa {

/// Minimal inverse condition value for a valid homography matrix
static const double ICOND_MIN = 0.1;

/// Return min/max/mean of row \a i of matrix \a A in \a m, \a M and \a avg.
void MinMaxAvgRow(const ModelEstimator::Mat& A, int i,
                  double& m, double& M, double& avg) {
  m = M = avg = A(i,0);
  const int n = A.ncol();
  for(int j=1; j<n; j++) {
    double a = A(i,j);
    if(m > a) m = a;
    if(M < a) M = a;
    avg += a;
  }
  avg /= static_cast<double>(n);
}

/// Constructor.
HomographyModel::HomographyModel(const Mat &x1, const Mat &x2, bool symError)
: ModelEstimator(x1, x2, symError), N1_(3,3), N2_(3,3) {
  double stat[12];
  MinMaxAvgRow(x1, 0, stat[0], stat[1],  stat[2]);
  MinMaxAvgRow(x1, 1, stat[3], stat[4],  stat[5]);
  MinMaxAvgRow(x2, 0, stat[6], stat[7],  stat[8]);
  MinMaxAvgRow(x2, 1, stat[9], stat[10], stat[11]);
  // Normalize both images by same factor, as our thresholds for SVD are based
  // on zoom around 1
  double w = std::max(stat[1]-stat[0], stat[7] -stat[6]);
  double h = std::max(stat[4]-stat[3], stat[10]-stat[9]);
  double dNorm = 1.0 / sqrt( static_cast<double>(w*h) );
  N1_ = Mat::eye(3); N2_ = Mat::eye(3);
  N1_(0,0) = N1_(1,1) = N2_(0,0) = N2_(1,1) = dNorm;
  N1_(0,2) = -stat[2]*dNorm; N1_(1,2) = -stat[5]*dNorm;
  N2_(0,2) = -stat[8]*dNorm; N2_(1,2) = -stat[11]*dNorm;
}

/// Unnormalize a given model (from normalized to image space).
void HomographyModel::Unnormalize(Model *model) const  {
  UnnormalizerI::Unnormalize(N1_, N2_, model);
}

/// Square of number.
inline double sqr(double x) {
  return (x*x);
}

/// If H is not orientation preserving at the point, the error is infinite.
/// For this test, it is assumed that det(H)>0.
/// \param[in] H The homography matrix.
/// \param[in] index The correspondence index
/// \param[out] side In which image is the error measured?
/// \return The square reprojection error.
double HomographyModel::Error(const Model &H, int index, int* side) const {
  double err = std::numeric_limits<double>::infinity();
  if(side) *side=1;
  libNumerics::vector<double> x(3);
  // Transfer error in image 2
  x(0) = x1_(0,index); x(1) = x1_(1,index); x(2) = 1.0;
  x = H * x;
  if(x(2)>0) {
    x /= x(2);
    err = sqr(x2_(0,index)-x(0)) + sqr(x2_(1,index)-x(1));
  }
  // Transfer error in image 1
  if(symError) { // ... but only if requested
    x(0) = x2_(0,index); x(1) = x2_(1,index); x(2) = 1.0;
    x = H.inv() * x;
    if(x(2)>0) {
      x /= x(2);
      double err1 = sqr(x1_(0,index)-x(0)) + sqr(x1_(1,index)-x(1));
      if(err1>err) { // Keep worse error
        err=err1;
        if(side) *side=0;
      }
    }
  }
  return err;
}

/// 2D homography estimation from point correspondences.
void HomographyModel::Fit(const std::vector<int> &indices,
                          std::vector<Model> *H) const {
  if(4 > indices.size())
    return;
  const int n = static_cast<int>( indices.size() );
  Mat A = Mat::zeros(n*2,9);
  for (int i = 0; i < n; ++i) {
    int j = indices[i];
    libNumerics::vector<double> x1(3), x2(3);
    x1(0)=x1_(0,j); x1(1)=x1_(1,j); x1(2)=1;
    x2(0)=x2_(0,j); x2(1)=x2_(1,j); x2(2)=1;
    x1 = N1_*x1;
    x2 = N2_*x2;
    j = 2*i;
    A(j,0) = x1(0);
    A(j,1) = x1(1);
    A(j,2) = 1.0;
    A(j,6) = -x2(0) * x1(0);
    A(j,7) = -x2(0) * x1(1);
    A(j,8) = -x2(0);

    ++j;
    A(j,3) = x1(0);
    A(j,4) = x1(1);
    A(j,5) = 1.0;
    A(j,6) = -x2(1) * x1(0);
    A(j,7) = -x2(1) * x1(1);
    A(j,8) = -x2(1);
  }

  libNumerics::vector<double> vecNullspace(9);
  if( libNumerics::SVD::Nullspace(A,&vecNullspace) )
  {
    Mat M(3,3);
    M.read(vecNullspace);
    if(libNumerics::SVD::InvCond(M) < ICOND_MIN)
      return;
    Unnormalize(&M);
    if(M.det() < 0)
      M = -M;
    M /= M(2,2);

    if(IsOrientationPreserving(indices,M) )
      H->push_back(M);
  }
}

/// Is \a H orientation preserving near match points of index in \a indices?
/// It is assumed that det(H)>0, otherwise the sign tests should be reversed.
bool HomographyModel::IsOrientationPreserving(const std::vector<int>&indices,
                                              const Mat& H) const {
  libNumerics::matrix<double> h=H.row(2);
  std::vector<int>::const_iterator it=indices.begin(), end=indices.end();
  for(; it != end; ++it)
    if(h(0)*x2_(0,*it)+h(1)*x2_(1,*it)+h(2) <= 0)
      return false;
  return true;
}

}  // namespace orsa

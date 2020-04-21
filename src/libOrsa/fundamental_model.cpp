/**
 * @file fundamental_model.cpp
 * @brief Fundamental matrix model
 * @author Pascal Monasse, Pierre Moulon
 * 
 * Copyright (c) 2011-2020 Pascal Monasse
 * Copyright (c) 2011-2015 Pierre Moulon
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

#define _USE_MATH_DEFINES //For Windows (M_PI)
#include "fundamental_model.hpp"
#include "libNumerics/cubicRoots.h"
#include "libNumerics/numerics.h"

typedef libNumerics::vector<double> Vec;

namespace orsa {

/// Min product of norms of two collinear vectors to be considered positive.
/// Below that threshold, one of the vectors is considered as too small to be
/// reliable. This is used in FundamentalModel::checkF.
static const double MIN_PRODUCT_NORMS=1e-5;

/// Constructor.
FundamentalModel::FundamentalModel(const Mat &x1, int w1, int h1,
                                   const Mat &x2, int w2, int h2,
                                   bool symmetricError)
: ModelEstimator(x1, w1, h1, x2, w2, h2, symmetricError) {}

/// Unnormalize a given model (from normalized to image space).
void FundamentalModel::Unnormalize(Model * model) const  {
  *model = N1_.t() * (*model) * N2_;
}

/// Build an nx9 matrix from point matches, where each row is equivalent to the
/// equation xT*F*x' = 0 for a single correspondence pair (x, x'). The domain of
/// the matrix is a 9 element vector corresponding to F. In other words, set up
/// the linear system
///   Af = 0,
/// where f is the F matrix as a 9-vector rather than a 3x3 matrix (row
/// major). If the points are well conditioned and there are 8 or more, then
/// the nullspace should be rank one. If the nullspace is two dimensional,
/// then the rank 2 constraint must be enforced to identify the appropriate F
///  matrix.
///
/// Do not resize the matrix A, assumed to have the appropriate size already.
static void EncodeEpipolarEquation(const ModelEstimator::Mat &x1,
                                   const ModelEstimator::Mat &x2,
                                   const std::vector<int> &indices,
                                   ModelEstimator::Mat *A) {
  for (size_t i = 0; i < indices.size(); ++i) {
    int j = indices[i];
    (*A)(i, 0) = x1(0,j) * x2(0,j);  // 0 represents x coords,
    (*A)(i, 1) = x1(0,j) * x2(1,j);  // 1 represents y coords.
    (*A)(i, 2) = x1(0,j);
    (*A)(i, 3) = x1(1,j) * x2(0,j);
    (*A)(i, 4) = x1(1,j) * x2(1,j);
    (*A)(i, 5) = x1(1,j);
    (*A)(i, 6) = x2(0,j);
    (*A)(i, 7) = x2(1,j);
    (*A)(i, 8) = 1.0;
  }
}

void FundamentalModel::Fit(const std::vector<int> &indices,
                           std::vector<Model> *Fs) const {
  assert(2 == x1_.nrow());
  assert(7 <= x1_.ncol());
  assert(x1_.nrow() == x2_.nrow());
  assert(x1_.ncol() == x2_.ncol());

  // Set up the homogeneous system Af = 0 from the equations x'T*F*x = 0.
  Mat A(indices.size(), 9);
  EncodeEpipolarEquation(x1_, x2_, indices, &A);

  if(indices.size() < 8) {
    algo7pt(A, Fs);
    checkF(indices, Fs);
  } else
    algo8pt(A, Fs);
}

/// \param F The fundamental matrix.
/// \param index The point correspondence.
/// \param side In which image is the error measured?
/// \return The square reprojection error.
double FundamentalModel::Error(const Model &F, int index, int* side) const {
  double xa = x1_(0,index), ya = x1_(1,index);
  double xb = x2_(0,index), yb = x2_(1,index);

  double a, b, c, d;
  // Transfer error in image 2
  if(side) *side=1;
  a = F(0,0) * xa + F(1,0) * ya + F(2,0);
  b = F(0,1) * xa + F(1,1) * ya + F(2,1);
  c = F(0,2) * xa + F(1,2) * ya + F(2,2);
  d = a*xb + b*yb + c;
  double err =  (d*d) / (a*a + b*b);
  // Transfer error in image 1
  if(symError) { // ... but only if requested
    a = F(0,0) * xb + F(0,1) * yb + F(0,2);
    b = F(1,0) * xb + F(1,1) * yb + F(1,2);
    c = F(2,0) * xb + F(2,1) * yb + F(2,2);
    d = a*xa + b*ya + c;
    double err1 =  (d*d) / (a*a + b*b);
    if(err1>err) {
      err = err1;
      if(side) *side=0;
    }
  }
  return err;
}

/// Determinant of 3x3 matrix expressed by its columns (v1 v2 v3).
inline double det3(const Vec& v1, const Vec& v2, const Vec& v3){
  return v1(0)*(v2(1)*v3(2)-v2(2)*v3(1))-
         v1(1)*(v2(0)*v3(2)-v2(2)*v3(0))+
         v1(2)*(v2(0)*v3(1)-v2(1)*v3(0));
}

/// 7-point algorithm.
/// \param A Matrix such that f is solution to Af=0
/// \param Fs Output vector of found F matrices (up to 3)
void FundamentalModel::algo7pt(const Mat& A, std::vector<Mat> *Fs) const {
  // Find two F matrices in the nullspace of A
  Mat F1(3,3), F2(3,3);
  libNumerics::SVD::Nullspace2_Remap33(A,F1,F2);
  F2 -= F1;

  // Use condition det(F)=0 to determine F: solve det(F1 + t*F2) = 0 for t
  Vec c1=F1.col(0), c2=F1.col(1), c3=F1.col(2);
  Vec d1=F2.col(0), d2=F2.col(1), d3=F2.col(2);
  double P[4] = {
    det3(c1,c2,c3),
    det3(d1,c2,c3)+det3(c1,d2,c3)+det3(c1,c2,d3),
    det3(c1,d2,d3)+det3(d1,c2,d3)+det3(d1,d2,c3),
    det3(d1,d2,d3)
  };
  if(std::abs(P[0]) > std::abs(P[3])) {
    libNumerics::swap(F1,F2);
    std::swap(P[0],P[3]);
    std::swap(P[1],P[2]);
  }

  // Find roots of polynomial P[3]*t^3 + P[2]*t^2 + P[1]*t + P[0]
  double roots[3];
  int num_roots = CubicRoots(P, roots);

  // Build fundamental matrix for each solution
  for (int s = 0; s < num_roots; ++s)
    Fs->push_back(F1 + roots[s]*F2);
}

/// 8-point algorithm.
/// \param A Matrix such that f is solution to Af=0
/// \param Fs Output vector of found F matrix (up to 1)
void FundamentalModel::algo8pt(const Mat& A, std::vector<Mat> *Fs) const {
  // Without rank constraint
  libNumerics::vector<double> vecNullspace(9);
  libNumerics::SVD::Nullspace(A, &vecNullspace, 1, 1);
  libNumerics::matrix<double> F(3,3);
  F.read(vecNullspace);

  // Force the rank 2 constraint
  libNumerics::matrix<double> F2(3,3);
  libNumerics::SVD::EnforceRank2_3x3(F, &F2);
  Fs->push_back(F2);
}

/// Return left epipole of matrix F. The returned vector is of norm 1.
static Vec leftEpipole(const ModelEstimator::Mat& F) {
    Vec e(3), e2(3);
    double norm=-1;
    for(int i=0; i<3; i++)
        for(int j=i+1; j<3; j++) {
            e2 = cross(F.col(i),F.col(j));
            double norm2=e2.qnorm();
            if(norm < norm2) {
                norm = norm2;
                e = e2;
            }
        }
    return (e/sqrt(norm));
}

/// Filter out F matrices that are not possible.
void FundamentalModel::checkF(const std::vector<int> &indices,
                              std::vector<Mat> *Fs) const {
  for(std::vector<Mat>::iterator F=Fs->begin(); F != Fs->end();) {
    Vec e = leftEpipole(*F);
    std::vector<int>::const_iterator it=indices.begin();
    do {
      int j=*it;
      Vec xL(x1_(0,j),x1_(1,j),1.0), xR(x2_(0,j),x2_(1,j),1.0);
      Vec exL = cross(e,xL);
      Vec FxR = *F * xR;
      double d = dot(exL,FxR);
      double t = sqrt(xL.qnorm()*xR.qnorm());
      if(d<=t*MIN_PRODUCT_NORMS) {
        if(it==indices.begin() && d<-t*MIN_PRODUCT_NORMS)
          e = -e;
        else
          break;
      }
    } while(++it != indices.end());
    if(it==indices.end())
      ++F;
    else
      F = Fs->erase(F);
  }
}

}  // namespace orsa

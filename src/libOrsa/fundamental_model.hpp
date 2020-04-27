/**
 * @file fundamental_model.cpp
 * @brief Fundamental matrix model
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

#ifndef FUNDAMENTAL_MODEL_H_
#define FUNDAMENTAL_MODEL_H_

#include "model_estimator.hpp"

namespace orsa {

///  Fundamental 7-point model, used for robust estimation.
///
/// See page 281 of book by Hartley-Zisserman.
/// The equation is \f$det(F_1 + \alpha F_2) = 0\f$.
class FundamentalModel : public ModelEstimator {
public:
  FundamentalModel(const Mat &x1, int w1, int h1,
                   const Mat &x2, int w2, int h2,
                   bool symError=false);

  /// 7 points are required to compute a fundamental matrix.
  int SizeSample() const { return  7;}

  /// Up to 3 fundamental matrices are computed from a sample of 7 points.
  int NbModels() const { return 3;}

  /// Distance used to distinguish inlier/outlier is to a line
  virtual bool DistToPoint() const { return false; }

  void Fit(const std::vector<int> &indices, std::vector<Model> *Fs) const;

  /// Square reprojection error for a given point through F.
  double Error(const Model &F, int index, int* side=0) const;
  
private:
  void Unnormalize(Model *model) const;
  void EpipolarEquation(const std::vector<int> &indices, Mat *A) const;
  void algo7pt(const Mat& A, std::vector<Mat> *Fs,
               const std::vector<int> &indices) const;
  void algo8pt(const Mat& A, std::vector<Mat> *Fs) const;
  bool checkF(const Mat &F, const std::vector<int> &indices) const;
};

}  // namespace orsa

#endif

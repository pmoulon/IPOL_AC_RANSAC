/**
 * @file model_estimator.hpp
 * @brief Model regression from sample matches.
 * @author Pascal Monasse, Pierre Moulon
 * 
 * Copyright (c) 2011 Pierre Moulon
 * Copyright (c) 2011,2020 Pascal Monasse
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

#ifndef MODEL_ESTIMATOR_H
#define MODEL_ESTIMATOR_H

#include <vector>
#include "libNumerics/matrix.h"

namespace orsa {

/// Generic class for parametric model estimation.
///
/// Subclasses must follow this interface:
///   1. SizeSample() Number correspondences necessary to compute a model.
///   2. NbModels() Number models from sample of SizeSample() correspondences.
///   3. DistToPoint() Residual is distance to a point or to a line?
///   4. Fit(const vector<size_t> &indices, vector<Model> *models)
///        Compute model(s) compatible with indexed correspondences.
///   5. Error(const Model &model, size_t index)
///        Reprojection square error for indexed correspondence.

class ModelEstimator {
 public:
  typedef libNumerics::matrix<double> Mat;
  typedef Mat Model;

  /// With true \c symError, the error is the maximum of dist(M(x1),x2) (side 1)
  /// and dist(x1,M'(x2)) (side 0). M' is the dual model of M (M'=M^-1 for
  /// homography, M'=M^T for fundamental).
  /// Return in \a side of method Error (if non-null pointer) the side this
  /// maximum is reached.
  bool symError;

  /// Constructor
  ModelEstimator(const Mat &x1, int w1, int h1,
                 const Mat &x2, int w2, int h2, bool symmetricError=false);
  virtual ~ModelEstimator() {}

  /// Number of data matches.
  int NbData() const { return x1_.ncol(); }

  /// Compute model from points.
  bool ComputeModel(const std::vector<int> &indices, Model *model) const;

  /// Minimum number of points required to compute a model.
  /// - homography -> 4
  /// - fundamental 7 pts -> 7
  /// - fundamental 8 pts -> 8
  virtual int SizeSample() const = 0;

  /// Maximum number of models possibly computed from a sample.
  /// - homography -> 1
  /// - fundamental 7 pts -> 3
  /// - fundamental 8 pts -> 1
  virtual int NbModels() const = 0;

  /// Indicate if distance used to distinguish inlier/outlier is to a point
  /// (true) or a line (false). Most regression routines, such as RANSAC, do not
  /// care about that, bur ORSA needs it.
  /// - homography -> true
  /// - fundamental -> false
  virtual bool DistToPoint() const = 0;

  /// Computes the square error of a correspondence wrt \a model.
  /// \param model The model to evaluate.
  /// \param index The point index stored in the Kernel.
  /// \param[out] side In which image is the error measured?
  virtual double Error(const Model &model, int index, int* side=0) const = 0;

  /// Computes the models associated to indexed sample.
  /// \param indices Indices of points to consider for model estimation.
  /// \param models  Estimated model(s) from sampled point.
  virtual void Fit(const std::vector<int> &indices,
                   std::vector<Model> *models) const = 0;

protected:
  Mat x1_; ///< Points in image 1
  Mat x2_; ///< Points in image 2
  Mat N1_; ///< Normalization for x1_
  Mat N2_; ///< Normalization for x2_ 
};

}  // namespace orsa

#endif

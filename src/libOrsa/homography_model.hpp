/**
 * @file homography_model.hpp
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

#ifndef HOMOGRAPHY_MODEL_H_
#define HOMOGRAPHY_MODEL_H_

#include "libOrsa/model_estimator.hpp"

namespace orsa {

/// Homography model used for robust estimation with ORSA algorithm.
class HomographyModel : public ModelEstimator {
public:
  HomographyModel(const Mat &x1, int w1, int h1,
                  const Mat &x2, int w2, int h2,
                  bool symmetricError=false);

  /// 4 point correspondences required to compute a homography.
  int SizeSample() const { return  4; }

  /// Only 1 homography can be estimated from a sample of 4 points.
  int NbModels() const { return 1; }

  /// Distance used to distinguish inlier/outlier is to a point
  virtual bool DistToPoint() const { return true; }

  /// Estimated homography satisfies the equation y = H x.
  void Fit(const std::vector<int> &indices, std::vector<Mat> *H) const;

  /// Square reprojection error for a given point through the model H.
  double Error(const Model &H, int index, int* side=0) const;

private:
  void Unnormalize(Model *model) const;
  bool IsOrientationPreserving(const std::vector<int> &indices,
                                 const Mat& H) const;
};

}  // namespace orsa

#endif

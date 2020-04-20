/**
 * @file ransac.hpp
 * @brief Model estimation by classical RANSAC algorithm.
 * @author Pascal Monasse
 * 
 * Copyright (c) 2020 Pascal Monasse
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

#ifndef RANSAC_H
#define RANSAC_H

#include "model_estimator.hpp"

namespace orsa {

/// Model estimation with ORSA algorithm.
class Ransac {
public:
  /// Constructor
  Ransac(const ModelEstimator* estimator);

  /// Generic implementation of RANSAC
  size_t run(std::vector<int> &vec_inliers,
             double precision,
             size_t nIterMax = 1000,
             ModelEstimator::Model *model = NULL,
             double beta = 0.95,
             bool bVerbose = false) const;

private:
  const ModelEstimator* estimator_;
};

}  // namespace orsa

#endif

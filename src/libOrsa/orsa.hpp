/**
 * @file orsa.hpp
 * @brief Model estimation by ORSA (aka AC-RANSAC) algorithm.
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

#ifndef ORSA_H
#define ORSA_H

#include <vector>
#include "model_estimator.hpp"

namespace orsa {

/// Model estimation with ORSA algorithm.
class Orsa {
 public:
  typedef libNumerics::vector<double> Vec;
  typedef libNumerics::matrix<double> Mat;
  typedef Mat Model;

  /// Constructor
  Orsa(const ModelEstimator* estimator, double alpha0Left, double alpha0Right);

  /// Generic implementation of ORSA (Optimized Random Sampling Algorithm). 
  double run(std::vector<int> &vec_inliers,
             size_t nIter = 1000,
             double *precision = NULL,
             ModelEstimator::Model *model = NULL,
             bool bVerbose = false) const;

  /// Enable the model refinement until convergence.
  void setRefineUntilConvergence(bool value);

  /// Return if convergence check is on or off.
  bool getRefineUntilConvergence() const;

private:
  const ModelEstimator* estimator_;
  double logalpha0_[2]; ///< Log probability of error<=1, set by subclass
  bool bConvergence;

private:
  /// Distance and associated index
  struct ErrorIndex {
    double error; ///< Square error
    int index; ///< Correspondence index
    int side;     ///< Error in image 1 (side=0) or 2 (side=1)?
    /// Constructor
    ErrorIndex(double e=0, int i=0, int s=0): error(e), index(i), side(s) {}
    bool operator<(const ErrorIndex& e) const { return (error<e.error); }
  };
  ErrorIndex bestNFA(const std::vector<ErrorIndex>& e,
                     double loge0, double maxThreshold,
                     const std::vector<float> & vec_logc_n,
                     const std::vector<float> & vec_logc_k) const;

  /// Iterative minimization NFA/RMSE.
  void refineUntilConvergence(const std::vector<float> & vec_logc_n,
                              const std::vector<float> & vec_logc_k,
                              double loge0,
                              double maxThreshold,
                              double minNFA,
                              ModelEstimator::Model *model,
                              bool bVerbose,
                              std::vector<int> & vec_inliers,
                              double & errorMax,
                              int & side) const;
};

}  // namespace orsa

#endif

/**
 * @file ransac.cpp
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

#include "ransac.hpp"
#include "sampling.hpp"
#include <algorithm>
#include <iostream>
#include <limits>
#include <cmath>

namespace orsa {

/// Constructor
Ransac::Ransac(const ModelEstimator* estimator)
: estimator_(estimator) {}

/// Find inliers of \a model passed as parameter (within \a precision).
void Ransac::find_inliers(const ModelEstimator::Model& model,
                          double precision, std::vector<int>& inliers) const {
  const int nData = estimator_->NbData();
  for(int i=0; i<nData; i++)
    if(estimator_->Error(model, i)<=precision)
      inliers.push_back(i);
}

/// Some information output
void display_info(size_t nInliers, size_t iter, size_t nIterMax,
                  const std::vector<int>& vec_sample){
  std::cout << " inliers=" << nInliers
            << " (iter=" << iter
            << ",iterMax=" << nIterMax;
  std::cout << ",sample=" << vec_sample.front();
  std::vector<int>::const_iterator it=vec_sample.begin();
  for(++it; it != vec_sample.end(); ++it)
      std::cout << ',' << *it;
  std::cout << ")" <<std::endl;
}

/// Generic implementation of RANSAC
size_t Ransac::run(std::vector<int> &vec_inliers,
                   double precision,
                   size_t nIterMax,
                   ModelEstimator::Model *model,
                   double beta,
                   bool bVerbose) const {
  assert(beta>0);
  if(beta > 1.0) {
    std::cerr << "RANSAC beta parameter adjusted to not exceed 1" <<std::endl;
    beta = 1.0;
  }
  double log1_beta = log(1-beta);
  precision *= precision;

  const int nData = estimator_->NbData();
  const int sizeSample = estimator_->SizeSample();

  vec_inliers.clear();
  size_t iter=0;
  for(; iter<nIterMax; iter++) {
    std::vector<int> vec_sample(sizeSample); // Sample indices
    UniformSample(sizeSample, nData, &vec_sample); // Get random sample
    std::vector<ModelEstimator::Model> vec_models;
    estimator_->Fit(vec_sample, &vec_models);
    for (size_t k = 0; k < vec_models.size(); ++k) {
      std::vector<int> inliers;
      find_inliers(vec_models[k], precision, inliers);
      if(vec_inliers.size() < inliers.size()) {
        if(model) *model = vec_models[k];
        std::swap(inliers,vec_inliers); // Avoid copy
        float denom = log(1-pow(vec_inliers.size()/(double)nData, sizeSample));
        if(denom<0) // Protect against 1-eps==1
          nIterMax = (size_t)std::min((double)nIterMax,ceil(log1_beta/denom));
        if(bVerbose)
          display_info(vec_inliers.size(), iter, nIterMax, vec_sample);
      }
    }
  }
  return iter;
}

}  // namespace orsa

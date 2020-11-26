/**
 * @file eval_model.cpp
 * @brief Model estimation with ORSA or RANSAC algorithm
 * @author Lionel Moisan, Pascal Monasse, Pierre Moulon
 * 
 * Copyright (c) 2011-2016 Lionel Moisan
 * Copyright (c) 2011-2016,2020 Pascal Monasse
 * Copyright (c) 2011-2016 Pierre Moulon
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

#include "eval_model.hpp"
#include "homography_model.hpp"
#include "fundamental_model.hpp"
#include "orsa.hpp"
#include "ransac.hpp"
#include <iostream>

namespace orsa {

/// Display RMSE/max error of inliers of model \a M.
static std::pair<double,double> display_stats(const ModelEstimator& model,
                                              const std::vector<int>& in,
                                              const ModelEstimator::Model& M) {
  std::vector<int>::const_iterator it=in.begin();
  double l2=0, linf=0;
  for(; it!=in.end(); ++it) {
    double e = model.Error(M, *it);
    l2 += e;
    if(linf < e)
      linf = e;
  }
  std::pair<double,double> err(sqrt(l2/in.size()),sqrt(linf));
  std::cout << "RMSE/max error: " << err.first << "/" << err.second <<std::endl;
  return err;
}

/// Generate a model estimator. Do not forget deletion.
template <typename Model>
ModelEstimator* build_model(const std::vector<Match>& vec_matchings) {
  return new Model(vec_matchings, true);
}

/// Refine model based on all inliers, and display statistics.
static void refine(ModelEstimator* model,
                   const std::vector<int>& vec_inliers,
                   ModelEstimator::Model *M) {
  std::cout << "Before refinement: ";
  std::pair<double,double> err = display_stats(*model, vec_inliers, *M);
  ModelEstimator::Model M2(3,3);
  if( model->ComputeModel(vec_inliers,&M2) ) // Re-estimate with all inliers
  {
    std::cout << "After  refinement: ";
    double maxBefore = err.second;
    if(display_stats(*model, vec_inliers, M2).first <= maxBefore)
      *M = M2;
    else
      std::cerr << "Warning: error after refinement is too large, thus ignored"
                <<std::endl;
  } else
    std::cerr << "Warning: error in refinement, result is suspect" <<std::endl;
}

/// Estimate the model using regular RANSAC and refinement.
/// \param[in] model Model estimator.
/// \param[in] precision Maximum inlier/outlier threshold (in pixels).
/// \param[in] nbIterMax Maximal number of iterations for RANSAC algorithm.
/// \param[in] beta Probability of one correct sample (to adjust iterations).
/// \param[out] M Model registering left on right image.
/// \param[out] vec_inliers Index of inliers in \a vec_matchings.
/// \return true if at least one viable sample (ie, producing a model) is found.
bool generic_ransac(ModelEstimator* model,
                    double precision, int nbIterMax, double beta,
                    libNumerics::matrix<double>& M,
                    std::vector<int>& vec_inliers) {
  if(model->NbData() >= model->SizeSample()) {
    Ransac ransac(model);
    nbIterMax = ransac.run(vec_inliers, precision, nbIterMax, &M, beta, true);
    std::cout << "Iterations: " << nbIterMax << std::endl;
    refine(model, vec_inliers, &M);
    return true;
  }
  std::cerr << "Error: RANSAC needs " << model->SizeSample()
            << " matches or more to proceed" <<std::endl;
  return false;
}

/// Estimate the model using ORSA method and refinement.
/// \param[in] model Model estimator.
/// \param[in] alpha0Left, alpha0Right Probabilities of error of 1 pixel.
/// \param[in,out] precision Maximum inlier/outlier threshold (in pixels).
/// \param[in] nbIter Maximal number of iterations for RANSAC algorithm.
/// \param[out] M model registering left on right image.
/// \param[out] vec_inliers Index of inliers in \a vec_matchings.
/// \return The status of the estimation. If no meaningful (NFA<1) model is
/// found, the \a M field is not filled and \a vec_inliers is empty.
bool generic_orsa(ModelEstimator* model, double alpha0Left, double alpha0Right,
                  double precision, int nbIter,
                  libNumerics::matrix<double>& M,
                  std::vector<int>& vec_inliers) {
  if(model->NbData() > model->SizeSample()) {
    Orsa orsa(model, alpha0Left, alpha0Right);
    if(orsa.run(vec_inliers, nbIter, &precision, &M, true)>0.0)
      return false;
    refine(model, vec_inliers, &M);
    return true;
  }
  std::cerr << "Error: ORSA needs " << model->SizeSample()+1
            << " matches or more to proceed" <<std::endl;
  return false;
}

/// Estimate the homography using regular RANSAC and refinement.
/// \param[in] vec_matchings List of correspondences.
/// \param[in] precision Maximum inlier/outlier threshold (in pixels).
/// \param[in] nbIterMax Maximal number of iterations for RANSAC algorithm.
/// \param[in] beta Probability of one correct sample (to adjust iterations).
/// \param[out] H Homography registering left on right image.
/// \param[out] vec_inliers Index of inliers in \a vec_matchings.
/// \return true if at least one viable sample (ie, producing a model) is found.
bool ransac_homography(const std::vector<Match>& vec_matchings,
                       double precision, int nbIterMax, double beta,
                       libNumerics::matrix<double>& H,
                       std::vector<int>& vec_inliers) {
  if(vec_matchings.empty()) return false;
  ModelEstimator* model = build_model<HomographyModel>(vec_matchings);
  bool ok = generic_ransac(model, precision, nbIterMax, beta, H, vec_inliers);
  delete model;
  return ok;
}

/// Estimate the fundamental matrix using regular RANSAC and refinement.
/// \param[in] vec_matchings List of correspondences.
/// \param[in] precision Maximum inlier/outlier threshold (in pixels).
/// \param[in] nbIterMax Maximal number of iterations for RANSAC algorithm.
/// \param[in] beta Probability of one correct sample (to adjust iterations).
/// \param[out] F Fundamental matrix between left and right image.
/// \param[out] vec_inliers Index of inliers in \a vec_matchings.
/// \return true if at least one viable sample (ie, producing a model) is found.
bool ransac_fundamental(const std::vector<Match>& vec_matchings,
                        double precision, int nbIterMax, double beta,
                        libNumerics::matrix<double>& F,
                        std::vector<int>& vec_inliers) {
  if(vec_matchings.empty()) return false;
  ModelEstimator* model = build_model<FundamentalModel>(vec_matchings);
  bool ok = generic_ransac(model, precision, nbIterMax, beta, F, vec_inliers);
  delete model;
  return ok;
}

/// Estimate the homography using ORSA method and refinement.
/// \param[in] vec_matchings List of correspondences.
/// \param[in] w1,h1 Dimensions of left image.
/// \param[in] w2,h2 Dimensions of right image.
/// \param[in,out] precision Maximum inlier/outlier threshold (in pixels).
/// \param[in] nbIter Maximal number of iterations for RANSAC algorithm.
/// \param[out] H Homography registering left on right image.
/// \param[out] vec_inliers Index of inliers in \a vec_matchings.
/// \return The status of the estimation. If no meaningful (NFA<1) model is
/// found, the \a H field is not filled and \a vec_inliers is empty.
bool orsa_homography(const std::vector<Match>& vec_matchings,
                     int w1,int h1, int w2,int h2,
                     double precision, int nbIter,
                     libNumerics::matrix<double>& F,
                     std::vector<int>& vec_inliers) {
  if(vec_matchings.empty()) return false;
  ModelEstimator* model = build_model<HomographyModel>(vec_matchings);
  double alpha0Left  = M_PI/(w1*(double)h1);
  double alpha0Right = M_PI/(w2*(double)h2);
  bool ok = generic_orsa(model, alpha0Left, alpha0Right, precision, nbIter,
                         F, vec_inliers);
  delete model;
  return ok;
}

/// Estimate the fundamental matrix using ORSA method and refinement.
/// If the mean error after refinement exceeds the max error of ORSA result,
/// the refinement is not applied.
/// \param[in] vec_matchings List of correspondences.
/// \param[in] w1,h1 Dimensions of left image.
/// \param[in] w2,h2 Dimensions of right image.
/// \param[in,out] precision Maximum inlier/outlier threshold (in pixels).
/// \param[in] nbIter Maximal number of iterations for RANSAC algorithm.
/// \param[out] F Fundamental matrix between left and right image.
/// \param[out] vec_inliers Index of inliers in \a vec_matchings.
/// \return The status of the estimation. If no meaningful (NFA<1) model is
/// found, the \a F field is not filled and \a vec_inliers is empty.
bool orsa_fundamental(const std::vector<Match>& vec_matchings,
                      int w1,int h1, int w2,int h2,
                      double precision, int nbIter,
                      libNumerics::matrix<double>& F,
                      std::vector<int>& vec_inliers) {
  if(vec_matchings.empty()) return false;
  ModelEstimator* model = build_model<FundamentalModel>(vec_matchings);
  double D, A; // Diameter and area of image
  D = sqrt(w1*(double)w1 + h1*(double)h1);
  A = w1*(double)h1;
  double alpha0Left  = 2.0*D/A;
  D = sqrt(w2*(double)w2 + h2*(double)h2);
  A = w2*(double)h2;
  double alpha0Right = 2.0*D/A;
  bool ok = generic_orsa(model, alpha0Left, alpha0Right, precision, nbIter,
                         F, vec_inliers);
  delete model;
  return ok;
}

} // namespace orsa

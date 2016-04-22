/**
 * @file orsa_homography.cpp
 * @brief Homography estimation with ORSA algorithm
 * @author Lionel Moisan, Pascal Monasse, Pierre Moulon
 * 
 * Copyright (c) 2011-2016 Lionel Moisan, Pascal Monasse, Pierre Moulon
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

#include "orsa_homography.hpp"
#include "homography_model.hpp"
#include "extras/libNumerics/homography.h"

#include <iostream>

namespace orsa {

/// Display average/max error of inliers of homography H.
static void display_stats(const std::vector<Match>& vec_matchings,
                          const std::vector<int>& vec_inliers,
                          libNumerics::matrix<double>& H) {
  std::vector<int>::const_iterator it=vec_inliers.begin();
  double l2=0, linf=0;
  for(; it!=vec_inliers.end(); ++it) {
    const Match& m=vec_matchings[*it];
    double x1=m.x1, y1=m.y1;
    TransformH(H, x1, y1);
    double e = (m.x2-x1)*(m.x2-x1) + (m.y2-y1)*(m.y2-y1);
    l2 += e;
    if(linf < e)
      linf = e;
  }
  std::cout << "Average/max error: "
            << sqrt(l2/vec_inliers.size()) << "/"
            << sqrt(linf) <<std::endl;
}

/// Estimate the homography using ORSA method and refinement.
/// \param[in] vec_matchings List of correspondences.
/// \param[in] w1,h1 Dimensions of left image.
/// \param[in] w2,h2 Dimensions of right image.
/// \param[in] precision Maximum inlier/outlier threshold (in pixels).
/// \param[in] nbIter Maximal number of iterations for RANSAC algorithm.
/// \param[out] H Homography registering left on right image.
/// \param[out] vec_inliers Index of inliers in \a vec_matchings.
/// \return The status of the estimation. If no meaningful (NFA<1) model is
/// found, the \a H field is not filled and \a vec_inliers is empty.
bool orsa_homography(const std::vector<Match>& vec_matchings,
                     int w1,int h1, int w2,int h2,
                     double precision, int nbIter,
                     libNumerics::matrix<double>& H,
                     std::vector<int>& vec_inliers)
{
  const int n = static_cast<int>( vec_matchings.size() );
  if(n < 5)
  {
      std::cerr << "Error: ORSA needs 5 matches or more to proceed" <<std::endl;
      return false;
  }
  libNumerics::matrix<double> xA(2,n), xB(2,n);

  for (int i=0; i < n; ++i)
  {
    xA(0,i) = vec_matchings[i].x1;
    xA(1,i) = vec_matchings[i].y1;
    xB(0,i) = vec_matchings[i].x2;
    xB(1,i) = vec_matchings[i].y2;
  }

  orsa::HomographyModel model(xA, w1, h1, xB, w2, h2, true);
  //model.setConvergenceCheck(true);

  if(model.orsa(vec_inliers, nbIter, &precision, &H, true)>0.0)
    return false;
  std::cout << "Before refinement: ";
  display_stats(vec_matchings, vec_inliers, H);
  if( model.ComputeModel(vec_inliers,&H) ) // Re-estimate with all inliers
  {
    std::cout << "After  refinement: ";
    display_stats(vec_matchings, vec_inliers, H);
  } else
    std::cerr << "Warning: error in refinement, result is suspect" <<std::endl;
  return true;
}

} // namespace orsa

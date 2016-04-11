/**
 * @file orsa_fundamental.cpp
 * @brief Fundamental matrix estimation
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

#include <cstdlib>
#include <ctime>

#include <iostream>
#include <vector>
#include <utility>

#include "libImage/image.hpp"
#include "libImage/image_io.hpp"

#include "extras/libNumerics/numerics.h"
#include "extras/sift/library.h"

#include "libOrsa/fundamental_model.hpp"

#include "demo/siftMatch.hpp"
#include "demo/cmdLine.h"
#include "demo/fundamental_graphical_output.hpp"

/// Number of random samples in ORSA
static const int ITER_ORSA=10000;

/// Display average/max error of inliers of fundamental matrix F.
static std::pair<double,double> display_stats(const std::vector<Match>& vec_matchings,
                                              const std::vector<int>& vec_inliers,
                                              const libNumerics::matrix<double>& F) {
    std::vector<int>::const_iterator it=vec_inliers.begin();
    double l2=0, linf=0;
    for(; it!=vec_inliers.end(); ++it) {
      const Match& m=vec_matchings[*it];
      double a = F(0,0) * m.x1 + F(1,0) * m.y1 + F(2,0);
      double b = F(0,1) * m.x1 + F(1,1) * m.y1 + F(2,1);
      double c = F(0,2) * m.x1 + F(1,2) * m.y1 + F(2,2);
      double d = a*m.x2 + b*m.y2 + c;
      double e =  (d*d) / (a*a + b*b);
      l2 += e;
      if(linf < e)
        linf = e;
    }
    std::pair<double,double> err(sqrt(l2/vec_inliers.size()),sqrt(linf));
    std::cout << "Average/max error: " << err.first << "/" << err.second <<std::endl;
    return err;
}

/// ORSA fundamental estimation
bool ORSA(const std::vector<Match>& vec_matchings, int w1,int h1, int w2,int h2,
          double precision,
          libNumerics::matrix<double>& F, std::vector<int>& vec_inliers)
{
  const int n = static_cast<int>( vec_matchings.size() );
  if(n < 7)
  {
    std::cerr << "Error: ORSA needs 7 matches or more to proceed" <<std::endl;
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

  orsa::FundamentalModel model(xA, w1, h1, xB, w2, h2);
  //model.setConvergenceCheck(true);

  if(model.orsa(vec_inliers, ITER_ORSA, &precision, &F, true)>0.0)
    return false;
  std::pair<double,double> err; // (RMSE,max)
  std::cout << "Before refinement: ";
  err = display_stats(vec_matchings, vec_inliers, F);
  libNumerics::matrix<double> F2(3,3);
  if( model.ComputeModel(vec_inliers,&F2) ) // Re-estimate with all inliers
  {
    std::cout << "After  refinement: ";
    double maxBefore = err.second;
    if(display_stats(vec_matchings, vec_inliers, F2).first <= maxBefore)
      F = F2;
    else
      std::cerr << "Warning: error after refinement is too large, thus ignored" <<std::endl;
  } else
    std::cerr << "Warning: error in refinement, result is suspect" <<std::endl;
  return true;
}

int main(int argc, char **argv)
{
  double precision=0;
  float fSiftRatio=0.6f;
  CmdLine cmd;
  cmd.add( make_option('p',precision, "prec") );
  cmd.add( make_option('s',fSiftRatio, "sift") );
  cmd.add( make_switch('r', "read") );
  try {
    cmd.process(argc, argv);
  } catch(const std::string& s) {
    std::cerr << s << std::endl;
    return 1;
  }
  if(argc!=5 && argc!=6 && argc!=7 && argc!=8) {
    std::cerr << "Usage: " << argv[0] << " imgInA imgInB "
              << "[-p|--prec precision] "
              << "[-s|--sift siftRatio] "
              << "[-r|--read] "
              << "allMatches.txt orsaMatches.txt "
              << "[imgInliers imgOutliers] [imgEpipolar]"
              << std::endl;
    return 1;
  }

  // Init random seed
  srand((unsigned int)time(0));

  Image<RGBColor> image1, image2;
  if(! libs::ReadImage(argv[1], &image1))
    return 1;
  if(! libs::ReadImage(argv[2], &image2))
    return 1;

  Image<unsigned char> image1Gray, image2Gray;
  libs::convertImage(image1, &image1Gray);
  libs::convertImage(image2, &image2Gray);

  // Find matches with SIFT or read correspondence file
  std::vector<Match> vec_matchings;
  if(cmd.used('r')) {
      if(Match::loadMatch(argv[3], vec_matchings))
          std::cout << "Read " <<vec_matchings.size()<< " matches" <<std::endl;
      else {
          std::cerr << "Failed reading matches from " << argv[3] <<std::endl;
          return 1;
      }
  } else
      SIFT(image1Gray, image2Gray, vec_matchings, fSiftRatio);

  // Remove duplicates (frequent with SIFT)
  rm_duplicates(vec_matchings);

  // Save match files
  if(! cmd.used('r') && ! Match::saveMatch(argv[3], vec_matchings)) {
    std::cerr << "Failed saving matches into " <<argv[3] <<std::endl;
    return 1;
  }

  int w1 = image1Gray.Width(), h1 = image1Gray.Height();
  int w2 = image2Gray.Width(), h2 = image2Gray.Height();

  // Estimation of fundamental matrix with ORSA
  libNumerics::matrix<double> F(3,3);
  std::vector<int> vec_inliers;
  bool ok = ORSA(vec_matchings, w1, h1, w2, h2, precision, F, vec_inliers);
  if(ok)
    std::cout << "F=" << F <<std::endl;

  std::vector<Match> good_match;
  std::vector<int>::const_iterator it = vec_inliers.begin();
  for(; it != vec_inliers.end(); it++)
    good_match.push_back(vec_matchings[*it]);
  if(! Match::saveMatch(argv[4], good_match)) {
    std::cerr << "Failed saving matches into " <<argv[4] <<std::endl;
    return 1;
  }

  if(argc>6) // Output images
  {
    const char *fileIn=0, *fileOut=0, *fileEpi=0;
    fileIn  = argv[5];
    fileOut = argv[6];
    if(argc!=7)
      fileEpi = argv[argc-1];
    fundamental_graphical_output(image1Gray, image2Gray,
                                 vec_matchings, vec_inliers, ok?&F:0,
                                 fileIn, fileOut, fileEpi);
  }

  if(! ok)
  {
    std::cerr << "Failed to estimate a model" << std::endl;
    return 1;
  }

  return 0;
}

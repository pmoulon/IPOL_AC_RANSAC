/**
 * @file show_epi.cpp
 * @brief Show epipolar matches between two images, useful for USAC
 * @author Pascal Monasse
 * 
 * Copyright (c) 2017 Pascal Monasse
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

#include <iostream>
#include <fstream>
#include <numeric>

#include "libImage/image_io.hpp"
#include "libOrsa/orsa_fundamental.hpp"

#include "siftMatch.hpp"
#include "fundamental_graphical_output.hpp"

using namespace libNumerics;

// Read calibration file in USAC format
bool read_calib(const char* sFile, matrix<double>& K1, matrix<double>& K2) {
    std::ifstream file(sFile);
    if(! file)
        return false;
    K1.fill(0); K1(2,2) = 1;
    K2.fill(0); K2(2,2) = 1;
    file >> K1(0,0) >> K1(0,1) >> K1(0,2)
                    >> K1(1,1) >> K1(1,2);
    file >> K2(0,0) >> K2(0,1) >> K2(0,2)
                    >> K2(1,1) >> K2(1,2);
        
    return !file.fail();
}

bool read_E(const char* sFile, matrix<double>& E) {
    std::ifstream file(sFile);
    if(! file)
        return false;
    file >> E(0,0) >> E(0,1) >> E(0,2)
         >> E(1,0) >> E(1,1) >> E(1,2)
         >> E(2,0) >> E(2,1) >> E(2,2);
    return !file.fail();
}

/// Return left epipole of matrix F. The returned vector is of norm 1.
static vector<double> leftEpipole(const matrix<double>& F) {
    vector<double> e(3), e2(3);
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
    return e/e(2);
}

int main(int argc, char **argv)
{
  if(argc != 4 && argc!=7) {
    std::cerr << "Usage: " << argv[0]
              << " [imgInA imgInB calib_matrices.txt] "
              << "E.txt "
              << "matches.txt "
              << "imgOut"
              << std::endl;
    return 1;
  }
  const char* sImage1 = "im1.jpg";
  const char* sImage2 = "im2.jpg";
  const char* sFileCalib = "calib_matrices.txt";
  char* sFileE = argv[argc-3];
  char* sFileMatch = argv[argc-2];
  char* sFileOut = argv[argc-1];
  if(argc > 4) {
      sImage1 = argv[1];
      sImage2 = argv[2];
      sFileCalib = argv[3];
  }

  // Read images
  Image<RGBColor> image1, image2;
  if(! (libs::ReadImage(sImage1, &image1) && libs::ReadImage(sImage2, &image2)))
    return 1;
  Image<unsigned char> image1Gray, image2Gray;
  libs::convertImage(image1, &image1Gray);
  libs::convertImage(image2, &image2Gray);

  // Read correspondence file
  std::vector<Match> matchings;
  if(Match::loadMatch(sFileMatch, matchings))
    std::cout << "Read " <<matchings.size()<< " matches" <<std::endl;
  else {
    std::cerr << "Failed reading matches from " << argv[3] <<std::endl;
    return 1;
  }

  // Read calibration file
  matrix<double> K1(3,3), K2(3,3);
  if(! read_calib(sFileCalib, K1, K2)) {
    std::cerr << "Error reading file " << sFileCalib << std::endl;
    return 1;
  }

  matrix<double> F(3,3);
  if(! read_E(sFileE, F)) {
    std::cerr << "Error reading file " << sFileE << std::endl;
    return 1;
  }
  F = K1.t().inv() * F * K2.inv();
  std::cout << "F=" << F <<std::endl;
  std::cout << "Epipoles: "
            << leftEpipole(F    ).copy(0,1) << " "
            << leftEpipole(F.t()).copy(0,1) << std::endl;

  std::vector<int> vec_inliers(matchings.size(), 1);
  if(!vec_inliers.empty()) vec_inliers[0]=0;
  std::partial_sum(vec_inliers.begin(), vec_inliers.end(), vec_inliers.begin());

  fundamental_graphical_output(image1Gray, image2Gray,
                               matchings, vec_inliers, &F,
                               0, 0, sFileOut);
  return 0;
}

/**
 * @file demo_orsa_fundamental.cpp
 * @brief Fundamental matrix estimation with ORSA algorithm
 * @author Lionel Moisan, Pascal Monasse, Pierre Moulon
 * 
 * Copyright (c) 2011-2018 Lionel Moisan, Pascal Monasse, Pierre Moulon
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

#include "libImage/image_io.hpp"
#include "libOrsa/orsa_fundamental.hpp"

#include "cmdLine.h"
#include "siftMatch.hpp"
#include "fundamental_graphical_output.hpp"

/// Number of random samples in ORSA
static const int ITER_ORSA=10000;

int main(int argc, char **argv)
{
  double precision=0;
  float fSiftRatio=0.6f;
  CmdLine cmd;
  cmd.add( make_option('p',precision, "prec")
           .doc("max precision (in pixels) of registration (0=arbitrary)") );
  cmd.add( make_option('s',fSiftRatio, "sift")
           .doc("SIFT distance ratio of descriptors") );
  cmd.add( make_switch('r', "read")
           .doc("Read file of matches allMatches.txt, do not use SIFT"));
  try {
    cmd.process(argc, argv);
  } catch(const std::string& s) {
    std::cerr << s << std::endl;
    return 1;
  }
  if(argc!=5 && argc!=6 && argc!=7 && argc!=8) {
    std::cerr << "Usage: " << argv[0] << " [options] imgInA imgInB "
              << "allMatches.txt orsaMatches.txt "
              << "[imgInliers imgOutliers] [imgEpipolar]\n"
              << "- imgInA, imgInB: the two input image (JPG or PNG format)\n"
              << "- allMatches.txt: output (input if option -r) text file of "
                 "format \"x1 y1 x2 y2\"\n"
              << "- orsaMatches.txt: output, but only with inliers.\n"
              << "- imgInliers (optional): output image showing inliers\n"
              << "- imgOutliers (optional): output image showing outliers and "
                 "their error\n"
              << "- imgEpipolar (optional): output image with small epipolar "
                 " lines\n"
              << "\tOptions:\n" << cmd;
      return 1;
  }

  // Init random seed
  srand((unsigned int)time(0));

  // Read images
  Image<RGBColor> image1, image2;
  if(! (libs::ReadImage(argv[1], &image1) && libs::ReadImage(argv[2], &image2)))
    return 1;
  Image<unsigned char> image1Gray, image2Gray;
  libs::convertImage(image1, &image1Gray);
  libs::convertImage(image2, &image2Gray);

  // Find matches with SIFT or read correspondence file
  std::vector<Match> matchings;
  if(cmd.used('r')) {
      if(Match::loadMatch(argv[3], matchings))
          std::cout << "Read " <<matchings.size()<< " matches" <<std::endl;
      else {
          std::cerr << "Failed reading matches from " << argv[3] <<std::endl;
          return 1;
      }
  } else
      SIFT(image1Gray, image2Gray, matchings, fSiftRatio);
  rm_duplicates(matchings); // Remove duplicates (frequent with SIFT)

  // Save match files
  if(! cmd.used('r') && ! Match::saveMatch(argv[3], matchings)) {
    std::cerr << "Failed saving matches into " <<argv[3] <<std::endl;
    return 1;
  }

  // Estimation of fundamental matrix with ORSA
  libNumerics::matrix<double> F(3,3);
  std::vector<int> vec_inliers;
  int w1 = image1Gray.Width(), h1 = image1Gray.Height();
  int w2 = image2Gray.Width(), h2 = image2Gray.Height();
  bool ok = orsa::orsa_fundamental(matchings, w1,h1,w2,h2, precision, ITER_ORSA,
                                   F, vec_inliers);
  if(ok)
    std::cout << "F=" << F <<std::endl;

  // Save inliers
  std::vector<Match> good_match;
  std::vector<int>::const_iterator it = vec_inliers.begin();
  for(; it != vec_inliers.end(); it++)
    good_match.push_back(matchings[*it]);
  if(! Match::saveMatch(argv[4], good_match))
  {
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
                                 matchings, vec_inliers, ok?&F:0,
                                 fileIn, fileOut, fileEpi);
  }

  if(! ok)
  {
    std::cerr << "Failed to estimate a model" << std::endl;
    return 1;
  }

  return 0;
}

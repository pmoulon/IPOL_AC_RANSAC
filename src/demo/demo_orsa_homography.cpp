/**
 * @file demo_orsa_homography.cpp
 * @brief Homographic image registration
 * @author Pascal Monasse, Pierre Moulon
 *
 * Copyright (c) 2011-2018 Lionel Moisan
 * Copyright (c) 2011-2018,2020 Pascal Monasse
 * Copyright (c) 2011-2018 Pierre Moulon
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
#include "libImage/image_crop.hpp"

#include "libOrsa/eval_model.hpp"
#include "libOrsa/libNumerics/homography.h"

#include "cmdLine.h"
#include "siftMatch.hpp"
#include "warping.hpp"
#include "homography_graphical_output.hpp"

/// Number of random samples in ORSA
static const int ITER_ORSA=10000;

int main(int argc, char **argv)
{
  double precision=0;
  float fSiftRatio=0.6f;
  double beta=0.95;
  CmdLine cmd;
  Geometry region; region.x0=region.y0=region.x1=region.y1=0;
  cmd.add( make_option('c',region,"cut")
           .doc("cut region of imagInA: wxh+x+y = rect [x,x+w] x [y,y+h]") );
  cmd.add( make_option('p',precision, "prec")
           .doc("max precision (in pixels) of registration (0=arbitrary)") );
  cmd.add( make_switch('r', "read")
           .doc("Read file of matches allMatches.txt, do not use SIFT") );
  cmd.add( make_option('s',fSiftRatio, "sift")
           .doc("SIFT distance ratio of descriptors") );
  cmd.add( make_option('b',beta,"beta")
           .doc("Beta iteration adjustment parameter (use RANSAC)") );
  try {
    cmd.process(argc, argv);
  } catch(const std::string& s) {
    std::cerr << s << std::endl;
    return 1;
  }
  if(argc!=5 && argc!=7 && argc!=8 && argc!=10) {
    std::cerr << "Usage: " << argv[0] << " [options] imgInA imgInB "
              << "allMatches.txt orsaMatches.txt "
              << "[imgInliers imgOutliers [imgMosaic "
              << "[imgMosaicA imgMosaicB]]\n"
              << "- imgInA, imgInB: the two input image (JPG or PNG format)\n"
              << "- allMatches.txt: output (input if option -r) text file of "
                 "format \"x1 y1 x2 y2\"\n"
              << "- orsaMatches.txt: output, but only with inliers.\n"
              << "- imgInliers (optional): output image showing inliers\n"
              << "- imgOutliers (optional): output image showing outliers and "
                 "their error\n"
              << "- imgMosaic (optional): output mosaic image\n"
              << "- imgMosaicA,imgMosaicB (optional): registered images\n"
              << "\tOptions:\n" << cmd;
    return 1;
  }

  // Init random seed
  time_t seed = time(0); // Replace by a fixed value to debug a reproducible run
  srand((unsigned int)seed);

  Image<RGBColor> image1, image2;
  if(! libs::ReadImage(argv[1], &image1))
    return 1;
  if(! libs::ReadImage(argv[2], &image2))
    return 1;

  bool bUseRegion = cmd.used('c');
  if(bUseRegion) { // Sanity check
    Geometry zone;
    zone.x0=zone.y0=0; zone.x1=int(image1.Width());zone.y1=int(image1.Height());
    if(! (region.x0<region.x1 && region.y0<region.y1 &&
          zone.inside(region.x0,region.y0) &&
          zone.inside(region.x1-1,region.y1-1))) {
      std::cout << "Invalid cut region " << region << std::endl;
      return 1;
    }
  }

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
      SIFT(image1Gray, image2Gray, vec_matchings, fSiftRatio,
           bUseRegion? &region: 0);

  if(bUseRegion) {
      Image<RGBColor> image = image1;
      Crop(image, region.x0, region.y0,
           region.x1-region.x0, region.y1-region.y0, image1);
      libs::convertImage(image1, &image1Gray);
  }

  // Remove duplicates (frequent with SIFT)
  rm_duplicates(vec_matchings);

  // Save match files
  if(! cmd.used('r') && ! Match::saveMatch(argv[3], vec_matchings)) {
    std::cerr << "Failed saving matches into " <<argv[3] <<std::endl;
    return 1;
  }

  bool bUseRansac = cmd.used('b');
  if(bUseRansac && !cmd.used('p')) {
    precision = 1;
    std::cerr << "No input for RANSAC threshold (option -p)."
              << " Using " << precision << " pixel" << std::endl;
  }

  const int w1=int(image1Gray.Width()), h1=int(image1Gray.Height());
  const int w2=int(image2Gray.Width()), h2=int(image2Gray.Height());

  // Estimation of homography with ORSA
  libNumerics::matrix<double> H(3,3);
  std::vector<int> vec_inliers;
  bool ok = false;
  if(bUseRansac)
    ok = orsa::ransac_homography(vec_matchings,w1,h1,w2,h2,precision,ITER_ORSA,
                                 beta, H, vec_inliers);
  else
    ok = orsa::  orsa_homography(vec_matchings,w1,h1,w2,h2,precision,ITER_ORSA,
                                 H, vec_inliers);
  if(ok) {
    H /= H(2,2);
    std::cout << "H=" << H <<std::endl;
  }

  std::vector<Match> good_match;
  std::vector<int>::const_iterator it = vec_inliers.begin();
  for(; it != vec_inliers.end(); it++)
    good_match.push_back(vec_matchings[*it]);
  if(! Match::saveMatch(argv[4], good_match)) {
    std::cerr << "Failed saving matches into " <<argv[4] <<std::endl;
    return 1;
  }

  // Sift de-duplicated output display
  if(argc>6)
    homography_matches_output(image1Gray, image2Gray,
                              vec_matchings, vec_inliers,
                              ok? &H:0,
                              argv[5], argv[6]);

  if(! ok) {
    std::cerr << "Failed to estimate a model" << std::endl;
    return 1;
  }

  // Mosaics
  if(argc>7) {
    Rect intersection;
    if(IntersectionBox(w1, h1, w2, h2, H, intersection) &&
       intersection.Width() > 0 && intersection.Height() > 0) {
      const char *fileReg1=0, *fileReg2=0;
      if(argc>9) { fileReg1=argv[8]; fileReg2=argv[9]; }
      homography_registration_output(image1, image2, H, intersection,
                                     argv[7], fileReg1, fileReg2);
    }
  }
  return 0;
}

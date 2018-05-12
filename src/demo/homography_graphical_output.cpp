/**
 * @file homography_graphical_output.cpp
 * @brief Graphical output to show homography estimation
 * @author Pascal Monasse, Pierre Moulon
 * 
 * Copyright (c) 2011-2018 Pascal Monasse, Pierre Moulon
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

#include "homography_graphical_output.hpp"
#include "graphical_output.hpp"
#include <algorithm>
#include <iostream>
#include "libOrsa/libNumerics/homography.h"
#include "libImage/image_drawing.hpp"
#include "libImage/image_io.hpp"
#include "warping.hpp"
#include "Rect.hpp"

using namespace libNumerics;

RGBColor COL_ERRH=YELLOW;

/// Return 3x3 translation matrix
libNumerics::matrix<double> translation(double dx, double dy) {
  libNumerics::matrix<double> T = libNumerics::matrix<double>::eye(3);
  T(0,2) = dx;
  T(1,2) = dy;
  return T;
}

/// Output inliers and outliers in image files \a fileIn and \a fileOut.
void homography_matches_output(const Image<unsigned char>& image1,
                               const Image<unsigned char>& image2,
                               const std::vector<Match>& vec_all,
                               std::vector<int> vec_in,
                               const libNumerics::matrix<double>* H,
                               const char* fileIn, const char* fileOut) {
  matrix<double> T1(3,3), T2(3,3);
  Image<RGBColor> in = concat_images(image1,image2,&T1,&T2);
  Image<RGBColor> out(in);

  // List of outliers
  std::vector<int> vec_out;
  complement(vec_all.size(), vec_in, &vec_out);

  // For outliers, show error vector from prediction to observation
  std::vector<int>::const_iterator it;
  if(H) // Otherwise, no prediction
    for(it=vec_out.begin(); it!=vec_out.end(); ++it)
      draw_match(vec_all.at(*it), COL_ERRH, &out, T2**H, T2);

  // Draw links for inliers and outliers in respective image
  for(it=vec_out.begin(); it!=vec_out.end(); ++it)
    draw_match(vec_all.at(*it), COL_OUT_LINK, &out, T1, T2);
  for(it=vec_in.begin(); it!=vec_in.end(); ++it)
    draw_match(vec_all.at(*it), COL_IN_LINK, &in, T1, T2);

  libs::WriteImage(fileIn,  in);
  libs::WriteImage(fileOut, out);
}

/// Output mosaic image and registered images to files.
void homography_registration_output(const Image<RGBColor>& image1,
                                    const Image<RGBColor>& image2,
                                    const libNumerics::matrix<double>& H,
                                    const Rect& rect,
                                    const char* fileMosaic,
                                    const char* fileReg1,
                                    const char* fileReg2) {
  int xc=(rect.left+rect.right)/2;
  int yc=(rect.top+rect.bottom)/2;
  size_t wM = std::max(image1.Width(), image2.Width());
  size_t hM = std::max(image1.Height(), image2.Height());
  int xo=static_cast<int>(wM/2);
  int yo=static_cast<int>(hM/2);
  libNumerics::matrix<double> T = translation(xo-xc, yo-yc);

  if(fileMosaic)
  {
    std::cout << "-- Render Mosaic -- " << std::endl;
    Image<RGBColor> imageMosaic(wM, hM);
    imageMosaic.fill(WHITE);
    Warp(image1, T*H, image2, T, imageMosaic);
    libs::WriteImage(fileMosaic, imageMosaic);
  }
  if(fileReg1)
  {
    std::cout << "-- Render Mosaic - Image 1 -- " << std::endl;
    Image<RGBColor> reg1(wM, hM, WHITE);
    Warp(image1, T*H, reg1);
    libs::WriteImage(fileReg1, reg1);
  }
  if(fileReg2)
  {
    std::cout << "-- Render Mosaic - Image 2 -- " << std::endl;
    Image<RGBColor> reg2(wM, hM, WHITE);
    Warp(image2, T, reg2);
    libs::WriteImage(fileReg2, reg2);
  }
}

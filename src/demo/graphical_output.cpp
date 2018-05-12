/**
 * @file graphical_output.cpp
 * @brief Graphical output to show fundamental matrix estimation
 * @author Lionel Moisan, Pascal Monasse, Pierre Moulon
 * 
 * Copyright (c) 2016-2018 Lionel Moisan, Pascal Monasse, Pierre Moulon
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

#include "graphical_output.hpp"
#include "warping.hpp"
#include "libImage/image_drawing.hpp"
#include "libImage/image_io.hpp"
#include <algorithm>
#include <numeric>

using namespace libNumerics;

/// Default colors
RGBColor COL_IN=GREEN;
RGBColor COL_IN_LINK=COL_IN;
RGBColor COL_OUT=RED;
RGBColor COL_OUT_LINK=COL_OUT;

/// Complement of values of \a in in range 0,1,...,size-1.
void complement(int size, std::vector<int>& in, std::vector<int>* out) {
  std::sort(in.begin(), in.end());
  // C++11: replace following three lines with std::iota
  std::vector<int> all(size, 1);
  if(! all.empty()) all[0]=0;
  std::partial_sum(all.begin(), all.end(), all.begin());
  std::set_difference(all.begin(), all.end(), in.begin(), in.end(),
                      std::back_inserter(*out));
}

/// Return 3x3 zoom-translation matrix
static matrix<double> zoomtrans(double z, double dx, double dy) {
    matrix<double> T = matrix<double>::eye(3);
    T(0,0)=T(1,1) = z;
    T(0,2) = dx;
    T(1,2) = dy;
    return T;
}

/// Concatenate images side by side.
/// \param image1,image2 Left and right images.
/// \param[out] T1,T2 Transforms coords from each image to concatenation.
/// \param horizontalLayout Layout images horizontally or vertically?
/// \param dim Dimension of output image (horizontal or vertical)
/// \return Concatenation image.
/// A negative value of dim means to keep the maximum dimension of the images.
Image<RGBColor> concat_images(const Image<unsigned char>& image1,
                              const Image<unsigned char>& image2,
                              matrix<double>* T1,
                              matrix<double>* T2,
                              bool horizontalLayout, int dim) {
  int w1 = image1.Width(), h1 = image1.Height();
  int w2 = image2.Width(), h2 = image2.Height();
  int w=std::max(w1,w2), h=std::max(h1,h2);
  if(dim>0) { if(horizontalLayout) w=dim; else h=dim; }
  float z = horizontalLayout? w/(float)(w1+w2): h/(float)(h1+h2);
  int wc=w, hc=h;
  if(horizontalLayout) {
    hc = int(z*h);
    *T1 = zoomtrans(z, 0,   (hc-h1*z)/2);
    *T2 = zoomtrans(z, w1*z,(hc-h2*z)/2);
  } else {
    wc = int(z*w);
    *T1 = zoomtrans(z, (wc-w1*z)/2, 0);
    *T2 = zoomtrans(z, (wc-w2*z)/2, h1*z);
  }
  Image<unsigned char> concat(wc, hc, 255);
  Warp(image1, *T1, image2, *T2, concat);

  // Split
  double xSplit1=0,ySplit1=0, xSplit2=0,ySplit2=image2.Height();
  TransformH(*T2, xSplit1,ySplit1);
  TransformH(*T2, xSplit2,ySplit2);
  Image<RGBColor> im;
  libs::convertImage(concat, &im);
  libs::DrawLine(xSplit1,ySplit1, xSplit2,ySplit2, BLUE, &im);

  return im;
}

/// Draw line or endpoints of the match in concatenated image.
void draw_match(const Match& m,  RGBColor col, Image<RGBColor>* im,
                const matrix<double>& T1, const matrix<double>& T2, bool link) {
    double x1=m.x1, y1=m.y1, x2=m.x2, y2=m.y2;
    TransformH(T1, x1,y1);
    TransformH(T2, x2,y2);
    if(link)
        libs::DrawLine((int)x1,(int)y1,(int)x2, (int)y2, col, im);
    else {
        libs::DrawCircle((int)x1,(int)y1, 2, col, im);
        libs::DrawCircle((int)x2,(int)y2, 2, col, im);
    }
}


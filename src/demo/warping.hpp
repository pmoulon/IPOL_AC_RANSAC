/**
 * @file warping.hpp
 * @brief Warp images
 * @author Pascal Monasse, Pierre Moulon
 * 
 * Copyright (c) 2011-2012 Pascal Monasse, Pierre Moulon
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

#ifndef WARPING_H
#define WARPING_H

#include <limits>
#include "libOrsa/libNumerics/numerics.h"
#include "libImage/sample.hpp"
#include "demo/Rect.hpp"

/// Apply homography transform.
/// Indicate if \a H is orientation preserving around the point.
bool TransformH(const libNumerics::matrix<double> &H, double &x, double &y)
{
  libNumerics::vector<double> X(3);
  X(0)=x; X(1)=y; X(2)=1.0;
  X = H*X;
  bool positive = (X(2)*H(2,2)>0);
  X /= X(2);
  x = X(0); y = X(1);
  return positive;
}

// Compute the common area of warped by homography image1 and image2.
bool IntersectionBox(int w1, int h1, int w2, int h2,
                     const libNumerics::matrix<double>& H, Rect &inter)
{
  int xCoord[4] = {0, w1-1, w1-1,    0};
  int yCoord[4] = {0,    0, h1-1, h1-1};

  Rect rect1(numeric_limits<int>::max(),
             numeric_limits<int>::max(),
             numeric_limits<int>::min(),
             numeric_limits<int>::min());
  for(int i=0; i<4; ++i)
  {
    double xT=xCoord[i], yT=yCoord[i];
    TransformH(H, xT, yT);
    rect1.growTo(xT,yT);
  }

  Rect rect2(0,0,w2-1,h2-1);
  return rect2.intersect(rect1, inter);
}

/// Warp an image given a homography
template <class Image>
void Warp(const Image &im, libNumerics::matrix<double> H, Image &out)
{
  const int wOut=static_cast<int>(out.Width());
  const int hOut=static_cast<int>(out.Height());
  H=H.inv();

#ifdef _OPENMP
  #pragma omp parallel for
#endif
  for(int j=0; j<hOut; ++j)
    for(int i=0; i<wOut; ++i)
    {
      double xT=i, yT=j;
      if(TransformH(H, xT, yT) && im.Contains(yT,xT))
        out(j,i) = libs::SampleLinear(im, (float)yT, (float)xT);
    }
}

/// Warp imageA and imageB with homographies.
/// Use backward mapping with bilinear sampling.
/// For destination pixel search which pixels from imA and imB contribute.
/// Perform a mean blending in the overlap zone, transfered original
/// value in the other part.
template<class Image>
void Warp(const Image &imA, libNumerics::matrix<double> HA,
          const Image &imB, libNumerics::matrix<double> HB,
          Image &out)
{
  const int wOut=static_cast<int>(out.Width());
  const int hOut=static_cast<int>(out.Height());
  HA = HA.inv();
  HB = HB.inv();

#ifdef _OPENMP
  #pragma omp parallel for
#endif
  for(int j=0; j < hOut; ++j)
    for(int i = 0; i < wOut; ++i)
    {
      double xA=i, yA=j;
      double xB=i, yB=j;
      bool bAContrib = TransformH(HA, xA, yA) && imA.Contains(yA,xA);
      bool bBContrib = TransformH(HB, xB, yB) && imB.Contains(yB,xB);

      if(bAContrib && bBContrib) //blending of ImageA and ImageB
        out(j,i) = typename Image::Tpixel(
          libs::SampleLinear(imA,(float)yA,(float)xA)*.5f +
          libs::SampleLinear(imB,(float)yB,(float)xB)*.5f);
      else if(bAContrib) //only ImageA
        out(j,i) = libs::SampleLinear(imA,(float)yA,(float)xA);
      else if(bBContrib) //only ImageB
        out(j,i) = libs::SampleLinear(imB,(float)yB,(float)xB);
    }
}

#endif

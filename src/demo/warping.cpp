/**
 * @file warping.cpp
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

#include "warping.hpp"
#include <limits>

/// Compute the common area of warped by homography image1 and image2.
bool IntersectionBox(int w1, int h1, int w2, int h2,
                     const libNumerics::matrix<double>& H, Rect &inter)
{
  int xCoord[4] = {0, w1-1, w1-1,    0};
  int yCoord[4] = {0,    0, h1-1, h1-1};

  Rect rect1(std::numeric_limits<int>::max(),
             std::numeric_limits<int>::max(),
             std::numeric_limits<int>::min(),
             std::numeric_limits<int>::min());
  for(int i=0; i<4; ++i)
  {
    double xT=xCoord[i], yT=yCoord[i];
    libNumerics::TransformH(H, xT, yT);
    rect1.growTo(xT,yT);
  }

  Rect rect2(0,0,w2-1,h2-1);
  return rect2.intersect(rect1, inter);
}

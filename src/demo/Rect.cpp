/**
 * @file Rect.hpp
 * @brief Image rectangle
 * @author Pascal Monasse
 * 
 * Copyright (c) 2013 Pascal Monasse <monasse@imagine.enpc.fr>
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

#include "Rect.hpp"
#include <algorithm>

/// Grow rectangle to include (x,y)
void Rect::growTo(double x, double y) {
  if(x < left)   left   = static_cast<int>( floor(x) );
  if(x > right)  right  = static_cast<int>( ceil(x) );
  if(y < top)    top    = static_cast<int>( floor(y) );
  if(y > bottom) bottom = static_cast<int>( ceil(y) );
}

/// True if the rectangle r2 intersect *this.
bool Rect::intersect(const Rect &r2) const {
  return !(left > r2.right || right < r2.left ||
           top > r2.bottom || bottom < r2.top);
}

/// Return the same as intersect(r2) and return the common intersection zone
bool Rect::intersect(const Rect &r2, Rect &inter) const {
  if( this->intersect(r2) ) {
    inter = Rect( std::max(left, r2.left)  , std::max(top, r2.top),
                  std::min(right, r2.right), std::min(bottom, r2.bottom));
    return true;
  }
  return false;
}


/// Intersection with line of equation aX+bY+c=0.
///
/// Resulting segment has endpoints (left,top) and (right,bottom), which are not
/// sorted anymore.
bool Rect::intersect(double a, double b, double c) {
  // Find half-planes of the four corners
  bool bTL = (a*left +b*top   +c >=0);
  bool bTR = (a*right+b*top   +c >=0);
  bool bBL = (a*left +b*bottom+c >=0);
  bool bBR = (a*right+b*bottom+c >=0);

  double x[2], y[2];
  int n=0;
  if(bTL != bBL) {
    x[n] = left;
    y[n] = -(a*left+c)/b;
    ++n;
  }
  if(bTR != bBR) {
    x[n] = right;
    y[n] = -(a*right+c)/b;
    ++ n;
  }
  if(bTL != bTR && n < 2) {
    x[n] = -(b*top+c)/a;
    y[n] = top;
    ++ n;
  }
  if(bBL != bBR && n < 2) {
    x[n] = -(b*bottom+c)/a;
    y[n] = bottom;
    ++ n;
  }
  if(n == 2) {
    top    = static_cast<int>( y[0] );
    bottom = static_cast<int>( y[1] );
    left   = static_cast<int>( x[0] );
    right  = static_cast<int>( x[1] );
    return true;
  }
  return false;
}

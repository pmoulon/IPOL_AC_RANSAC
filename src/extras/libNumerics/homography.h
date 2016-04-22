/**
 * @file homography.h
 * @brief Homography
 * @author Pascal Monasse
 * 
 * Copyright (c) 2010 Pascal Monasse
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

#ifndef HOMOGRAPHY_H
#define HOMOGRAPHY_H

#include "matrix.h"

namespace libNumerics {

/// Apply homography transform.
/// Indicate if \a H is orientation preserving around the point.
template <typename T>
bool TransformH(const libNumerics::matrix<T> &H, T &x, T &y)
{
  libNumerics::vector<T> X(3);
  X(0)=x; X(1)=y; X(2)=1.0;
  X = H*X;
  bool positive = (X(2)*H(2,2)>0);
  X /= X(2);
  x = X(0); y = X(1);
  return positive;
}

} // namespace libNumerics

#endif

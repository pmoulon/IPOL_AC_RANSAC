/**
 * @file Rect.hpp
 * @brief Image rectangle
 * @author Lionel Moisan, Pascal Monasse, Pierre Moulon
 * 
 * Copyright (c) 2011 Lionel Moisan, Pascal Monasse, Pierre Moulon
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

#ifndef LIBS_RECT_H
#define LIBS_RECT_H

#include <cmath>

/// Minimal class for a rectangle
class Rect
{
  public:
    int top, bottom;
    int left, right;

    /// Rectangle constructor
    Rect(int l=0, int t=0, int r=0, int b=0)
    : top(t), bottom(b), left(l), right(r)
    {}

    void growTo(double x, double y);

    bool intersect(const Rect &r2) const;
    bool intersect(const Rect &r2, Rect &inter) const;
    bool intersect(double a, double b, double c);

    /// Return width of the rectangle
    int Width() const {return right - left; }
    /// Return height of the rectangle
    int Height()const {return bottom - top; }
};

#endif //LIBS_RECT_H

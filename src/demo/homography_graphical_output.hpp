/**
 * @file homography_graphical_output.hpp
 * @brief Graphical output to show homography estimation
 * @author Pascal Monasse, Pierre Moulon
 * 
 * Copyright (c) 2011-2016 Pascal Monasse, Pierre Moulon
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

#ifndef HOMOGRAPHY_GRAPHICAL_OUTPUT_H
#define HOMOGRAPHY_GRAPHICAL_OUTPUT_H

#include <vector>
#include "extras/libMatch/match.h"
#include "extras/libNumerics/matrix.h"
#include "libImage/image.hpp"
#include "libImage/pixelTypes.hpp"
class Rect;

void homography_matches_output(const Image<unsigned char>& image1,
                               const Image<unsigned char>& image2,
                               const std::vector<Match>& vec_all,
                               const std::vector<int>& vec_in,
                               const libNumerics::matrix<double>* H,
                               const char* fileIn, const char* fileOut);

void homography_registration_output(const Image<RGBColor>& image1,
                                    const Image<RGBColor>& image2,
                                    const libNumerics::matrix<double>& H,
                                    const Rect& rect,
                                    const char* fileMosaic,
                                    const char* fileReg1,
                                    const char* fileReg2);

#endif

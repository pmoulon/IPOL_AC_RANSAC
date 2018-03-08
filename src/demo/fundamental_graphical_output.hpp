/**
 * @file fundamental_graphical_output.hpp
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

#ifndef FUNDAMENTAL_GRAPHICAL_OUTPUT_H
#define FUNDAMENTAL_GRAPHICAL_OUTPUT_H

#include <vector>
#include "libOrsa/match.hpp"
#include "libOrsa/libNumerics/matrix.h"
#include "libImage/image.hpp"
#include "libImage/pixelTypes.hpp"
#include "Rect.hpp"

/// Global colors, having default values (green, red, yellow). Can be changed.
extern RGBColor COL_IN,COL_IN_LINK, COL_OUT,COL_OUT_LINK, COL_EPI,COL_EPI_DIST;

Image<RGBColor> concat_images(const Image<unsigned char>& image1,
                              const Image<unsigned char>& image2,
                              libNumerics::matrix<double>* T1,
                              libNumerics::matrix<double>* T2,
                              bool horizontalLayout=true);

void draw_match(const Match& m,  RGBColor col, Image<RGBColor>* im,
                const libNumerics::matrix<double>& T1,
                const libNumerics::matrix<double>& T2,
                bool link=true);

bool draw_epi(const libNumerics::matrix<double>& F, const Match& m,
              const libNumerics::matrix<double>& T, double halfLength,
              Image<RGBColor>* out, RGBColor col, const RGBColor* colDist=0);
bool draw_epi(const libNumerics::matrix<double>& F, const Match& m,
              const libNumerics::matrix<double>& T, Rect R,
              Image<RGBColor>* out, RGBColor col, const RGBColor* colDist=0);

void fundamental_graphical_output(const Image<unsigned char>& image1,
                                  const Image<unsigned char>& image2,
                                  const std::vector<Match>& vec_all,
                                  std::vector<int> vec_in,
                                  const libNumerics::matrix<double>* F,
                                  const char* fileIn, const char* fileOut,
                                  const char* fileEpi);

#endif

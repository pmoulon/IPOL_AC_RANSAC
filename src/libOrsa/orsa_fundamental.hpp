/**
 * @file orsa_fundamental.hpp
 * @brief Fundamental matrix estimation with ORSA algorithm
 * @author Lionel Moisan, Pascal Monasse, Pierre Moulon
 * 
 * Copyright (c) 2011-2016 Lionel Moisan, Pascal Monasse, Pierre Moulon
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

#ifndef ORSA_FUNDAMENTAL_H
#define ORSA_FUNDAMENTAL_H

#include <vector>
#include "libNumerics/matrix.h"

class Match;

namespace orsa {
/// ORSA fundamental matrix estimation.
bool orsa_fundamental(const std::vector<Match>& vec_matchings,
                      int w1,int h1, int w2,int h2,
                      double precision, int nbIter,
                      libNumerics::matrix<double>& F,
                      std::vector<int>& vec_inliers);
}

#endif

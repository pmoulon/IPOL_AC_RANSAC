/**
 * @file sampling.hpp
 * @brief Uniform sampling.
 * @author Lionel Moisan, Pascal Monasse, Pierre Moulon
 *
 * Copyright (c) 2007 Lionel Moisan
 * Copyright (c) 2010-2011,2020 Pascal Monasse
 * Copyright (c) 2010-2011 Pierre Moulon
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

#ifndef SAMPLING_H
#define SAMPLING_H

#include <vector>

namespace orsa {

/// Sample a set of  \a sizeSample indices in range 0:size-1.
void UniformSample(int sizeSample, int size, std::vector<int> *sample);

/// Sample a set of \a sizeSample indices among \a vec_index.
void UniformSample(int sizeSample,
                   const std::vector<int> &vec_index,
                   std::vector<int> *sample);

}

#endif

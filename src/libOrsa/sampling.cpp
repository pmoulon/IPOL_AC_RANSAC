/**
 * @file sampling.cpp
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

#include "sampling.hpp"
#include <cstdlib>

namespace orsa {

/// Get a (sorted) random sample of size X in [0:n-1]
static void random_sample(std::vector<int> &k, int X, int n)
{
  for(int i=0; i<X; i++) {
    int r = (rand()>>3)%(n-i), j;
    for(j=0; j<i && r>=k[j]; j++)
      r++;
    int j0 = j;
    for(j=i; j > j0; j--)
      k[j]=k[j-1];
    k[j0] = r;
  }
}

/// Sample a set of  \a sizeSample indices in range 0:size-1.
void UniformSample(int sizeSample, int size, std::vector<int> *sample) {
  sample->clear();
  sample->resize(sizeSample);
  random_sample(*sample, sizeSample, size);
}

/// Sample a set of \a sizeSample indices among \a vec_index.
void UniformSample(int sizeSample,
                   const std::vector<int> &vec_index,
                   std::vector<int> *sample) {
  UniformSample(sizeSample, static_cast<int>(vec_index.size()),sample);
  for(int i = 0; i < sizeSample; ++i)
    (*sample)[i] = vec_index[ (*sample)[i] ];
}

}

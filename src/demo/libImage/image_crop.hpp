//Copyright (C) 2012 Pascal Monasse
//
//This program is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//This program is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with this program.  If not, see <http://www.gnu.org/licenses/>.

#ifndef LIBS_IMAGE_IMAGE_CROP_H_
#define LIBS_IMAGE_IMAGE_CROP_H_

#include "libImage/image.hpp"

/// Horizontal concatenation of images
template < class Image >
void Crop(const Image & image, int x0, int y0, int w, int h, Image & Out)
{
  Out = Image(w,h);

  // Fill with original data from imageA.
  for(size_t i = 0; i < (size_t)h; ++i)
    for(size_t j = 0; j < (size_t)w; ++j)
      Out(i,j) = image(y0+i,x0+j);
}

#endif // LIBS_IMAGE_IMAGE_CROP_H_

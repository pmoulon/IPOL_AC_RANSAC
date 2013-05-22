//Copyright (C) 2011 Pierre Moulon, Pascal Monasse
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

#ifndef LIBS_IMAGE_IMAGE_CONVERTER_H
#define LIBS_IMAGE_IMAGE_CONVERTER_H

#include "libImage/image.hpp"
#include "libImage/pixelTypes.hpp"
#include <cassert>

namespace libs{

template<typename Tin, typename Tout>
inline void convert(const Tin& in, Tout& out)
{
  out = static_cast<Tout>(in);
}

template<>
inline void convert<RGBColor, unsigned char>(const RGBColor& in, unsigned char& out)
{
  out = static_cast<unsigned char>(0.2127*in.r + 0.7152*in.g + 0.0722*in.b);
}

template<>
inline void convert<unsigned char,RGBColor>(const unsigned char& in, RGBColor& out)
{
  out.r = out.g = out.b = in;
}

template<typename ImageIn, typename ImageOut>
void convertImage(const ImageIn & imaIn, ImageOut * imaOut) {

  (*imaOut) = ImageOut(imaIn.Width(), imaIn.Height());
  for(size_t j = 0; j < imaIn.Height(); ++j)
    for(size_t i = 0; i < imaIn.Width(); ++i)  {
      convert(imaIn(j,i), (*imaOut)(j,i));
    }
}

} // namespace libs

#endif  // LIBS_IMAGE_IMAGE_CONVERTER_H

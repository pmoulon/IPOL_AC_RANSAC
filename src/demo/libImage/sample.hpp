// Copyright (c) 2007-2011 libmv authors.
// 
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to
// deal in the Software without restriction, including without limitation the
// rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
// sell copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
// 
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
// IN THE SOFTWARE.
	
#ifndef LIBS_IMAGE_SAMPLE_H_
#define LIBS_IMAGE_SAMPLE_H_

#include "libImage/image.hpp"
#include <cmath>

namespace libs {

/// Nearest neighbor interpolation.
template<typename T>
inline T SampleNearest(const Image<T> &image,
                       float y, float x) {
  const int i = int(round(y));
  const int j = int(round(x));
  return image(i, j);
}

static inline void LinearInitAxis(float fx, int width,
                                  int *x1, int *x2,
                                  float *dx1, float *dx2) {
  const int ix = int(fx);
  if (ix < 0) {
    *x1 = 0;
    *x2 = 0;
    *dx1 = 1;
    *dx2 = 0;
  } else if (ix > width-2) {
    *x1 = width-1;
    *x2 = width-1;
    *dx1 = 1;
    *dx2 = 0;
  } else {
    *x1 = ix;
    *x2 = *x1 + 1;
    *dx1 = *x2 - fx;
    *dx2 = 1 - *dx1;
  }
}

/// Linear interpolation.
template<typename T>
inline T SampleLinear(const Image<T> &image, float y, float x) {
  int x1, y1, x2, y2;
  float dx1, dy1, dx2, dy2;

  LinearInitAxis(y, int(image.Height()), &y1, &y2, &dy1, &dy2);
  LinearInitAxis(x, int(image.Width()),  &x1, &x2, &dx1, &dx2);

  const T im11 = image(y1, x1);
  const T im12 = image(y1, x2);
  const T im21 = image(y2, x1);
  const T im22 = image(y2, x2);

  return T(( im11 * dx1 + im12 * dx2) * dy1 +
           ( im21 * dx1 + im22 * dx2) * dy2 );
}

/// Linear interpolation.
/// RGBColor specialization.
template<>
inline RGBColor SampleLinear<RGBColor>(const Image<RGBColor> &image, float y, float x) {
  int x1, y1, x2, y2;
  float dx1, dy1, dx2, dy2;

  LinearInitAxis(y, int(image.Height()), &y1, &y2, &dy1, &dy2);
  LinearInitAxis(x, int(image.Width()),  &x1, &x2, &dx1, &dx2);

  const RGBColor im11rgb = image(y1, x1);
  const RGBColor im12rgb = image(y1, x2);
  const RGBColor im21rgb = image(y2, x1);
  const RGBColor im22rgb = image(y2, x2);

  return RGBColor(
    (unsigned char)(((float)im11rgb.r * dx1 + (float)im12rgb.r * dx2) * dy1 + dy2 * ((float)im21rgb.r * dx1 + (float)im22rgb.r * dx2)),
    (unsigned char)(((float)im11rgb.g * dx1 + (float)im12rgb.g * dx2) * dy1 + dy2 * ((float)im21rgb.g * dx1 + (float)im22rgb.g * dx2)),
    (unsigned char)(((float)im11rgb.b * dx1 + (float)im12rgb.b * dx2) * dy1 + dy2 * ((float)im21rgb.b * dx1 + (float)im22rgb.b * dx2)));
}

}  // namespace libs

#endif  // LIBS_IMAGE_SAMPLE_H_

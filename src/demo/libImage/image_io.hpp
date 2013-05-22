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

#ifndef LIBS_IMAGE_IMAGE_IMAGE_IO_H
#define LIBS_IMAGE_IMAGE_IMAGE_IO_H

#include "libImage/image.hpp"
#include "libImage/image_converter.hpp"
#include "libImage/pixelTypes.hpp"
#include <vector>

namespace libs {

typedef unsigned char uchar;
typedef Image<uchar> ByteImage;

enum Format {
  Pnm, Png, Jpg, Unknown
};

Format GetFormat(const char *c);

/// Try to load the given file in the image<T> image.
template<class T>
int ReadImage(const char *, Image<T> *);

/// Open an png image with unsigned char as memory target.
/// The memory point must be null as input.
int ReadImage(const char *, std::vector<unsigned char> *, int * w, int * h, int * depth);
int ReadPng(const char *, std::vector<unsigned char> *, int * w, int * h, int * depth);
int ReadPngStream(FILE *, std::vector<unsigned char> *, int * w, int * h, int * depth);

template<class T>
int WriteImage(const char *, const Image<T> &);

int WriteImage(const char *, const std::vector<unsigned char> & array, int w, int h, int depth);
int WritePng(const char *, const std::vector<unsigned char> & array, int w, int h, int depth);
int WritePngStream(FILE *,  const std::vector<unsigned char> & array, int w, int h, int depth);

template<class T>
int WriteJpg(const char *, const Image<T> &, int quality=90);
int WriteJpg(const char *, const std::vector<unsigned char> & array, int w, int h, int depth, int quality=90);
int WriteJpgStream(FILE *, const std::vector<unsigned char> & array, int w, int h, int depth, int quality=90);


/// Open a jpg image with unsigned char as memory target.
/// The memory point must be null as input.
int ReadJpg(const char *, std::vector<unsigned char> *, int * w, int * h, int * depth);
int ReadJpgStream(FILE *, std::vector<unsigned char> *, int * w, int * h, int * depth);

int ReadPnm(const char *, std::vector<unsigned char> *, int * w, int * h, int * depth);
int ReadPnmStream(FILE *, std::vector<unsigned char> *, int * w, int * h, int * depth);

int WritePnm(const char *, const std::vector<unsigned char> & array, int w, int h, int depth);
int WritePnmStream(FILE *,  const std::vector<unsigned char> & array, int w, int h, int depth);

template<>
inline int ReadImage(const char * path, Image<unsigned char> * im)
{
  std::vector<unsigned char> ptr;
  int w, h, depth;

  int res = ReadImage(path, &ptr, &w, &h, &depth);

  if (res == 1) {
    if(depth == 1)
    {
      (*im) = Image<unsigned char>(w, h, &ptr[0]);
    } else if(depth == 3)
    {
      Image<RGBColor> color(w, h, (const RGBColor*) &ptr[0]);
      convertImage(color, im);
    } else
      res = 0; // Do not know how to convert to gray
  }
  return res;
}

template<>
inline int ReadImage(const char * path, Image<RGBColor> * im)
{
  std::vector<unsigned char> ptr;
  int w, h, depth;

  int res = ReadImage(path, &ptr, &w, &h, &depth);

  if (res == 1) {
    if(depth == 3)
    {
      (*im) = Image<RGBColor>(w, h, (const RGBColor*) &ptr[0]);
    } else if(depth == 1)
    {
      Image<unsigned char> gray(w, h, &ptr[0]);
      convertImage(gray, im);
    } else
      res = 0; // Do not know how to convert to color
  }
  return res;
}

//--------
//-- Image Writing
//--------

template<class T>
int WriteImage(const char * filename, const Image<T> & im)
{
  const unsigned char * ptr = (unsigned char*)(im.data());
  int depth = sizeof( T ) / sizeof(unsigned char);
  std::vector<unsigned char> array( ptr , ptr + im.Width()*im.Height()*depth );
  const int w = static_cast<int>( im.Width() );
  const int h = static_cast<int>( im.Height()); 

  int res = WriteImage(filename, array, w, h, depth);
  
  return res;
}

template<class T>
int WriteJpg(const char * filename, const Image<T> & im, int quality)
{
  const unsigned char * ptr = (unsigned char*)(im.data());
  const int w = static_cast<int>( im.Width() );
  const int h = static_cast<int>( im.Height());
  const int depth = sizeof( T ) / sizeof(unsigned char);

  std::vector<unsigned char> array( ptr , ptr + w*h*depth );
  

  int res = WriteJpg(filename, array, w, h, depth, quality);

  return res;
}

}  // namespace libs

#endif  // LIBS_IMAGE_IMAGE_IMAGE_IO_H

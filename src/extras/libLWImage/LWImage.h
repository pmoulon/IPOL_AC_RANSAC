/**
 * @file LWImage.h
 * @brief Light weight image structure
 * @author Pascal Monasse
 *
 * Copyright (c) 2010 Pascal Monasse
 * All rights reserved.
 *
 * This program is free software: you can use, modify and/or
 * redistribute it under the terms of the simplified BSD License. You
 * should have received a copy of this license along this program. If
 * not, see <http://www.opensource.org/licenses/bsd-license.html>.
 */

#ifndef LWIMAGE_H
#define LWIMAGE_H

#include <cstdlib>
#include <string.h>

/// A lightweight image container around an array. It is not responsible for
/// allocation/free of the array.
template <typename T>
class LWImage
{
public:
    LWImage();
    LWImage(T* i_data, int i_w, int i_h, int i_comps=1);
    ~LWImage() {}

    bool valid(int x, int y) const;
    int sizeBuffer() const { return w*h*comps; }
    int step() const;
    int stepComp() const;
    T* pixel(int x, int y);
    const T* pixel(int x, int y) const;
    const T* pixel_ext(int x, int y) const;

    T* data; ///< Array of pixels.
    int w; ///< Width of image.
    int h; ///< Height of image.
    int comps; ///< Components per pixel.
    bool planar; ///< Are components separated (true) or contiguous (false)
};

/// Utility function, avoiding the need to precise the type:
/// \code
///   void f(LWImage<unsigned char>&);
///   unsigned char* data = new unsigned char[10*10];
///   f( make_image(data, 10, 10) );
/// \endcode
template <typename T>
LWImage<T> make_image(T* data, int w, int h, int comps=1)
{
    return LWImage<T>(data, w, h, comps);
}

/// Do not forget free
template <typename T>
LWImage<T> alloc_image(int w, int h, int comps=1)
{
    return LWImage<T>((T*)malloc(w*h*comps*sizeof(T)), w, h, comps);
}

/// Do not forget free
template <typename T>
LWImage<T> alloc_image(const LWImage<T>& im)
{
    LWImage<T> out = alloc_image<T>(im.w, im.h, im.comps);
    out.planar = im.planar;
    //for(int i=im.sizeBuffer()-1; i>=0; i--)
    //    out.data[i] = im.data[i];
    memcpy(out.data, im.data, im.sizeBuffer()*sizeof(float));
    return out;
}

#include "LWImage.cpp"

#endif

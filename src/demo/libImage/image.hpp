//Copyright (C) 2011 Pierre Moulon
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

#ifndef LIBS_IMAGE_IMAGE_H_
#define LIBS_IMAGE_IMAGE_H_

#include <iostream>
#include <cstring>
using namespace std;

//-- Container for a 2D image
//-- This class ensure that the image have a width and a height
//-- and a 2D array of T data.
//-
//-- Data is saved in row major format
//-- Pixel access is done with operator(y,x)
//  [2/3/2011 pierre MOULON]
//---------------------------

template <typename T>
class Image
{

public:
  typedef T Tpixel;	//-- Store the pixel data type

  //------------------------------
  //-- Image construction method
	Image(){ _data = NULL; _width = _height = 0; }

  Image(size_t width, size_t height, const T val = T())
  {
    _height = height;
    _width = width;
	const size_t n=_width*_height;
    _data = new T[n];
    for(size_t i=0; i<n; i++)
        _data[i]=val;
  }

  Image(size_t width, size_t height, const T * imdata){
    _height = height;
    _width = width;
    _data = new T[_width * _height];
    memcpy( _data, imdata, _width * _height * sizeof(T) );
  }

  Image(const Image<T> & I){
    _data = NULL; _width = _height = 0;
    (*this) = I;
  }

  Image& operator=(const Image<T> & I) {

    if( this != &I)
    {
      _height = I._height;
      _width = I._width;
      if(_data)
         delete [] _data;
      _data = new T[_width * _height];
      memcpy( _data, I._data, _width * _height * sizeof(T) );
    }
    return (*this);
  }

  virtual ~Image(){
    if (_data)
      delete [] _data;
    _width = _height = 0;
  }
  //-- Image construction method
  //------------------------------


  //------------------------------
  //-- accessors/getters methods
  /// Retrieve the width of the image
  size_t Width()  const { return _width; }
  /// Retrieve the height of the image
  size_t Height() const { return _height; }
  /// Return the depth in byte of the pixel (unsigned char will return 1)
  size_t Depth() const  { return sizeof(Tpixel);}

  /// image memory access
  T * data() const {return _data;}
  /// random pixel access
  inline const T & operator()(size_t y, size_t x) const {return _data[ _width * y + x];}
  /// random pixel access
  inline T & operator()(size_t y, size_t x) {return _data[ _width * y + x];}

  inline bool operator==(const Image<T> & img) const
  {
    if(_width == img._width &&
             _height == img._height)
      return (0 == memcmp(_data, img._data, Width()*Height()*sizeof(T)));
    else
      return false;
  }
  //-- accessors/getters methods
  //------------------------------

  /// Tell if a point is inside the image.
  bool Contains(size_t y, size_t x) const {
    return x < _width && y < _height;
  }

    /// Tell if a point is inside the image.
  bool Contains(double y, double x) const {
    return (x>=0 && x<_width && y>=0 && y<_height);
  }

  //-- Utils
  //
  void fill(const T & val)
  {
    memset(_data, val, sizeof(T) * _width * _height);
  }

  friend ostream & operator<<(ostream & os, const Image<T> & im) {
    os << im.Width() << " " << im.Height();
    for (size_t j = 0; j < im.Height(); ++j)
    {
      for (size_t i = 0; i < im.Width(); ++i)
        os << im(j,i) << " ";
      os << endl;
    }
    return os;
  }

protected :
  //-- Internal data :
  // None
  T * _data;
  size_t _width;
  size_t _height;
};

#endif // LIBS_IMAGE_IMAGE_H_

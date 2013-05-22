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

#ifndef LIBS_IMAGE_PIXELTYPES_H_
#define LIBS_IMAGE_PIXELTYPES_H_

#include <iostream>
using namespace std;

/// RGB templated pixel type
template <class T>
class RGBClass
{
public:
	T r, g, b;

  //------------------------------
  //-- construction method
  RGBClass(){ r = g = b = T(0);}
  RGBClass(T red,T green,T blue){
    r = red;
    g = green;
    b = blue;
  }
  //explicit RGBClass(const Base & val) : Base(val) {}
  explicit inline RGBClass(const T t[3]) {
    r = t[0];
    g = t[1];
    b = t[2];
  }
  explicit inline RGBClass(const T val){
    r = g = b = val;
  }

  /// Return gray
  inline operator T() const { return T(0.3*r+0.59*g+0.11*b);}

  friend ostream & operator<<(ostream & os, const RGBClass & col) {
    os << col.r << " " << col.g << " " << col.b;
    return os;
  }

  template<class Z>
  inline RGBClass operator /(const Z & val) const
  {
    return RGBClass(T((Z)(r) / val),
                    T((Z)(g) / val),
                    T((Z)(b) / val) );
  }

  template<class Z>
  inline RGBClass operator *(const Z & val) const
  {
  	return RGBClass(T((Z)r * val),
                    T((Z)g * val),
                    T((Z)b * val) );
  }

  template<class Z>
  inline RGBClass operator +(const Z & val) const
  {
    return RGBClass(T((Z)(r) + val),
                    T((Z)(g) + val),
                    T((Z)(b) + val) );
  }

  inline RGBClass operator +(const RGBClass<T> & val) const
  {
    return RGBClass(T(r + val.r),
                    T(g + val.g),
                    T(b + val.b) );
  }
};
typedef RGBClass<unsigned char> RGBColor;

const RGBColor WHITE(255,255,255);
const RGBColor BLACK(0,0,0);
const RGBColor BLUE(0,0,255);
const RGBColor RED(255,0,0);
const RGBColor GREEN(0,255,0);
const RGBColor YELLOW(255,255,0);
const RGBColor CYAN(0,255,255);
const RGBColor MAGENTA(255,0,255);

#endif

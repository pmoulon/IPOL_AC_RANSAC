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

#include "image.hpp"
#include "pixelTypes.hpp"
#include "testing/testing.h"

#include <iostream>
using namespace std;

TEST(Image, Basis)
{
  //-- Gray(int) Image creation
  Image<int> imaGray(10,10);
  imaGray(1,1) = 1; //-- Pixel modification
  imaGray(2,2) = 2;
  imaGray(5,0) = 2;
  
  cout << imaGray << endl << endl;
  //-- Get raw ptr to image data :
  const int * ptr = imaGray.data();
  ((int*)ptr)[0] = 2;
  fill(((int*)ptr+9*10),((int*)ptr+10*10),2);
  cout << "After" << endl << imaGray;

  // Construction by re-copy
  Image<int> imageGray2(imaGray);

  Image<int> imageGray3;
  imageGray3 = imaGray;
  
  //-- RGB Image creation
  Image<RGBColor> imaRGB(10,10);
  imaRGB(0,0) = RGBColor(0,1,2);
}
TEST(Image, PixelTypes)
{
  RGBColor  a(BLACK);
  // RGBColor  c(0); // Not accepted because can cause bad pixel affectation value (mixed type...)
  // The following issue must used : (at your own risk)
  RGBColor  b(static_cast<unsigned char>(0));
}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */

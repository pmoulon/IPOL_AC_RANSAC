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

#include "libImage/image.hpp"
#include "libImage/pixelTypes.hpp"
#include "libImage/image_drawing.hpp"
#include "testing/testing.h"
using namespace libs;

#define _USE_MATH_DEFINES // For Windows (M_PI)
#include <math.h>



// Horizontal / Vertical scanlines
// Assert that pixels was drawn at the good place
TEST(ImageDrawing, Scanlines) {

  const int w = 10, h = 10;
  Image<unsigned char> image(h,w);
  image.fill(0);

  // horizontal scanline
  //  __________
  //  |         |
  //  |__here___|
  //  |         |
  //  |_________|
  const int y = 5;
  DrawLine( 0, y, w-1, y, 255, &image);
  for(int i=0; i < w; ++i)
    EXPECT_EQ( image(y,i), 255);

  image.fill(0);

  // Vertical scanline
  //  __________
  //  |    h|   |
  //  |    e|   |
  //  |    r|   |
  //  |____e|___|
  const int x = 5;
  DrawLine( x, 0, x, h-1, 255, &image);
  for (int i = 0; i < h; ++i)
    EXPECT_EQ(image(i,y), 255);
}

TEST(ImageDrawing, Scanlines_RGB) {

  const int w = 10, h = 10;
  Image<RGBColor> image(h,w);
  image.fill(RGBColor(BLACK));

  // horizontal scanline
  //  __________
  //  |         |
  //  |__here___|
  //  |         |
  //  |_________|
  const int y = 5;
  DrawLine( 0, y, w-1, y, RGBColor(GREEN), &image);
  for(int i=0; i < w; ++i)
    EXPECT_EQ( image(y,i), RGBColor(GREEN));

  image.fill(RGBColor(BLACK));

  // Vertical scanline
  //  __________
  //  |    h|   |
  //  |    e|   |
  //  |    r|   |
  //  |____e|___|
  const int x = 5;
  DrawLine( x, 0, x, h-1, RGBColor(YELLOW), &image);
  for (int i = 0; i < h; ++i)
    EXPECT_EQ(image(i,y), RGBColor(YELLOW));
}

// Lines with a given angle +/-45°
// Assert that pixels was drawn at the good place
TEST(ImageDrawing, Lines45) {

  const int w = 10, h = 10;
  Image<unsigned char> image(h,w);
  image.fill(0);

  //  _____
  //  |\  |
  //  | \ |
  //  |__\|

  DrawLine(0, 0, w-1, h-1, 255, &image);
  for (int i = 0; i < w; ++i)
    EXPECT_EQ(image(i,i), 255);

  image.fill(0);

  //  _____
  //  |  / |
  //  | /  |
  //  |/___|_
  DrawLine(0, h-1, w-1, 0, 255, &image);
  for (int i = 0; i < h; ++i)
    EXPECT_EQ(image(h-1-i,i), 255);
}

// Draw a circle in an image and assert that all the points are
// at a distance equal to the radius.
TEST(ImageDrawing, Circle) {

  Image<unsigned char> image(10,10);
  image.fill(0);

  const int radius = 3;
  const int x = 5, y = 5;

  DrawCircle(x, y, radius, (unsigned char)255, &image);

  // Distance checking :
  for ( int j = 0; j < image.Height(); ++j)
  for ( int i = 0; i < image.Width(); ++i) {
    if (image(j, i) == 255)  {
      const float distance =  sqrt((float)((j-y)*(j-y) + (i-x)*(i-x)));
      EXPECT_NEAR(radius, distance, 1.0f);
      // Due to discretisation we cannot expect better precision
    }
  }
}

// Draw an ellipse with the two radius equal each other...
// in an image and assert that all the points are
// at a distance equal to the radius.
TEST(ImageDrawing, Ellipse) {

  Image<unsigned char> image(10,10);
  image.fill(0);

  const int radius = 3, angle = 0;
  const int x = 5, y = 5;

  DrawEllipse(x, y, radius, radius, (unsigned char)255, &image, (double)angle);

  // Distance checking :
  for ( int j = 0; j < image.Height(); ++j)
  for ( int i = 0; i < image.Width(); ++i) {
    if (image(j, i) == 255)  {
      const float distance =  sqrt((float)((j-y)*(j-y) + (i-x)*(i-x)));
      EXPECT_NEAR(radius, distance, 1.0f);
      // Due to discretisation we cannot expect better precision
    }
  }
}

// Draw an ellipse with the two radius and rotated ...
// in an image and assert that all the points are
// within the given radius.
TEST(ImageDrawing, RotatedEllipse) {

  Image<unsigned char> image(30,30);
  image.fill(0);

  const int radius = 6;
  const int x = 10, y = 10;

  DrawEllipse(x, y, radius, radius/2, (unsigned char)255, &image, M_PI/4.0);

  // Distance checking :
  for ( int j = 0; j < image.Height(); ++j)
  for ( int i = 0; i < image.Width(); ++i) {
    if (image(j, i) == 255)  {
      const float distance =  sqrt((float)((j-y)*(j-y) + (i-x)*(i-x)));
      EXPECT_EQ( radius+1 >= distance && radius/2.0-1 <= distance, true);
      // Due to discretization we cannot expect better precision
      // Use +-1 to avoid rasterization error.
    }
  }
}

/// Assert that the DrawLine function do not crash
/// when one point is outside the image
TEST(ImageDrawing, DrawLine_PointOutsideTheImage) {

  Image<unsigned char> image(30,30);
  image.fill(0);

  const int radius = 20;
  int x = 15, y = 15;

  // Distance checking :
  for (double i=0; i < 2.0*3.14; i+=3.14/12)
  {
    int x1 = static_cast<int>( cos(i) * radius );
    int y1 = static_cast<int>( sin(i) * radius );
    DrawLine( x, y, x+x1, y+y1, 255, &image);
  }
  // Translate :
  x += 15/2;
  for (double i=0; i < 2.0*3.14; i+=3.14/12)
  {
    int x1 = static_cast<int>( cos(i) * radius );
    int y1 = static_cast<int>( sin(i) * radius );
    DrawLine( x, y, x+x1, y+y1, 255, &image);
  }
  // Translate :
  x += 15/2;
  for (double i=0; i < 2.0*3.14; i+=3.14/12)
  {
    int x1 = static_cast<int>( cos(i) * radius );
    int y1 = static_cast<int>( sin(i) * radius );
    DrawLine( x, y, x+x1, y+y1, 255, &image);
  }

  //Point totally outside the image
  x = y = -100;
  for (double i=0; i < 2.0*3.14; i+=3.14/12)
  {
    int x1 = static_cast<int>( cos(i) * radius );
    int y1 = static_cast<int>( sin(i) * radius );
    DrawLine( x, y, x+x1, y+y1, 255, &image);
  }
  //WriteImage( image, "toto.png");
}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */

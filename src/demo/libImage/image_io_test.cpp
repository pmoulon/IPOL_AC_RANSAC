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

#include <cstdio>
#include <iostream>
#include <string>

#include "libImage/image.hpp"
#include "libImage/pixelTypes.hpp"
#include "libImage/image_io.hpp"
#include "testing/testing.h"

using namespace libs;
using std::string;

TEST(ReadJpg, Jpg_Color) {
  Image<RGBColor> image;
  string jpg_filename = string(THIS_SOURCE_DIR) + "/image_test/two_pixels_color.jpg";
  EXPECT_TRUE(ReadImage(jpg_filename.c_str(), &image));
  EXPECT_EQ(2, image.Width());
  EXPECT_EQ(1, image.Height());
  EXPECT_EQ(3, image.Depth());
  EXPECT_EQ(image(0,0), RGBColor(255, 125, 11));
  EXPECT_EQ(image(0,1), RGBColor( 20, 127, 255));
}

TEST(ReadJpg, Jpg_Monochrome) {
  Image<unsigned char> image;
  string jpg_filename = string(THIS_SOURCE_DIR) + "/image_test/two_pixels_monochrome.jpg";
  EXPECT_TRUE(ReadImage(jpg_filename.c_str(), &image));
  EXPECT_EQ(2, image.Width());
  EXPECT_EQ(1, image.Height());
  EXPECT_EQ(1, image.Depth());
  EXPECT_EQ(image(0,0), (unsigned char)255);
  EXPECT_EQ(image(0,1), (unsigned char)0);
}

TEST(ReadPng, Png_Monochrome) {
  Image<unsigned char> image;
  string png_filename = string(THIS_SOURCE_DIR) + "/image_test/two_pixels_monochrome.png";
  EXPECT_TRUE(ReadImage(png_filename.c_str(), &image));
  EXPECT_EQ(2, image.Width());
  EXPECT_EQ(1, image.Height());
  EXPECT_EQ(1, image.Depth());
  EXPECT_EQ(image(0,0), (unsigned char)255);
  EXPECT_EQ(image(0,1), (unsigned char)0);
}

TEST(GetFormat, filenames) {
  EXPECT_EQ(GetFormat("something.jpg"), libs::Jpg);
  EXPECT_EQ(GetFormat("something.png"), libs::Png);
  EXPECT_EQ(GetFormat("something.pnm"), libs::Pnm);
  EXPECT_EQ(GetFormat("/some/thing.JpG"), libs::Jpg);
  EXPECT_EQ(GetFormat("/some/thing.pNG"), libs::Png);
  EXPECT_EQ(GetFormat("some/thing.PNm"), libs::Pnm);
  EXPECT_EQ(GetFormat(".s/o.m/e.t/h.i/n.g.JPG"), libs::Jpg);
  EXPECT_EQ(GetFormat(".s/o.m/e.t/h.i/n.g.PNG"), libs::Png);
  EXPECT_EQ(GetFormat(".s/o.m/e.t/h.i/n.g.PNM"), libs::Pnm);
}

TEST(ImageIOTest, Png_Out) {
  Image<unsigned char> image(1,2);
  image(0,0) = 255;
  image(1,0) = 0;
  string out_filename = ("test_write_png.png");
  EXPECT_TRUE(WriteImage(out_filename.c_str(), image));

  Image<unsigned char> read_image;
  EXPECT_TRUE(ReadImage(out_filename.c_str(), &read_image));
  EXPECT_TRUE(read_image == image);
  remove(out_filename.c_str());
}

TEST(ImageIOTest, InvalidFiles) {
  Image<unsigned char> image;
  string filename = string(THIS_SOURCE_DIR) + "/donotexist.jpg";
  EXPECT_FALSE(ReadImage(filename.c_str(), &image));
  EXPECT_FALSE(ReadImage("hopefully_unexisting_file", &image));
  remove(filename.c_str());
}

TEST(ImageIOTest, Jpg) {
  Image<unsigned char> image(1,2);
  image(0,0) = 255;
  image(1,0) = 0;
  string filename = ("test_write_jpg.jpg");
  EXPECT_TRUE(WriteJpg(filename.c_str(), image, 100));

  Image<unsigned char> read_image;
  EXPECT_TRUE(ReadImage(filename.c_str(), &read_image));
  EXPECT_TRUE(read_image == image);
  remove(filename.c_str());
}

TEST(ReadPnm, Pgm) {
  Image<unsigned char> image;
  string pgm_filename = string(THIS_SOURCE_DIR) + "/image_test/two_pixels.pgm";
  EXPECT_TRUE(ReadImage(pgm_filename.c_str(), &image));
  EXPECT_EQ(2, image.Width());
  EXPECT_EQ(1, image.Height());
  EXPECT_EQ(1, image.Depth());
  EXPECT_EQ(image(0,0), (unsigned char)255);
  EXPECT_EQ(image(0,1), (unsigned char)0);
}

TEST(ReadPnm, PgmComments) {
  Image<unsigned char> image;
  string pgm_filename = string(THIS_SOURCE_DIR) + "/image_test/two_pixels_gray.pgm";
  EXPECT_TRUE(ReadImage(pgm_filename.c_str(), &image));
  EXPECT_EQ(2, image.Width());
  EXPECT_EQ(1, image.Height());
  EXPECT_EQ(1, image.Depth());
  EXPECT_EQ(image(0,0), (unsigned char)255);
  EXPECT_EQ(image(0,1), (unsigned char)0);
}


TEST(ImageIOTest, Pgm) {
  Image<unsigned char> image(1,2);
  image(0,0) = 255;
  image(1,0) = 0;
  string out_filename = "test_write_pnm.pgm";
  EXPECT_TRUE(WriteImage(out_filename.c_str(),image));

  Image<unsigned char> read_image;
  EXPECT_TRUE(ReadImage(out_filename.c_str(), &read_image));
  EXPECT_TRUE(read_image == image);
  remove(out_filename.c_str());
}

TEST(ReadPnm, Ppm) {
  Image<RGBColor> image;
  string ppm_filename = string(THIS_SOURCE_DIR) + "/image_test/two_pixels.ppm";
  EXPECT_TRUE(ReadImage(ppm_filename.c_str(), &image));
  EXPECT_EQ(2, image.Width());
  EXPECT_EQ(1, image.Height());
  EXPECT_EQ(3, image.Depth());
  EXPECT_EQ(image(0,0), RGBColor( (unsigned char)255));
  EXPECT_EQ(image(0,1), RGBColor( (unsigned char)0));
}

TEST(ImageIOTest, Ppm) {
  Image<RGBColor> image(1,2);
  image(0,0) = RGBColor((unsigned char)255);
  image(1,0) = RGBColor((unsigned char)0);
  string out_filename = "test_write_pnm.ppm";
  EXPECT_TRUE(WriteImage(out_filename.c_str(), image));

  Image<RGBColor> read_image;
  EXPECT_TRUE(ReadImage(out_filename.c_str(), &read_image));
  EXPECT_TRUE(read_image == image);
  remove(out_filename.c_str());
}


/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */

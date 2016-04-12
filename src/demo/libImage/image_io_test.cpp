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
#include "CppUnitLite/TestHarness.h"

using namespace libs;
using std::string;

TEST(ReadJpg, Jpg_Color) {
  Image<RGBColor> image;
  string jpg_filename = string(THIS_SOURCE_DIR) + "/image_test/two_pixels_color.jpg";
  CHECK(ReadImage(jpg_filename.c_str(), &image));
  CHECK_EQUAL(2, image.Width());
  CHECK_EQUAL(1, image.Height());
  CHECK_EQUAL(3, image.Depth());
  CHECK_EQUAL(image(0,0), RGBColor(255, 125, 11));
  CHECK_EQUAL(image(0,1), RGBColor( 20, 127, 255));
}

TEST(ReadJpg, Jpg_Monochrome) {
  Image<unsigned char> image;
  string jpg_filename = string(THIS_SOURCE_DIR) + "/image_test/two_pixels_monochrome.jpg";
  CHECK(ReadImage(jpg_filename.c_str(), &image));
  CHECK_EQUAL(2, image.Width());
  CHECK_EQUAL(1, image.Height());
  CHECK_EQUAL(1, image.Depth());
  CHECK_EQUAL(image(0,0), (unsigned char)255);
  CHECK_EQUAL(image(0,1), (unsigned char)0);
}

TEST(ReadPng, Png_Monochrome) {
  Image<unsigned char> image;
  string png_filename = string(THIS_SOURCE_DIR) + "/image_test/two_pixels_monochrome.png";
  CHECK(ReadImage(png_filename.c_str(), &image));
  CHECK_EQUAL(2, image.Width());
  CHECK_EQUAL(1, image.Height());
  CHECK_EQUAL(1, image.Depth());
  CHECK_EQUAL(image(0,0), (unsigned char)255);
  CHECK_EQUAL(image(0,1), (unsigned char)0);
}

TEST(GetFormat, filenames) {
  CHECK_EQUAL(GetFormat("something.jpg"), libs::Jpg);
  CHECK_EQUAL(GetFormat("something.png"), libs::Png);
  CHECK_EQUAL(GetFormat("something.pnm"), libs::Pnm);
  CHECK_EQUAL(GetFormat("/some/thing.JpG"), libs::Jpg);
  CHECK_EQUAL(GetFormat("/some/thing.pNG"), libs::Png);
  CHECK_EQUAL(GetFormat("some/thing.PNm"), libs::Pnm);
  CHECK_EQUAL(GetFormat(".s/o.m/e.t/h.i/n.g.JPG"), libs::Jpg);
  CHECK_EQUAL(GetFormat(".s/o.m/e.t/h.i/n.g.PNG"), libs::Png);
  CHECK_EQUAL(GetFormat(".s/o.m/e.t/h.i/n.g.PNM"), libs::Pnm);
}

TEST(ImageIOTest, Png_Out) {
  Image<unsigned char> image(1,2);
  image(0,0) = 255;
  image(1,0) = 0;
  string out_filename = ("test_write_png.png");
  CHECK(WriteImage(out_filename.c_str(), image));

  Image<unsigned char> read_image;
  CHECK(ReadImage(out_filename.c_str(), &read_image));
  CHECK(read_image == image);
  remove(out_filename.c_str());
}

TEST(ImageIOTest, InvalidFiles) {
  Image<unsigned char> image;
  string filename = string(THIS_SOURCE_DIR) + "/donotexist.jpg";
  CHECK(! ReadImage(filename.c_str(), &image));
  CHECK(! ReadImage("hopefully_unexisting_file", &image));
  remove(filename.c_str());
}

TEST(ImageIOTest, Jpg) {
  Image<unsigned char> image(1,2);
  image(0,0) = 255;
  image(1,0) = 0;
  string filename = ("test_write_jpg.jpg");
  CHECK(WriteJpg(filename.c_str(), image, 100));

  Image<unsigned char> read_image;
  CHECK(ReadImage(filename.c_str(), &read_image));
  CHECK(read_image == image);
  remove(filename.c_str());
}

TEST(ReadPnm, Pgm) {
  Image<unsigned char> image;
  string pgm_filename = string(THIS_SOURCE_DIR) + "/image_test/two_pixels.pgm";
  CHECK(ReadImage(pgm_filename.c_str(), &image));
  CHECK_EQUAL(2, image.Width());
  CHECK_EQUAL(1, image.Height());
  CHECK_EQUAL(1, image.Depth());
  CHECK_EQUAL(image(0,0), (unsigned char)255);
  CHECK_EQUAL(image(0,1), (unsigned char)0);
}

TEST(ReadPnm, PgmComments) {
  Image<unsigned char> image;
  string pgm_filename = string(THIS_SOURCE_DIR) + "/image_test/two_pixels_gray.pgm";
  CHECK(ReadImage(pgm_filename.c_str(), &image));
  CHECK_EQUAL(2, image.Width());
  CHECK_EQUAL(1, image.Height());
  CHECK_EQUAL(1, image.Depth());
  CHECK_EQUAL(image(0,0), (unsigned char)255);
  CHECK_EQUAL(image(0,1), (unsigned char)0);
}


TEST(ImageIOTest, Pgm) {
  Image<unsigned char> image(1,2);
  image(0,0) = 255;
  image(1,0) = 0;
  string out_filename = "test_write_pnm.pgm";
  CHECK(WriteImage(out_filename.c_str(),image));

  Image<unsigned char> read_image;
  CHECK(ReadImage(out_filename.c_str(), &read_image));
  CHECK(read_image == image);
  remove(out_filename.c_str());
}

TEST(ReadPnm, Ppm) {
  Image<RGBColor> image;
  string ppm_filename = string(THIS_SOURCE_DIR) + "/image_test/two_pixels.ppm";
  CHECK(ReadImage(ppm_filename.c_str(), &image));
  CHECK_EQUAL(2, image.Width());
  CHECK_EQUAL(1, image.Height());
  CHECK_EQUAL(3, image.Depth());
  CHECK_EQUAL(image(0,0), RGBColor( (unsigned char)255));
  CHECK_EQUAL(image(0,1), RGBColor( (unsigned char)0));
}

TEST(ImageIOTest, Ppm) {
  Image<RGBColor> image(1,2);
  image(0,0) = RGBColor((unsigned char)255);
  image(1,0) = RGBColor((unsigned char)0);
  string out_filename = "test_write_pnm.ppm";
  CHECK(WriteImage(out_filename.c_str(), image));

  Image<RGBColor> read_image;
  CHECK(ReadImage(out_filename.c_str(), &read_image));
  CHECK(read_image == image);
  remove(out_filename.c_str());
}


/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */

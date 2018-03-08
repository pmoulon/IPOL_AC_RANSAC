/**
 * @file fundamental_graphical_output.cpp
 * @brief Graphical output to show fundamental matrix estimation
 * @author Lionel Moisan, Pascal Monasse, Pierre Moulon
 * 
 * Copyright (c) 2016-2018 Lionel Moisan, Pascal Monasse, Pierre Moulon
 * All rights reserved.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "fundamental_graphical_output.hpp"

#include <algorithm>
#include <numeric>

#include "demo/warping.hpp"
#include "libImage/image_drawing.hpp"
#include "libImage/image_io.hpp"

using namespace libNumerics;

/// Default colors
RGBColor COL_IN=GREEN;
RGBColor COL_IN_LINK=COL_IN;
RGBColor COL_OUT=RED;
RGBColor COL_OUT_LINK=COL_OUT;
RGBColor COL_EPI=YELLOW;
RGBColor COL_EPI_DIST=YELLOW; ///< Distance to epipolar line

/// Inverse match
inline Match inv(const Match& m) {
    return Match(m.x2,m.y2,m.x1,m.y1);
}

/// Return 3x3 zoom-translation matrix
static matrix<double> zoomtrans(double z, double dx, double dy) {
    matrix<double> T = matrix<double>::eye(3);
    T(0,0)=T(1,1) = z;
    T(0,2) = dx;
    T(1,2) = dy;
    return T;
}

/// Draw line or endpoints of the match in concatenated image.
void draw_match(const Match& m,  RGBColor col, Image<RGBColor>* im,
                const matrix<double>& T1, const matrix<double>& T2, bool link) {
    double x1=m.x1, y1=m.y1, x2=m.x2, y2=m.y2;
    TransformH(T1, x1,y1);
    TransformH(T2, x2,y2);
    if(link)
        libs::DrawLine((int)x1,(int)y1,(int)x2, (int)y2, col, im);
    else {
        libs::DrawCircle((int)x1,(int)y1, 2, col, im);
        libs::DrawCircle((int)x2,(int)y2, 2, col, im);
    }
}

static vector<double> epipolar_line(const matrix<double>& F, const Match& m) {
  vector<double> epi(m.x1, m.y1, 1);
  double qnorm = epi.qnorm();
  epi = F * epi;
  if(epi.qnorm() < 1.0e-11f*qnorm) {
      std::cerr << "Warning: not drawing line close to epipole " << m;
      epi.fill(0);
      return epi;
  }
  return epi / sqrt(epi(0)*epi(0) + epi(1)*epi(1));
}

static void project_on_epi(double& x, double& y, const vector<double>& epi) {
  double a = epi(0), b = epi(1), c = epi(2);
  double lambda = a*x+b*y+c;
  x -= lambda*a;
  y -= lambda*b;
}

/// Draw portion of epipolar line associated to (m.x1,m.y1).
///
/// The center of the segment is the projection of (m.x2,m.y2) on the line.
bool draw_epi(const matrix<double>& F, const Match& m,
              const matrix<double>& T, double halfLength,
              Image<RGBColor>* out, RGBColor col, const RGBColor* colDist)
{
  vector<double> epi=epipolar_line(F, m);
  if(epi.qnorm()==0) return false;
  double xp=m.x2, yp=m.y2;
  project_on_epi(xp,yp,epi);
  if(colDist) {
      double xi=m.x2, yi=m.y2, xf=xp, yf=yp;
      TransformH(T, xi, yi); TransformH(T, xf, yf);
      libs::DrawLine((int)xi,(int)yi,(int)xf,(int)yf, *colDist, out);
  }
  double a = epi(0), b = epi(1);
  double x1=xp-b*halfLength, y1=yp+a*halfLength;
  double x2=xp+b*halfLength, y2=yp-a*halfLength;
  TransformH(T, x1,y1); TransformH(T, x2,y2);
  libs::DrawLine((int)x1, (int)y1, (int)x2, (int)y2, col, out);
  return true;
}

/// Show epipolar line corresponding to (m.x1,m.y1) and orthogonal projection
/// of (m.x2,m.y2) on this line.
///
/// Epipolar line is clipped inside rectangle @R and coordinates are transformed
/// by homography @H before display.
bool draw_epi(const matrix<double>& F, const Match& m,
              const matrix<double>& T, Rect R,
              Image<RGBColor>* out, RGBColor col, const RGBColor* colDist)
{
  vector<double> epi=epipolar_line(F, m);
  if(epi.qnorm()==0) return false;
  double xp=m.x2, yp=m.y2;
  project_on_epi(xp,yp,epi);
  if(colDist) {
      double xi=m.x2, yi=m.y2, xf=xp, yf=yp;
      TransformH(T, xi, yi); TransformH(T, xf, yf);
      libs::DrawLine((int)xi,(int)yi,(int)xf,(int)yf, *colDist, out);
  }
  // Draw epipolar line
  if( R.intersect(epi(0),epi(1),epi(2)) ) {
    double x1=R.left, y1=R.top, x2=R.right, y2=R.bottom;
    TransformH(T, x1, y1);
    TransformH(T, x2, y2);
    libs::DrawLine((int)x1,(int)y1,(int)x2,(int)y2, col, out);
  }
  return true;
}

/// Return image bounding box.
inline Rect rectImage(const Image<unsigned char>& im) {
    return Rect(0,0,im.Width(),im.Height());
}

/// Return length of image diagonal.
inline double diag(const Image<unsigned char>& im) {
    double w=im.Width(), h=im.Height();
    return sqrt(w*w+h*h);
}

/// Concatenate images side by side.
/// \param image1,image2 Left and right images.
/// \param[out] T1,T2 Transforms coords from each image to concatenation.
/// \param horizontalLayout Layout images horizontally or vertically?
/// \return Concanetaion image.
Image<RGBColor> concat_images(const Image<unsigned char>& image1,
                              const Image<unsigned char>& image2,
                              matrix<double>* T1,
                              matrix<double>* T2,
                              bool horizontalLayout) {
  int w1 = image1.Width(), h1 = image1.Height();
  int w2 = image2.Width(), h2 = image2.Height();
  int w=std::max(w1,w2), h=std::max(h1,h2);
  float z = horizontalLayout? w/(float)(w1+w2): h/(float)(h1+h2);
  int wc = w, hc=h;
  if(horizontalLayout) {
    hc = int(z*h);
    *T1 = zoomtrans(z, 0,   (hc-h1*z)/2);
    *T2 = zoomtrans(z, w1*z,(hc-h2*z)/2);
  } else {
    wc = int(z*w);
    *T1 = zoomtrans(z, (wc-w1*z)/2, 0);
    *T2 = zoomtrans(z, (wc-w2*z)/2, h1*z);
  }
  Image<unsigned char> concat(wc, hc, 255);
  Warp(image1, *T1, image2, *T2, concat);

  // Split
  double xSplit1=0,ySplit1=0, xSplit2=0,ySplit2=image2.Height();
  TransformH(*T2, xSplit1,ySplit1);
  TransformH(*T2, xSplit2,ySplit2);
  Image<RGBColor> im;
  libs::convertImage(concat, &im);
  libs::DrawLine(xSplit1,ySplit1, xSplit2,ySplit2, BLUE, &im);

  return im;
}

/// Graphical output for estimated fundamental matrix.
/// \param image1,image2 Left and right images.
/// \param vec_all All correspondences.
/// \param vec_in Index in vec_all of inliers.
/// \param F Estimated fundamental matrix (null pointer if failure).
/// \param fileIn File name of inlier correspondences (can be null).
/// \param fileOut File name of outlier correspondences (can be null).
/// \param fileEpi File name of output epipolar lines (can be null).
void fundamental_graphical_output(const Image<unsigned char>& image1,
                                  const Image<unsigned char>& image2,
                                  const std::vector<Match>& vec_all,
                                  std::vector<int> vec_in,
                                  const matrix<double>* F,
                                  const char* fileIn,
                                  const char* fileOut,
                                  const char* fileEpi)
{
  // List of outliers
  std::vector<int> vec_out;
  {
    std::sort(vec_in.begin(), vec_in.end());
    // C++11: replace following three lines with std::iota
    std::vector<int> all(vec_all.size(), 1);
    if(! all.empty()) all[0]=0;
    std::partial_sum(all.begin(), all.end(), all.begin());
    std::set_difference(all.begin(), all.end(), vec_in.begin(), vec_in.end(),
                        std::back_inserter(vec_out));
  }

  matrix<double> T1(3,3), T2(3,3);
  Image<RGBColor> im = concat_images(image1,image2,&T1,&T2);

  std::vector<int>::const_iterator it;
  if(fileIn) { // Draw links for inliers
    Image<RGBColor> in(im);
    for(it=vec_in.begin(); it!=vec_in.end(); ++it)
      draw_match(vec_all.at(*it), COL_IN_LINK, &in, T1, T2, true);
    libs::WriteImage(fileIn, in);
  }
  if(fileOut) { // Draw epipolar lines and links for outliers
    Image<RGBColor> out(im);
    if(F) {
      Rect R1=rectImage(image1), R2=rectImage(image2);
      for(it=vec_out.begin(); it!=vec_out.end(); ++it) {
        const Match& m=vec_all[*it];
        draw_epi(*F,     inv(m), T1, R1, &out, COL_EPI, &COL_EPI_DIST);
        draw_epi(F->t(),     m,  T2, R2, &out, COL_EPI, &COL_EPI_DIST);
      }
    }
    for(it=vec_out.begin(); it!=vec_out.end(); ++it)
      draw_match(vec_all[*it], COL_OUT_LINK, &out, T1, T2, true);
    libs::WriteImage(fileOut, out);
  }
  if(fileEpi) {
    const double margin1 = diag(image1) *.02;
    const double margin2 = diag(image2) *.02;
    // Draw matches with inlier/outlier color
    for(it=vec_in.begin();  it!=vec_in.end();  ++it)
      draw_match(vec_all.at(*it), COL_IN,  &im, T1, T2, false);
    for(it=vec_out.begin(); it!=vec_out.end(); ++it)
      draw_match(vec_all[*it]   , COL_OUT, &im, T1, T2, false);

    if(F) // Draw small epipolar lines for inliers
      for(it=vec_in.begin();  it!=vec_in.end();  ++it) {
        const Match& m=vec_all.at(*it);
        draw_epi(*F,     inv(m), T1, margin1, &im, COL_EPI);
        draw_epi(F->t(),     m,  T2, margin2, &im, COL_EPI);
      }
    libs::WriteImage(fileEpi, im);
  }
}

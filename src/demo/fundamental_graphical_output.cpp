/**
 * @file fundamental_graphical_output.cpp
 * @brief Graphical output to show fundamental matrix estimation
 * @author Lionel Moisan, Pascal Monasse, Pierre Moulon
 * 
 * Copyright (c) 2016 Lionel Moisan, Pascal Monasse, Pierre Moulon
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

#include "demo/warping.hpp"
#include "libImage/image_drawing.hpp"
#include "libImage/image_io.hpp"

/// Inverse match
inline Match inv(const Match& m)
{
    return Match(m.x2,m.y2,m.x1,m.y1);
}

/// Return 3x3 zoom-translation matrix
libNumerics::matrix<double> zoomtrans(double z, double dx, double dy) {
    libNumerics::matrix<double> T = libNumerics::matrix<double>::eye(3);
    T(0,0)=T(1,1) = z;
    T(0,2) = dx;
    T(1,2) = dy;
    return T;
}

/// Draw endpoints of the match in concatenated image.
void draw_points(const Match& m,  RGBColor col, Image<RGBColor>* im,
                 const libNumerics::matrix<double>& T1,
                 const libNumerics::matrix<double>& T2) {
    double x1=m.x1, y1=m.y1, x2=m.x2, y2=m.y2;
    TransformH(T1, x1,y1);
    TransformH(T2, x2,y2);
    libs::DrawCircle((int)x1,(int)y1, 2, col, im);
    libs::DrawCircle((int)x2,(int)y2, 2, col, im);
}

/// Show epipolar line corresponding to (m.x1,m.y1) and orthogonal projection
/// of (m.x2,m.y2) on this line.
///
/// Epipolar line is clipped inside rectangle @R and coordinates are transformed
/// by homography @H before display.
void display_error(const libNumerics::matrix<double>& F,
                   const Match& m,
                   const libNumerics::matrix<double>& H, Rect R,
                   Image<RGBColor>* out)
{
  static const RGBColor col=YELLOW;
  double x1=m.x1, y1=m.y1, x2=m.x2, y2=m.y2;
  //Epipolar line on second image
  libNumerics::vector<double> epi(x1, y1, 1);
  double qnorm = epi.qnorm();
  epi = F * epi;
  if(epi.qnorm() < 1.0e-5f*qnorm) {
      cerr << "Warning: not drawing line close to epipole " << m;
      return;
  }
  epi /= sqrt(epi(0)*epi(0) + epi(1)*epi(1));

  //Draw segment to projection on epipolar line
  double a = epi(0), b = epi(1), c = epi(2);
  double lambda = a*x2 + b*y2 + c;
  x1 = x2 - lambda*a;
  y1 = y2 - lambda*b;
  TransformH(H, x1, y1);
  TransformH(H, x2, y2);
  libs::DrawLine((int)x1,(int)y1,(int)x2,(int)y2, col, out);

  // Draw epipolar line
  if( R.intersect(a,b,c) ) {
    x1 = R.left; y1 = R.top; x2= R.right; y2 = R.bottom;
    TransformH(H, x1, y1);
    TransformH(H, x2, y2);
    libs::DrawLine((int)x1,(int)y1,(int)x2,(int)y2, col, out);
  }
}

/// Output inlier and oulier matches in image files.
void display_match(const std::vector<Match>& vec_matchings,
                   std::vector<int>& vec_inliers,
                   const libNumerics::matrix<double>* F,
                   const libNumerics::matrix<double>& H1,
                   const libNumerics::matrix<double>& H2,
                   const Rect& R1, const Rect& R2,
                   Image<RGBColor>* in, Image<RGBColor>* out)
{
  std::sort(vec_inliers.begin(), vec_inliers.end());

  // For outliers, show vector (yellow) from prediction to observation
  std::vector<int>::const_iterator it = vec_inliers.begin();
  std::vector<Match>::const_iterator m = vec_matchings.begin();
  if(F) // Otherwise, no prediction
    for(int i=0; m != vec_matchings.end(); ++m, ++i) {
      if(it != vec_inliers.end() && i==*it)
        ++it;
      else { //Outlier
        display_error(*F,     inv(*m), H1, R1, out);  
        display_error(F->t(),     *m,  H2, R2, out);
      }
    }

  // Show link for inliers (green) and outliers (red)
  it = vec_inliers.begin();
  m = vec_matchings.begin();
  for(int i=0; m != vec_matchings.end(); ++m, ++i)
  {
    Image<RGBColor>* im=out;
    RGBColor col=RED;
    if(it != vec_inliers.end() && i==*it) {
      ++it;
      im=in;
      col=GREEN;
    }
    double x1=m->x1, y1=m->y1, x2=m->x2, y2=m->y2;
    TransformH(H1, x1, y1);
    TransformH(H2, x2, y2);
    libs::DrawLine((int)x1,(int)y1,(int)x2, (int)y2, col, im);
  }
}

/// Draw portion of epipolar line associated to (m.x1,m.y1).
///
/// The center of the segment is the projection of (m.x2,m.y2) on the line.
void draw_small_epi(const libNumerics::matrix<double>& F,
                    const Match& m,
                    const libNumerics::matrix<double>& T, double halfLength,
                    Image<RGBColor>* out)
{
  libNumerics::vector<double> epi(m.x1, m.y1, 1);
  double qnorm = epi.qnorm();
  epi = F * epi;
  if(epi.qnorm() < 1.0e-11f*qnorm) {
      cerr << "Warning: not drawing line close to epipole " << m;
      return;
  }
  epi /= sqrt(epi(0)*epi(0) + epi(1)*epi(1));
  double a = epi(0), b = epi(1), c = epi(2);
  double lambda = (a*m.x2+b*m.y2+c);
  double x1 = m.x2 - lambda*a - b*halfLength;
  double x2 = m.x2 - lambda*a + b*halfLength;
  double y1 = m.y2 - lambda*b + a*halfLength;
  double y2 = m.y2 - lambda*b - a*halfLength;
  TransformH(T, x1,y1);
  TransformH(T, x2,y2);
  libs::DrawLine((int)x1, (int)y1, (int)x2, (int)y2, YELLOW, out);
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
                                  const libNumerics::matrix<double>* F,
                                  const char* fileIn,
                                  const char* fileOut,
                                  const char* fileEpi)
{
  int w1 = image1.Width(), h1 = image1.Height();
  int w2 = image2.Width(), h2 = image2.Height();
  int w=std::max(w1,w2), h=std::max(h1,h2);
  float z = w/(float)(w1+w2); //Set width as max of two images
  Image<unsigned char> concat(w, int(z*h), 255);
  libNumerics::matrix<double> T1=zoomtrans(z, 0,   (concat.Height()-h1*z)/2);
  libNumerics::matrix<double> T2=zoomtrans(z, w1*z,(concat.Height()-h2*z)/2);
  Warp(image1, T1, image2, T2, concat);

  Image<RGBColor> im;
  if(fileIn || fileOut)
  {
    libs::convertImage(concat, &im);
    libs::DrawLine(int(w1*z),0, int(w1*z),int(concat.Height()), BLUE, &im);
    Image<RGBColor> out(im);
    Rect R1(0,0,w1,h1), R2(0,0,w2,h2);
    display_match(vec_all,vec_in, F, T1,T2, R1,R2, &im,&out);
    if(fileIn)  libs::WriteImage(fileIn, im);
    if(fileOut) libs::WriteImage(fileOut, out);
  }
  if(fileEpi)
  {
    libs::convertImage(concat, &im);
    libs::DrawLine(int(w1*z),0, int(w1*z),int(concat.Height()), BLUE, &im);
    const double margin1 = sqrt((double)w1*w1+h1*h1) *.02;
    const double margin2 = sqrt((double)w2*w2+h2*h2) *.02;
    // Draw in red outliers
    for(size_t index = 0; index < vec_all.size(); ++index)
      if(std::find(vec_in.begin(),vec_in.end(),index) == vec_in.end())
        draw_points(vec_all[index], RED, &im, T1, T2);
    // Draw in green inliers and small epipolar lines
    for(std::vector<int>::const_iterator it=vec_in.begin();
        it!=vec_in.end(); ++it)
    {
      const Match& m = vec_all[*it];
      draw_points(m, GREEN, &im, T1, T2);
      if(F)
      {
        draw_small_epi(*F,     inv(m), T1, margin1, &im);
        draw_small_epi(F->t(),     m,  T2, margin2, &im);
      }
    }
    libs::WriteImage(fileEpi, im);
  }
}

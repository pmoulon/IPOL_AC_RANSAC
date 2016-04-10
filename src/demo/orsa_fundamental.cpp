/**
 * @file orsa_fundamental.cpp
 * @brief Fundamental matrix estimation
 * @author Lionel Moisan, Pascal Monasse, Pierre Moulon
 * 
 * Copyright (c) 2011-2015 Lionel Moisan, Pascal Monasse, Pierre Moulon
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

#include <cstdlib>
#include <ctime>

#include <iostream>
#include <vector>
#include <algorithm>
#include <utility>

#include "libImage/image.hpp"
#include "libImage/image_io.hpp"
#include "libImage/image_drawing.hpp"

#include "extras/libNumerics/numerics.h"
#include "extras/sift/library.h"

#include "libOrsa/fundamental_model.hpp"

#include "demo/siftMatch.hpp"
#include "demo/warping.hpp"
#include "demo/cmdLine.h"

/// Number of random samples in ORSA
static const int ITER_ORSA=10000;

/// Display average/max error of inliers of fundamental matrix F.
static std::pair<double,double> display_stats(const std::vector<Match>& vec_matchings,
                                              const std::vector<int>& vec_inliers,
                                              const libNumerics::matrix<double>& F) {
    std::vector<int>::const_iterator it=vec_inliers.begin();
    double l2=0, linf=0;
    for(; it!=vec_inliers.end(); ++it) {
      const Match& m=vec_matchings[*it];
      double a = F(0,0) * m.x1 + F(1,0) * m.y1 + F(2,0);
      double b = F(0,1) * m.x1 + F(1,1) * m.y1 + F(2,1);
      double c = F(0,2) * m.x1 + F(1,2) * m.y1 + F(2,2);
      double d = a*m.x2 + b*m.y2 + c;
      double e =  (d*d) / (a*a + b*b);
      l2 += e;
      if(linf < e)
        linf = e;
    }
    std::pair<double,double> err(sqrt(l2/vec_inliers.size()),sqrt(linf));
    std::cout << "Average/max error: " << err.first << "/" << err.second <<std::endl;
    return err;
}

/// ORSA fundamental estimation
bool ORSA(const std::vector<Match>& vec_matchings, int w1,int h1, int w2,int h2,
          double precision,
          libNumerics::matrix<double>& F, std::vector<int>& vec_inliers)
{
  const int n = static_cast<int>( vec_matchings.size() );
  if(n < 7)
  {
    std::cerr << "Error: ORSA needs 7 matches or more to proceed" <<std::endl;
    return false;
  }
  libNumerics::matrix<double> xA(2,n), xB(2,n);

  for (int i=0; i < n; ++i)
  {
    xA(0,i) = vec_matchings[i].x1;
    xA(1,i) = vec_matchings[i].y1;
    xB(0,i) = vec_matchings[i].x2;
    xB(1,i) = vec_matchings[i].y2;
  }

  orsa::FundamentalModel model(xA, w1, h1, xB, w2, h2);
  //model.setConvergenceCheck(true);

  if(model.orsa(vec_inliers, ITER_ORSA, &precision, &F, true)>0.0)
    return false;
  std::pair<double,double> err; // (RMSE,max)
  std::cout << "Before refinement: ";
  err = display_stats(vec_matchings, vec_inliers, F);
  libNumerics::matrix<double> F2(3,3);
  if( model.ComputeModel(vec_inliers,&F2) ) // Re-estimate with all inliers
  {
    std::cout << "After  refinement: ";
    double maxBefore = err.second;
    if(display_stats(vec_matchings, vec_inliers, F2).first <= maxBefore)
      F = F2;
    else
      std::cerr << "Warning: error after refinement is too large, thus ignored" <<std::endl;
  } else
    std::cerr << "Warning: error in refinement, result is suspect" <<std::endl;
  return true;
}

/// Inverse match
Match inv(const Match& m)
{
    return Match(m.x2,m.y2,m.x1,m.y1);
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
  epi = F*epi;
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
  epi = F * epi;
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
  matchingslist::const_iterator m = vec_matchings.begin();
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

/// Return 3x3 zoom-translation matrix
libNumerics::matrix<double> zoomtrans(double z, double dx, double dy) {
    libNumerics::matrix<double> T = libNumerics::matrix<double>::eye(3);
    T(0,0)=T(1,1) = z;
    T(0,2) = dx;
    T(1,2) = dy;
    return T;
}


/// Draw endpoints of the match in concatenated image.
static void draw_points(const Match& m,  RGBColor col, Image<RGBColor>* im,
                        const libNumerics::matrix<double>& T1,
                        const libNumerics::matrix<double>& T2) {
    double x1=m.x1, y1=m.y1, x2=m.x2, y2=m.y2;
    TransformH(T1, x1,y1);
    TransformH(T2, x2,y2);
    libs::DrawCircle((int)x1,(int)y1, 2, col, im);
    libs::DrawCircle((int)x2,(int)y2, 2, col, im);
}


int main(int argc, char **argv)
{
  double precision=0;
  float fSiftRatio=0.6f;
  CmdLine cmd;
  cmd.add( make_option('p',precision, "prec") );
  cmd.add( make_option('s',fSiftRatio, "sift") );
  cmd.add( make_switch('r', "read") );
  try {
    cmd.process(argc, argv);
  } catch(const std::string& s) {
    std::cerr << s << std::endl;
    return 1;
  }
  if(argc!=5 && argc!=6 && argc!=7 && argc!=8) {
    std::cerr << "Usage: " << argv[0] << " imgInA imgInB "
              << "[-p|--prec precision] "
              << "[-s|--sift siftRatio] "
              << "[-r|--read] "
              << "allMatches.txt orsaMatches.txt "
              << "[imgInliers imgOutliers] [imgEpipolar]"
              << std::endl;
    return 1;
  }

  // Init random seed
  srand((unsigned int)time(0));

  Image<RGBColor> image1, image2;
  if(! libs::ReadImage(argv[1], &image1))
    return 1;
  if(! libs::ReadImage(argv[2], &image2))
    return 1;

  Image<unsigned char> image1Gray, image2Gray;
  libs::convertImage(image1, &image1Gray);
  libs::convertImage(image2, &image2Gray);

  // Find matches with SIFT or read correspondence file
  std::vector<Match> vec_matchings;
  if(cmd.used('r')) {
      if(Match::loadMatch(argv[3], vec_matchings))
          std::cout << "Read " <<vec_matchings.size()<< " matches" <<std::endl;
      else {
          std::cerr << "Failed reading matches from " << argv[3] <<std::endl;
          return 1;
      }
  } else
      SIFT(image1Gray, image2Gray, vec_matchings, fSiftRatio);

  // Remove duplicates (frequent with SIFT)
  rm_duplicates(vec_matchings);

  // Save match files
  if(! cmd.used('r') && ! Match::saveMatch(argv[3], vec_matchings)) {
    std::cerr << "Failed saving matches into " <<argv[3] <<std::endl;
    return 1;
  }

  int w1 = image1Gray.Width(), h1 = image1Gray.Height();
  int w2 = image2Gray.Width(), h2 = image2Gray.Height();

  // Estimation of fundamental matrix with ORSA
  libNumerics::matrix<double> F(3,3);
  std::vector<int> vec_inliers;
  bool ok = ORSA(vec_matchings, w1, h1, w2, h2, precision, F, vec_inliers);
  if(ok)
    std::cout << "F=" << F <<std::endl;

  std::vector<Match> good_match;
  std::vector<int>::const_iterator it = vec_inliers.begin();
  for(; it != vec_inliers.end(); it++)
    good_match.push_back(vec_matchings[*it]);
  if(! Match::saveMatch(argv[4], good_match)) {
    std::cerr << "Failed saving matches into " <<argv[4] <<std::endl;
    return 1;
  }

  if(argc>=6) // Output images
  {
    int w=std::max(w1,w2), h=std::max(h1,h2);
    float z = w/(float)(w1+w2); //Set width as max of two images
    Image<unsigned char> concat(w, int(z*h), 255);
    libNumerics::matrix<double> T1=zoomtrans(z, 0,   (concat.Height()-h1*z)/2);
    libNumerics::matrix<double> T2=zoomtrans(z, w1*z,(concat.Height()-h2*z)/2);
    Warp(image1Gray, T1, image2Gray, T2, concat);

    if(argc>6) // Show inliers and outliers
    {
      Image<RGBColor> in;
      libs::convertImage(concat, &in);
      libs::DrawLine(int(w1*z),0, int(w1*z),int(concat.Height()), BLUE, &in);
      Image<RGBColor> out(in);
      Rect R1(0,0,w1,h1), R2(0,0,w2,h2);
      display_match(vec_matchings,vec_inliers, ok?&F:0, T1,T2, R1,R2, &in,&out);
      libs::WriteImage(argv[5], in);
      libs::WriteImage(argv[6], out);
    }
    if(argc!=7) // Show (mini) epipolar lines
    {
      Image<RGBColor> im;
      libs::convertImage(concat, &im);
      libs::DrawLine(int(w1*z),0, int(w1*z),int(concat.Height()), BLUE, &im);
      const double margin1 = sqrt((double)w1*w1+h1*h1) *.02;
      const double margin2 = sqrt((double)w2*w2+h2*h2) *.02;
      // Draw in red outliers
      for(size_t index = 0; index < vec_matchings.size(); ++index)
        if(std::find(vec_inliers.begin(),vec_inliers.end(),index) == vec_inliers.end())
          draw_points(vec_matchings[index], RED, &im, T1, T2);
      // Draw in green inliers and small epipolar lines
      for(it=vec_inliers.begin(); it!=vec_inliers.end(); ++it)
      {
        const Match& m = vec_matchings[*it];
        draw_points(m, GREEN, &im, T1, T2);
        draw_small_epi(F,     inv(m), T1, margin1, &im);
        draw_small_epi(F.t(),     m,  T2, margin2, &im);
      }
      libs::WriteImage(argv[argc-1], im);
    }
  }

  if(! ok)
  {
    std::cerr << "Failed to estimate a model" << std::endl;
    return 1;
  }

  return 0;
}

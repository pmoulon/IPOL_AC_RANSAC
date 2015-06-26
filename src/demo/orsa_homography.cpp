/**
 * @file orsa_homography.cpp
 * @brief Homographic image registration
 * @author Pascal Monasse, Pierre Moulon
 * 
 * Copyright (c) 2011-2012 Pascal Monasse, Pierre Moulon
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
#include <cmath>

#include "libImage/image_io.hpp"
#include "libImage/image_drawing.hpp"
#include "libImage/image_crop.hpp"

#include "libOrsa/homography_model.hpp"

#include "extras/libNumerics/numerics.h"
#include "extras/sift/library.h"

#include "demo/siftMatch.hpp"
#include "demo/warping.hpp"
#include "demo/cmdLine.h"

/// Number of random samples in ORSA
static const int ITER_ORSA=10000;

/// Display average/max error of inliers of homography H.
static void display_stats(const std::vector<Match>& vec_matchings,
                          const std::vector<int>& vec_inliers,
                          libNumerics::matrix<double>& H) {
  std::vector<int>::const_iterator it=vec_inliers.begin();
  double l2=0, linf=0;
  for(; it!=vec_inliers.end(); ++it) {
    const Match& m=vec_matchings[*it];
    double x1=m.x1, y1=m.y1;
    TransformH(H, x1, y1);
    double e = (m.x2-x1)*(m.x2-x1) + (m.y2-y1)*(m.y2-y1);
    l2 += e;
    if(linf < e)
      linf = e;
  }
  std::cout << "Average/max error: "
            << sqrt(l2/vec_inliers.size()) << "/"
            << sqrt(linf) <<std::endl;
}

/// ORSA homography estimation
bool ORSA(const std::vector<Match>& vec_matchings, int w1,int h1, int w2,int h2,
          double precision,
          libNumerics::matrix<double>& H, std::vector<int>& vec_inliers)
{
  const int n = static_cast<int>( vec_matchings.size() );
  if(n < 5)
  {
      std::cerr << "Error: ORSA needs 5 matches or more to proceed" <<std::endl;
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

  orsa::HomographyModel model(xA, w1, h1, xB, w2, h2, true);
  //model.setConvergenceCheck(true);

  if(model.orsa(vec_inliers, ITER_ORSA, &precision, &H, true)>0.0)
    return false;
  std::cout << "Before refinement: ";
  display_stats(vec_matchings, vec_inliers, H);
  if( model.ComputeModel(vec_inliers,&H) ) // Re-estimate with all inliers
  {
    std::cout << "After  refinement: ";
    display_stats(vec_matchings, vec_inliers, H);
  } else
    std::cerr << "Warning: error in refinement, result is suspect" <<std::endl;
  return true;
}

/// Output inlier and oulier matches in image files.
void display_match(const std::vector<Match>& vec_matchings,
                   std::vector<int>& vec_inliers,
                   const libNumerics::matrix<double>* H,
                   const libNumerics::matrix<double>& H1,
                   const libNumerics::matrix<double>& H2,
                   Image<RGBColor>& in, Image<RGBColor>& out)
{
  std::sort(vec_inliers.begin(), vec_inliers.end());

  // For outliers, show vector (yellow) from prediction to observation
  const RGBColor col=YELLOW;
  std::vector<int>::const_iterator it = vec_inliers.begin();
  matchingslist::const_iterator m = vec_matchings.begin();
  if(H) // Otherwise, no prediction
    for(int i=0; m != vec_matchings.end(); ++m, ++i) {
      if(it != vec_inliers.end() && i==*it)
        ++it;
      else { //Outlier
          double x1=m->x1, y1=m->y1, x2=m->x2, y2=m->y2;
          TransformH(H2 * *H, x1, y1);
          TransformH(H2, x2, y2);
          libs::DrawLine((int)x1,(int)y1,(int)x2,(int)y2, col,&out);
      }
    }

  // Show link for inliers (green) and outliers (red)
  it = vec_inliers.begin();
  m = vec_matchings.begin();
  for(int i=0; m != vec_matchings.end(); ++m, ++i)
  {
    Image<RGBColor>* im=&out;
    RGBColor col=RED;
    if(it != vec_inliers.end() && i==*it) {
      ++it;
      im=&in;
      col=GREEN;
    }
    double x1=m->x1, y1=m->y1, x2=m->x2, y2=m->y2;
    TransformH(H1, x1, y1);
    TransformH(H2, x2, y2);
    libs::DrawLine((int)x1,(int)y1,(int)x2, (int)y2, col, im);
  }
}

/// Return 3x3 translation matrix
libNumerics::matrix<double> translation(double dx, double dy) {
    libNumerics::matrix<double> T = libNumerics::matrix<double>::eye(3);
    T(0,2) = dx;
    T(1,2) = dy;
    return T;
}

/// Return 3x3 zoom-translation matrix
libNumerics::matrix<double> zoomtrans(double z, double dx, double dy) {
    libNumerics::matrix<double> T = libNumerics::matrix<double>::eye(3);
    T(0,0)=T(1,1) = z;
    T(0,2) = dx;
    T(1,2) = dy;
    return T;
}

int main(int argc, char **argv)
{
    CmdLine cmd;
    Geometry region;
    double precision=0;
    cmd.add( make_option('c',region,"cut") );
    cmd.add( make_option('p',precision, "prec") );
    cmd.add( make_switch('r', "read") );
    float fSiftRatio=0.6f;
    cmd.add( make_option('s',fSiftRatio, "sift") );
    try {
        cmd.process(argc, argv);
    } catch(const std::string& s) {
        std::cerr << s << std::endl;
        return 1;
    }
  if(argc!=5 && argc!=7 && argc!=8 && argc!=10) {
    std::cerr << "Usage: " << argv[0] << ' '
              << "[-c|--cut geom] "
              << "[-p|--prec precision] "
              << "[-s|--sift siftRatio] "
              << "[-r|--read] "
              << "imgInA imgInB "
              << "allMatches.txt orsaMatches.txt "
              << "[imgInliers imgOutliers [imgMosaic "
              << "[imgMosaicA imgMosaicB]]"
              << std::endl;
    std::cerr << " geom is the cut region: wxh+x+y" << std::endl
              << "  = rectangle of corners (x,y) and (x+w,y+h)" << std::endl;
    return 1;
  }

  // Init random seed
  time_t seed = time(0); // Replace by a fixed value to debug a reproducible run
  srand((unsigned int)seed);

  Image<RGBColor> image1, image2;
  if(! libs::ReadImage(argv[1], &image1))
    return 1;
  if(! libs::ReadImage(argv[2], &image2))
    return 1;

  bool bUseRegion = cmd.used('c');
  if(bUseRegion) { // Sanity check
    Geometry zone;
    zone.x0=zone.y0=0; zone.x1=int(image1.Width()); zone.y1=int(image1.Height());
    if(! (region.x0<region.x1 && region.y0<region.y1 &&
          zone.inside(region.x0,region.y0) &&
          zone.inside(region.x1-1,region.y1-1))) {
      std::cout << "Invalid cut region " << region << std::endl;
      return 1;
    }
  }

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
      SIFT(image1Gray, image2Gray, vec_matchings, fSiftRatio,
           bUseRegion? &region: 0);

  if(bUseRegion) {
      Image<RGBColor> image = image1;
      Crop(image, region.x0, region.y0,
           region.x1-region.x0, region.y1-region.y0, image1);
      libs::convertImage(image1, &image1Gray);
  }

  // Remove duplicates (frequent with SIFT)
  rm_duplicates(vec_matchings);

  // Save match files
  if(! cmd.used('r') && ! Match::saveMatch(argv[3], vec_matchings)) {
    std::cerr << "Failed saving matches into " <<argv[3] <<std::endl;
    return 1;
  }

  const int w1=int(image1Gray.Width()), h1=int(image1Gray.Height());
  const int w2=int(image2Gray.Width()), h2=int(image2Gray.Height());

  // Estimation of homography with ORSA
  libNumerics::matrix<double> H(3,3);
  std::vector<int> vec_inliers;
  bool ok = ORSA(vec_matchings, w1, h1, w2, h2, precision, H, vec_inliers);
  if(ok)
  {
    H /= H(2,2);
    std::cout << "H=" << H <<std::endl;
  }

  std::vector<Match> good_match;
  std::vector<int>::const_iterator it = vec_inliers.begin();
  for(; it != vec_inliers.end(); it++)
    good_match.push_back(vec_matchings[*it]);
  if(! Match::saveMatch(argv[4], good_match)) {
    std::cerr << "Failed saving matches into " <<argv[4] <<std::endl;
    return 1;
  }

  // Sift de-duplicated output display
  if(argc>6)
  {
    int w = std::max(w1,w2);
    float z = w/(float)(w1+w2); //Set width as max of two images
    Image<unsigned char> concat(w, int(z*std::max(h1,h2)), 255);
    libNumerics::matrix<double> T1=zoomtrans(z, 0,   (concat.Height()-h1*z)/2);
    libNumerics::matrix<double> T2=zoomtrans(z, w1*z,(concat.Height()-h2*z)/2);
    Warp(image1Gray, T1, image2Gray, T2, concat);
    Image<RGBColor> in;
    libs::convertImage(concat, &in);
    libs::DrawLine(int(w1*z),0, int(w1*z),int(concat.Height()), BLUE, &in);
    Image<RGBColor> out(in);
    display_match(vec_matchings, vec_inliers, ok? &H:0, T1, T2, in, out);
    libs::WriteImage(argv[5], in);
    libs::WriteImage(argv[6], out);
  }
  if(! ok)
  {
    std::cerr << "Failed to estimate a model" << std::endl;
    return 1;
  }

  // Mosaics
  if(argc>7)
  {
    cout << "-- Render Mosaic -- " << endl;

    Rect intersection;
    if(IntersectionBox(w1, h1, w2, h2, H, intersection) &&
      intersection.Width() > 0 && intersection.Height() > 0)
    {
      int xc=(intersection.left+intersection.right)/2;
      int yc=(intersection.top+intersection.bottom)/2;
      size_t wM = std::max(image1.Width(), image2.Width());
      size_t hM = std::max(image1.Height(), image2.Height());
      int xo=static_cast<int>(wM/2);
      int yo=static_cast<int>(hM/2);
      libNumerics::matrix<double> T = translation(xo-xc, yo-yc);

      Image<RGBColor> imageMosaic(wM, hM);
      imageMosaic.fill(WHITE);
      Warp(image1, T*H, image2, T, imageMosaic);
      libs::WriteImage(argv[7], imageMosaic);

      if (argc>9) // Registered images
      {
        Image<RGBColor> WarpingA(wM, hM, WHITE);
        Image<RGBColor> WarpingB(wM, hM, WHITE);

        cout << "-- Render Mosaic - Image A -- " << endl;
        Warp(image1, T*H, WarpingA);

        cout << "-- Render Mosaic - Image B -- " << endl;
        H = libNumerics::matrix<double>::eye(3);
        Warp(image2, T, WarpingB);

        libs::WriteImage(argv[8], WarpingA);
        libs::WriteImage(argv[9], WarpingB);
      }
    }
  } 
  return 0;
}

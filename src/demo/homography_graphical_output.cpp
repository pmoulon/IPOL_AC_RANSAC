/**
 * @file homography_graphical_output.cpp
 * @brief Graphical output to show homography estimation
 * @author Pascal Monasse, Pierre Moulon
 * 
 * Copyright (c) 2011-2016 Pascal Monasse, Pierre Moulon
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

#include "homography_graphical_output.hpp"
#include <algorithm>
#include "libImage/image_drawing.hpp"
#include "libImage/image_io.hpp"
#include "warping.hpp"

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
  std::vector<Match>::const_iterator m = vec_matchings.begin();
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

void homography_matches_output(const Image<unsigned char>& image1,
                               const Image<unsigned char>& image2,
                               const std::vector<Match>& vec_all,
                               std::vector<int> vec_in,
                               const libNumerics::matrix<double>* H,
                               const char* fileIn, const char* fileOut) {
    cout << "-- Render Mosaic -- " << endl;
    const int w1=image1.Width(), h1=image1.Height();
    const int w2=image2.Width(), h2=image2.Height();
    int w = std::max(w1,w2);    // Set width as max of two images
    float z = w/(float)(w1+w2); // Keep aspect ratio
    Image<unsigned char> concat(w, int(z*std::max(h1,h2)), 255);

    libNumerics::matrix<double> T1=zoomtrans(z, 0,   (concat.Height()-h1*z)/2);
    libNumerics::matrix<double> T2=zoomtrans(z, w1*z,(concat.Height()-h2*z)/2);
    Warp(image1, T1, image2, T2, concat);

    Image<RGBColor> in;
    libs::convertImage(concat, &in);
    libs::DrawLine(int(w1*z),0, int(w1*z),int(concat.Height()), BLUE, &in);
    Image<RGBColor> out(in);
    display_match(vec_all, vec_in, H, T1, T2, in, out);

    libs::WriteImage(fileIn,  in);
    libs::WriteImage(fileOut, out);
}

void homography_registration_output(const Image<RGBColor>& image1,
                                    const Image<RGBColor>& image2,
                                    const libNumerics::matrix<double>& H,
                                    const Rect& rect,
                                    const char* fileMosaic,
                                    const char* fileReg1,
                                    const char* fileReg2) {
    int xc=(rect.left+rect.right)/2;
    int yc=(rect.top+rect.bottom)/2;
    size_t wM = std::max(image1.Width(), image2.Width());
    size_t hM = std::max(image1.Height(), image2.Height());
    int xo=static_cast<int>(wM/2);
    int yo=static_cast<int>(hM/2);
    libNumerics::matrix<double> T = translation(xo-xc, yo-yc);

    if(fileMosaic)
    {
        Image<RGBColor> imageMosaic(wM, hM);
        imageMosaic.fill(WHITE);
        Warp(image1, T*H, image2, T, imageMosaic);
        libs::WriteImage(fileMosaic, imageMosaic);
    }
    if(fileReg1 || fileReg2)
    {
        Image<RGBColor> reg1(wM, hM, WHITE);
        Image<RGBColor> reg2(wM, hM, WHITE);

        cout << "-- Render Mosaic - Image 1 -- " << endl;
        Warp(image1, T*H, reg1);

        cout << "-- Render Mosaic - Image 2 -- " << endl;
        Warp(image2, T, reg2);

        libs::WriteImage(fileReg1, reg1);
        libs::WriteImage(fileReg2, reg2);
    }
}

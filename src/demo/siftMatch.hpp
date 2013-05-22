/**
 * @file siftMatch.hpp
 * @brief SIFT extraction and matching
 * @author Lionel Moisan, Pascal Monasse, Pierre Moulon
 * 
 * Copyright (c) 2011 Lionel Moisan, Pascal Monasse, Pierre Moulon
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

#ifndef SIFT_MATCH_H
#define SIFT_MATCH_H

#include "extras/sift/demo_lib_sift.h"
#include <istream>
#include <ostream>

/// Rectangle in image
class Geometry {
public:
    int x0, x1, y0, y1;
    bool inside(double x, double y) const {
        int ix=static_cast<int>(x), iy=static_cast<int>(y);
        return (x0<=ix && ix<x1 && y0<=iy && iy<y1);
    }
    bool inside(keypoint k) const {
        return inside(k.x, k.y);
    }
};

/// Output @geo.
std::ostream& operator<<(std::ostream& str, Geometry& geo) {
    return str <<geo.x1-geo.x0 <<'x'<<geo.y1-geo.y0 <<'+'<<geo.x0 <<'+'<<geo.y0;
}

/// Input @geo. Format: wxh+x0+y0, eg 100x100+0+0
std::istream& operator>>(std::istream& str, Geometry& geo) {
    char c;
    str >> geo.x1 >> c;
    if(str.fail() || c!='x') return str;
    str >> geo.y1 >> c;
    if(str.fail() || c!='+') return str;
    str >> geo.x0 >> c;
    if(str.fail() || c!='+') return str;
    str >> geo.y0;
    geo.x1 += geo.x0;
    geo.y1 += geo.y0;
    return str;
}

/// Functor to test inclusion in region
class OutRegion {
public:
    OutRegion(Geometry* r): rect(r) {}
    bool operator()(const keypoint& k) const {
        return !rect->inside(k);
    }
private:
    Geometry* rect;
};

/// SIFT matches
static void SIFT(const Image<unsigned char> &im1,
                 const Image<unsigned char> &im2,
                 std::vector<Match>& vec_matchings,
                 float fMatchRatio=0.6f, Geometry* rect=0) {
    //Convert images to float
    Image<float> If1, If2;
    libs::convertImage(im1, &If1);
    libs::convertImage(im2, &If2);

    siftPar param;
    default_sift_parameters(param);
    param.MatchRatio = fMatchRatio;
    param.DoubleImSize=0;

    // Keypoints in 1st image
    keypointslist keyp1;
    compute_sift_keypoints(If1.data(), keyp1, int(If1.Width()), int(If1.Height()), param);
    std::cout<< "sift:: 1st image: " << keyp1.size() << " keypoints";
    if(rect) { // Remove points outside region
        keypointslist::iterator
            oldEnd=keyp1.end(),
            newEnd=std::remove_if(keyp1.begin(), oldEnd, OutRegion(rect));
        std::cout << " (but remove " << std::distance(newEnd,oldEnd)
                  << " outside " << *rect << ')';
        keyp1.erase(newEnd,oldEnd);
        // Translate points
        for(keypointslist::iterator it=keyp1.begin(); it!=keyp1.end(); ++it) {
            it->x -= rect->x0;
            it->y -= rect->y0;
        }
    }
    std::cout << std::endl;

    // Keypoints in 2nd image
    keypointslist keyp2;
    compute_sift_keypoints(If2.data(), keyp2, int(If2.Width()), int(If2.Height()), param);
    std::cout<< "sift:: 2nd image: " << keyp2.size() << " keypoints"<<std::endl;

    // Find putatives matches
    compute_sift_matches(keyp1, keyp2, vec_matchings, param);
    std::cout << "sift:: matches: " << vec_matchings.size() <<std::endl;
}

/// Remove multiple "same position" matches
static void rm_duplicates(std::vector<Match>& m) {
  std::sort(m.begin(), m.end());
  std::vector<Match>::iterator end = std::unique(m.begin(), m.end());
  if(end != m.end()) {
    std::cout << "Remove " << std::distance(end, m.end())
      << "/" << m.size() << " duplicate matches, "
      << "keeping " << std::distance(m.begin(), end) <<std::endl;
    m.erase(end, m.end());
  }
}

#endif

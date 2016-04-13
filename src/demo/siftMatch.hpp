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

#include "third_party/sift_anatomy/lib_sift_anatomy.h"
#include "third_party/sift_anatomy/lib_matching.h"
#include "libOrsa/match.hpp"

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
        return inside(k.y, k.x); // Weird SiftAnatomy's coordinate system...
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

/// This should be in SiftAnatomy.
static struct sift_keypoints* sift_anatomy(const float* x, int w, int h,
                                           const struct sift_parameters* p)
{
    struct sift_keypoints *kk[6];
    for(int i = 0; i < 6; i++)
        kk[i] = sift_malloc_keypoints();
    struct sift_scalespace *ss[4];

    struct sift_keypoints* k = sift_anatomy(x, w, h, p, ss, kk);

    for(int i=0; i<6; i++)
        sift_free_keypoints(kk[i]);
    for(int i=0; i<4; i++)
        sift_free_scalespace(ss[i]);

    return k;
}


/// SIFT matches
static void SIFT(const Image<unsigned char> &im1,
                 const Image<unsigned char> &im2,
                 std::vector<Match>& vec_matchings,
                 float fMatchRatio=0.6f, Geometry* rect=0) {
    //Convert images to float
    Image<float> If1, If2;
    libs::convertImage(im1, &If1);
    libs::convertImage(im2, &If2);

    for(size_t y=0; y<If1.Height(); y++)
        for(size_t x=0; x<If1.Width(); x++)
            If1(y,x)/=256.0f;
    for(size_t y=0; y<If2.Height(); y++)
        for(size_t x=0; x<If2.Width(); x++)
            If2(y,x)/=256.0f;

    sift_parameters* param = sift_assign_default_parameters();
    // Start from octave 0
    param->delta_min = 1.0;
    param->sigma_min = 1.6;

    sift_keypoints* keyp1 = sift_anatomy(If1.data(), If1.Width(), If1.Height(),
                                         param);
    std::cout<< "sift:: 1st image: " << keyp1->size << " keypoints";
    if(rect) { // Remove points outside region
        int j=0;
        for(int i=0; i<keyp1->size; i++)
            if(rect->inside(*keyp1->list[i]))
                keyp1->list[j++] = keyp1->list[i];
            else
                sift_free_keypoint(keyp1->list[i]);
        if(j<keyp1->size)
            std::cout <<" (remove "<< keyp1->size-j <<" outside "<< *rect <<')';
        keyp1->size = j;
        // Translate points
        for(int i=0; i<keyp1->size; i++) {
            keyp1->list[i]->y -= rect->x0; // SiftAnatomy's weird coordinates
            keyp1->list[i]->x -= rect->y0;
        }
    }
    std::cout << std::endl;

    sift_keypoints* keyp2 = sift_anatomy(If2.data(), If2.Width(), If2.Height(),
                                         param);
    std::cout<< "sift:: 2nd image: " << keyp2->size << " keypoints"<<std::endl;

    // Find putatives matches
    sift_keypoints* m1 = sift_malloc_keypoints();
    sift_keypoints* m2 = sift_malloc_keypoints();
    sift_keypoints* unused = sift_malloc_keypoints();
    matching(keyp1, keyp2, m1, m2, unused, fMatchRatio, 1);
    vec_matchings.clear();
    for(int i=0; i<m1->size; i++) {
        // Adapt to SiftAnatomy's weird coordinate system...
        Match m(m1->list[i]->y, m1->list[i]->x, m2->list[i]->y, m2->list[i]->x);
        vec_matchings.push_back(m);
    }
    std::cout << "sift:: matches: " << vec_matchings.size() <<std::endl;
    free(param);
    sift_free_keypoints(keyp1);
    sift_free_keypoints(keyp2);
    sift_free_keypoints(m1);
    sift_free_keypoints(m2);
    sift_free_keypoints(unused);    
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

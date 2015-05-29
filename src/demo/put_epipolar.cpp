/**
 * @file put_epipolar.cpp
 * @brief Write a point or epipolar line in transparent image
 * @author Pascal Monasse
 *
 * Copyright (c) 2013 Pascal Monasse
 * All rights reserved.
 *
 * This program is free software: you can use, modify and/or
 * redistribute it under the terms of the simplified BSD License. You
 * should have received a copy of this license along this program. If
 * not, see <http://www.opensource.org/licenses/bsd-license.html>.
 */

#include "Rect.hpp"
#include "cmdLine.h"
#include "libImage/image_io.hpp"
#include "libImage/image_drawing.hpp"
#include "extras/libNumerics/matrix.h"
#include <sstream>
#include <fstream>

static const RGBA REDA(RGBColor(255,0,0));

int main(int argc, char** argv) {
    std::string fileF;
    CmdLine cmd;
    cmd.add( make_switch('t',"transpose") );
    cmd.add( make_option('f',fileF,"fmatrix") );
    try {
        cmd.process(argc,argv);
    } catch(const std::string& c) {
        std::cerr << c << std::endl;
        return 1;
    }
    if(argc!=4 && argc!=5) {
        std::cerr << "Usage: " << argv[0]
                  << " [-f|--fmatrix fileF.txt [-t|--transpose]]"
                  << " x y img [imgOut]" << std::endl;
        return 1;
    }

    int x,y;
    if(! (std::stringstream(argv[1])>>x).eof()) {
        std::cerr << "Unable to read " << argv[1] << " as integer" << std::endl;
        return 1;
    }
    if(! (std::stringstream(argv[2])>>y).eof()) {
        std::cerr << "Unable to read " << argv[2] << " as integer" << std::endl;
        return 1;
    }

    Image<RGBA> in;
    if(! libs::ReadImage(argv[3], &in)) {
        std::cerr << "Error reading image file " << argv[3] << std::endl;
        return 1;
    }

    if(! fileF.empty()) {
        libNumerics::matrix<float> F(3,3);
        std::fstream file(fileF.c_str());
        if( (file>>F).fail() ) {
            std::cerr << "Error reading file " << fileF << std::endl;
            return 1;
        }
        if(cmd.used('t'))
            F = F.t();
        libNumerics::vector<float> l(3);
        l(0) = static_cast<float>(x);
        l(1) = static_cast<float>(y);
        l(2) = 1.0f;
        l = F*l;
        Rect R(0,0,in.Width(),in.Height());
        if( R.intersect(l(0),l(1),l(2)) )
            libs::DrawLine((int)R.left,(int)R.top, (int)R.right,(int)R.bottom,
                           REDA, &in);
    } else {
        if(cmd.used('t')) {
            std::cerr << "Error: option -t must be used with -f" << std::endl;
            return 1;
        }
        libs::DrawCircle(x,y, 2, REDA, &in);
    }

    if(! libs::WriteImage(argv[argc-1], in)) {
        std::cerr << "Error writing file " << argv[2] << std::endl;
        return 1;
    }
    return 0;
}

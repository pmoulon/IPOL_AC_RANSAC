ORSA - A RANSAC variant with a contrario elimination of outliers

Lionel Moisan lionel.moisan@parisdescartes.fr MAP5, Universite Paris Descartes
Pierre Moulon pmoulon@gmail.com IMAGINE/LIGM, Universite Paris Est
Pascal Monasse monasse@imagine.enpc.fr IMAGINE/LIGM, Universite Paris Est

The generic framework is applied to two model estimations (IPOL publications):
- homography matrix  (https://doi.org/10.5201/ipol.2012.mmm-oh)
- fundamental matrix (https://doi.org/10.5201/ipol.2016.147)

Future releases and updates:
http://imagine.enpc.fr/~moulonp/AC_Ransac.html

- Description:
This provides a variant, called ORSA, of the well known RANSAC method for
model parameter estimation. It is based on an a contrario criterion of
inlier/outlier discrimination, is parameter free and relies on an optimized
random sampling procedure.

- Licensing: See LICENSE.txt file
- Build, usage: See BUILD.txt file.

- Reviewed files in IPOL:
src/libOrsa/libNumerics/cubicRoots.h
src/libOrsa/sampling.{hpp,cpp}
src/libOrsa/model_estimator.{hpp,cpp}
src/libOrsa/homography_model.{hpp,cpp}
src/libOrsa/fundamental_model.{hpp.cpp}
src/libOrsa/orsa_homography.{hpp,cpp}
src/libOrsa/orsa_fundamental.{hpp,cpp}
src/demo/demo_orsa_homography.cpp
src/demo/demo_orsa_fundamental.cpp
src/demo/homography_graphical_output.{hpp,cpp}
src/demo/fundamental_graphical_output.{hpp,cpp}
src/demo/put_epipolar.cpp
src/demo/Rect.{hpp,cpp}

- The SIFT algorithm used by the demo is based on
  sift_anatomy_20141201
from the article
  Ives Rey Otero, and Mauricio Delbracio, Anatomy of the SIFT Method,
  Image Processing On Line, 4 (2014), pp. 370–396.
  https://doi.org/10.5201/ipol.2014.82

- Reusing the libOrsa library:
This library itself is self-contained and has no external dependencies. To reuse
it in your own program, just copy the folder along with its sub-folder
libNumerics in your program. You may need to remove the UNIT_TEST instructions
in the CMakeLists.txt files if you do not want to take along the third-party
CppUnitLite. The simplest usage is through the function
- Homography:
  bool orsa_homography(const std::vector<Match>& vec_matchings,
                       int w1,int h1, int w2,int h2,
                       double precision, int nbIter,
                       libNumerics::matrix<double>& H,
                       std::vector<int>& vec_inliers);
A lower-level function, just applying ORSA with no refinement:
  HomographyModel model(xA, w1, h1, xB, w2, h2);
  Orsa orsa(&model, M_PI/w1/h1, M_PI/w2/h2);
  orsa.run(vec_inliers, 10000, 0, &H);
- Fundamental:
  bool orsa_fundamental(const std::vector<Match>& vec_matchings,
                        int w1,int h1, int w2,int h2,
                        double precision, int nbIter,
                        libNumerics::matrix<double>& F,
                        std::vector<int>& vec_inliers);
A lower-level function, just applying ORSA with no refinement:
  FundamentalModel model(xA, w1, h1, xB, w2, h2);
  Orsa orsa(&model, 2*hypot(w1,h1)/(w1*h1), 2*hypot(w2,h2)/(w2*h2));
  orsa.run(vec_inliers, 10000, 0, &F);

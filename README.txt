Fundamental Matrix of a Stereo Pair, with A Contrario Elimination of Outliers
Lionel Moisan lionel.moisan@parisdescartes.fr MAP5, Universite Paris Descartes
Pierre Moulon pmoulon@gmail.com IMAGINE/LIGM, Universite Paris Est
Pascal Monasse monasse@imagine.enpc.fr IMAGINE/LIGM, Universite Paris Est
Version 1.0rc1 released on 2015/07/21
Future releases and updates:
http://imagine.enpc.fr/~moulonp/AC_Ransac.html

- Description:
This provides a variant, called ORSA, of the well known RANSAC method for
model parameter estimation. It is based on an a contrario criterion of
inlier/outlier discrimination, is parameter free and relies on an optimized
random sampling procedure.
The generic approach is applied to the robust estimation of the fundamental
matrix relating two images. For the general ORSA principle, see the companion
article:
http://dx.doi.org/10.5201/ipol.2012.mmm-oh

- Licensing: See LICENSE.txt file
- Build, usage: See BUILD.txt file.

- Reviewed files in IPOL:
src/libOrsa/libNumerics/cubicRoots.h
src/libOrsa/fundamental_model.hpp
src/libOrsa/fundamental_model.cpp
src/libOrsa/orsa_fundamental.hpp
src/libOrsa/orsa_fundamental.cpp
src/demo/demo_orsa_fundamental.cpp
src/demo/fundamental_graphical_output.hpp
src/demo/fundamental_graphical_output.cpp
src/demo/put_epipolar.cpp
src/demo/Rect.hpp
src/demo/Rect.cpp

- The SIFT algorithm used by the demo is based on
  sift_anatomy_20141201
from the article
  Ives Rey Otero, and Mauricio Delbracio, Anatomy of the SIFT Method,
  Image Processing On Line, 4 (2014), pp. 370–396.
  http://dx.doi.org/10.5201/ipol.2014.82

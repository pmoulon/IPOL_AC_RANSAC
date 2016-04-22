Name: OrsaHomography
Version: 20160422

Long name:
Automatic homographic registration of a pair of images, with a contrario elimination of outliers

Brief description:
An a contrario criterion to estimate the rigid transform (homography)
registering two images using a parameter free variant of RANSAC, ORSA.

IPOL publication:
http://www.ipol.im/pub/algo/mmm_orsa_homography/

Future releases and updates:
http://imagine.enpc.fr/~moulonp/AC_Ransac.html

Authors:
Lionel Moisan lionel.moisan[AT]parisdescartes.fr MAP5, Universite Paris Descartes
Pierre Moulon pmo[AT]mikrosimage.eu IMAGINE/LIGM, Universite Paris Est & Mikros Image
Pascal Monasse monasse[AT]imagine.enpc.fr IMAGINE/LIGM, Universite Paris Est

Keywords: a contrario method, image matching, homography, robust estimation

Comments:
This provides a variant, called ORSA, of the well known RANSAC method for
model parameter estimation. It is based on an a contrario criterion of
inlier/outlier discrimination, is parameter free and relies on an optimized
random sampling procedure.
The generic approach is explained and applied to the robust estimation of a
homography registering two images.

Licensing: See LICENSE.txt file

Build, usage: See BUILD.txt file.

- The SIFT algorithm used by the demo is based on
  sift_anatomy_20141201
from the article
  Ives Rey Otero, and Mauricio Delbracio, Anatomy of the SIFT Method,
  Image Processing On Line, 4 (2014), pp. 370–396.
  http://dx.doi.org/10.5201/ipol.2014.82

- Reusing the libOrsa library:
This library itself is self-contained and has no external dependencies. To reuse
it in your own program, just copy the folder along with its sub-folder
libNumerics in your program. You may need to remove the UNIT_TEST instructions
in the CMakeLists.txt files if you do not want to take along the third-party
CppUnitLite. The simplest usage is through the function
  bool orsa_homography(const std::vector<Match>& vec_matchings,
                       int w1,int h1, int w2,int h2,
                       double precision, int nbIter,
                       libNumerics::matrix<double>& H,
                       std::vector<int>& vec_inliers);
A lower-level function, just applying ORSA with no refinement:
  HomographyModel model(xA, w1, h1, xB, w2, h2);
  model.orsa(vec_inliers, 10000, 0, &H);

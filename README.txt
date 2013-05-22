Name: OrsaHomography
Version: 20120515

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

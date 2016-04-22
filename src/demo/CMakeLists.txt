ADD_SUBDIRECTORY(libImage)

FIND_PACKAGE(OpenMP) 
if (OPENMP_FOUND) 
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif (OPENMP_FOUND)

# Executable demo to estimate a homography matrix from SIFT matches and
# display the corresponding mosaic.
ADD_EXECUTABLE(demo_orsa_homography
               demo_orsa_homography.cpp
               homography_graphical_output.cpp homography_graphical_output.hpp
               Rect.cpp Rect.hpp
               siftMatch.hpp
               warping.cpp warping.hpp
               cmdLine.h)
TARGET_LINK_LIBRARIES(demo_orsa_homography Orsa image sift)

IF(UNIX)
    SET_TARGET_PROPERTIES(demo_orsa_homography PROPERTIES
                          COMPILE_FLAGS "-std=c++98 -Wall -Wextra -Werror")
ENDIF(UNIX)

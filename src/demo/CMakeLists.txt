ADD_SUBDIRECTORY(libImage)

FIND_PACKAGE(OpenMP) 
if (OPENMP_FOUND) 
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif (OPENMP_FOUND)

# Executable demo to estimate the fundamental matrix from SIFT matches
ADD_EXECUTABLE(demo_orsa_fundamental
               demo_orsa_fundamental.cpp
               fundamental_graphical_output.cpp
               Rect.cpp)
TARGET_LINK_LIBRARIES(demo_orsa_fundamental Orsa image SiftAnatomy)

# Utility to overlay points and lines on images (used by IPOL demo)
ADD_EXECUTABLE(put_epipolar put_epipolar.cpp Rect.cpp)
TARGET_LINK_LIBRARIES(put_epipolar image)

IF(UNIX)
    SET_TARGET_PROPERTIES(demo_orsa_fundamental PROPERTIES
                          COMPILE_FLAGS "-std=c++98 -Wall -Wextra -Werror")
    SET_TARGET_PROPERTIES(put_epipolar PROPERTIES
                          COMPILE_FLAGS "-std=c++98 -Wall -Wextra -Werror")
ENDIF(UNIX)

# Show small epipolar lines
ADD_EXECUTABLE(show_epi show_epi.cpp fundamental_graphical_output.cpp Rect.cpp)
TARGET_LINK_LIBRARIES(show_epi Orsa image SiftAnatomy)

ADD_LIBRARY(SiftAnatomy
            lib_description.h lib_description.cpp
            lib_scalespace.h lib_scalespace.cpp
            lib_sift.h lib_sift.cpp
            lib_keypoint.h lib_keypoint.cpp
            lib_discrete.h lib_discrete.cpp
            lib_matching.h lib_matching.cpp
            lib_util.h lib_util.cpp
            lib_sift_anatomy.h lib_sift_anatomy.cpp)

IF(UNIX)
    SET_TARGET_PROPERTIES(SiftAnatomy PROPERTIES
                          COMPILE_FLAGS "-std=c++98 -Wall -Wextra")
ENDIF(UNIX)

ADD_LIBRARY(SiftAnatomy
            lib_description.h lib_description.c
            lib_scalespace.h lib_scalespace.c
            lib_sift.h lib_sift.c
            lib_keypoint.h lib_keypoint.c
            lib_discrete.h lib_discrete.c
            lib_matching.h lib_matching.c
            lib_util.h lib_util.c
            lib_sift_anatomy.h lib_sift_anatomy.c)

IF(UNIX)
    SET_TARGET_PROPERTIES(SiftAnatomy PROPERTIES
                          COMPILE_FLAGS "-std=c99 -Wall")
ENDIF(UNIX)

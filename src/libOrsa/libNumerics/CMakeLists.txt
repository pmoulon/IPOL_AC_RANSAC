project(libNumerics)

set(SRC 
    cubicRoots.h
    homography.h
    matrix.cpp matrix.h
    numerics.cpp numerics.h
    svd.cpp
    vector.cpp)

add_library(Numerics ${SRC})

if(UNIX)
    set_target_properties(Numerics PROPERTIES
                          COMPILE_FLAGS "-std=c++98 -Wall -Wextra")
endif()

#----------
# UNIT TEST
#----------
UNIT_TEST(cubicRoots Numerics)
UNIT_TEST(svd Numerics)

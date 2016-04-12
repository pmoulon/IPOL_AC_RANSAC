PROJECT(libNumerics)

SET(SRC 
        cubicRoots.h
        ccmath_svd.cpp
        matrix.cpp matrix.h
        numerics.cpp numerics.h
        vector.cpp)

ADD_LIBRARY(Numerics ${SRC})

IF(UNIX)
    SET_TARGET_PROPERTIES(Numerics PROPERTIES
                          COMPILE_FLAGS "-std=c++98 -Wall -Wextra -Werror")
ENDIF(UNIX)

#----------
# UNIT TEST
#----------
UNIT_TEST(cubicRoots Numerics)
UNIT_TEST(ccmath_svd Numerics)

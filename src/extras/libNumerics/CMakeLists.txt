PROJECT(libNumerics)

SET(SRC 
        ccmath_svd.cpp
        matrix.cpp matrix.h
        numerics.cpp numerics.h
        vector.cpp)

ADD_LIBRARY(Numerics ${SRC})

IF(UNIX)
    SET_TARGET_PROPERTIES(Numerics PROPERTIES
                          COMPILE_FLAGS "-std=c++98 -Wall -Wextra -Werror")
ENDIF(UNIX)

INSTALL(FILES CMakeLists.txt ${SRC} DESTINATION src/extras/libNumerics)

#----------
# UNIT TEST
#----------
UNIT_TEST(ccmath_svd Numerics)
INSTALL(FILES ccmath_svd_test.cpp DESTINATION src/extras/libNumerics)

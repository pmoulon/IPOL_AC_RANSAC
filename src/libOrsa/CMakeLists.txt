SET(SRC conditioning.hpp conditioning.cpp
        match.hpp
        model_estimator.hpp model_estimator.cpp
        homography_model.hpp homography_model.cpp
        fundamental_model.hpp fundamental_model.cpp
        sampling.hpp sampling.cpp
        ransac.hpp ransac.cpp
        orsa.hpp orsa.cpp
        eval_model.hpp eval_model.cpp)

ADD_LIBRARY(Orsa ${SRC})

IF(UNIX)
    SET_TARGET_PROPERTIES(Orsa PROPERTIES
                          COMPILE_FLAGS "-std=c++98 -Wall -Wextra")
ENDIF(UNIX)

IF(MSVC)
    ADD_DEFINITIONS(-D_CRT_SECURE_NO_WARNINGS)
ENDIF(MSVC)

ADD_SUBDIRECTORY(libNumerics)

TARGET_LINK_LIBRARIES(Orsa Numerics)

#----------
# UNIT TEST
#----------

UNIT_TEST(fundamental_model Orsa)
UNIT_TEST(fundamental_orsa Orsa)
UNIT_TEST(homography_model Orsa)
UNIT_TEST(homography_orsa Orsa)

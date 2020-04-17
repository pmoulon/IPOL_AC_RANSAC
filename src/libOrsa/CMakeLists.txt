SET(SRC conditioning.hpp conditioning.cpp
        match.hpp
        model_estimator.hpp model_estimator.cpp
        homography_model.hpp homography_model.cpp
        fundamental_model.hpp fundamental_model.cpp
        orsa.hpp orsa.cpp
        orsa_homography.hpp orsa_homography.cpp
        orsa_fundamental.hpp orsa_fundamental.cpp)

ADD_LIBRARY(Orsa ${SRC})

IF(UNIX)
    SET_TARGET_PROPERTIES(Orsa PROPERTIES
                          COMPILE_FLAGS "-std=c++98 -Wall -Wextra -Werror")
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

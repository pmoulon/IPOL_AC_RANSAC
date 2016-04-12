SET(SRC conditioning.hpp conditioning.cpp
	    orsa_model.hpp orsa_model.cpp
        homography_model.hpp homography_model.cpp
        orsa_homography.hpp orsa_homography.cpp
        match.hpp)

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

# UNIT TEST -> HOMOGRAPHY

UNIT_TEST(homography_model "Orsa")
UNIT_TEST(homography_orsa "Orsa")

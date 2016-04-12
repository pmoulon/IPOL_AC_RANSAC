SET(SRC conditioning.hpp conditioning.cpp
	    orsa_model.hpp orsa_model.cpp
        fundamental_model.hpp fundamental_model.cpp
        orsa_fundamental.hpp orsa_fundamental.cpp
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

UNIT_TEST(fundamental_model "Orsa")
UNIT_TEST(fundamental_orsa "Orsa")

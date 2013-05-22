SET(SRC conditioning.hpp conditioning.cpp
	    orsa_model.hpp orsa_model.cpp
#        fundamental_model.hpp fundamental_model.cpp
        homography_model.hpp homography_model.cpp)

ADD_LIBRARY(Orsa ${SRC})

IF(UNIX)
    SET_TARGET_PROPERTIES(Orsa PROPERTIES
                          COMPILE_FLAGS "-std=c++98 -Wall -Wextra -Werror")
ENDIF(UNIX)

TARGET_LINK_LIBRARIES(Orsa Numerics)

INSTALL(FILES CMakeLists.txt ${SRC} DESTINATION src/libOrsa)

# ADD_SUBDIRECTORY(numeric)

#----------
# UNIT TEST
#----------

# UNIT TEST -> HOMOGRAPHY

UNIT_TEST(homography_model "Orsa")
UNIT_TEST(homography_orsa "Orsa")

INSTALL(FILES homography_model_test.cpp
              homography_orsa_test.cpp
        DESTINATION src/libOrsa)

PROJECT(sift)

SET(SRC library.cpp library.h
        domain.cpp domain.h
        splines.cpp splines.h
        filter.cpp filter.h
        demo_lib_sift.cpp demo_lib_sift.h
        )

INCLUDE_DIRECTORIES(..)

ADD_LIBRARY(sift ${SRC})
IF(UNIX)
    SET_TARGET_PROPERTIES(sift PROPERTIES
                          COMPILE_FLAGS "-std=c++98 -Wall -Wextra -Werror")
ENDIF(UNIX)

INSTALL(FILES CMakeLists.txt ${SRC} DESTINATION src/extras/sift)

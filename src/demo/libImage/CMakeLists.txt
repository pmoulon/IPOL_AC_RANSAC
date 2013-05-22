
ADD_DEFINITIONS(-DTHIS_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")

SET(SRC image.hpp
        sample.hpp
        pixelTypes.hpp
        image_converter.hpp
        image_crop.hpp
        image_io.hpp
        image_io.cpp
        image_drawing.hpp)

ADD_LIBRARY(image ${SRC})

TARGET_LINK_LIBRARIES(image png jpeg)

FIND_PACKAGE(PNG)
IF(NOT PNG_FOUND)
    SET(PNG_INCLUDE_DIR ../../third_party/png)
ENDIF(NOT PNG_FOUND)

FIND_PACKAGE(JPEG)
IF(NOT JPEG_FOUND)
    ADD_DEFINITIONS(-DHAVE_PROTOTYPES)
    SET(JPEG_INCLUDE_DIR ../../third_party/jpeg)
ENDIF(NOT JPEG_FOUND)

ADD_DEFINITIONS(${PNG_DEFINITIONS})
INCLUDE_DIRECTORIES(${PNG_INCLUDE_DIR})
INCLUDE_DIRECTORIES(${JPEG_INCLUDE_DIR})

IF(UNIX)
    SET_TARGET_PROPERTIES(image PROPERTIES
                          COMPILE_FLAGS "-std=c++98 -Wall -Wextra -Werror")
ENDIF(UNIX)

IF(WIN32)
    ADD_DEFINITIONS(-D_CRT_SECURE_NO_WARNINGS)
ENDIF(WIN32)

UNIT_TEST(image "image")
UNIT_TEST(image_drawing "image")
UNIT_TEST(image_io "image")

INSTALL(FILES CMakeLists.txt ${SRC} DESTINATION src/demo/libImage)
INSTALL(FILES image_test.cpp image_drawing_test.cpp image_io_test.cpp
        DESTINATION src/demo/libImage)
SET(IMG image_test/two_pixels_color.jpg image_test/two_pixels_monochrome.jpg
        image_test/two_pixels.ppm image_test/two_pixels_color.png
        image_test/two_pixels_monochrome.png image_test/two_pixels_gray.pgm
        image_test/two_pixels.pgm)
INSTALL(FILES ${IMG} DESTINATION src/demo/libImage/image_test)

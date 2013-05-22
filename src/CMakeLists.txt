CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

PROJECT(mmm_orsa C CXX)

ENABLE_TESTING()

INCLUDE_DIRECTORIES(
  ./demo/
  ./third_party/
  ./
)

IF (WIN32)
  INCLUDE_DIRECTORIES(
    third_party/jpeg
    third_party/zlib
    third_party/png)
ENDIF (WIN32)

MACRO (UNIT_TEST NAME EXTRA_LIBS)
  ADD_EXECUTABLE(${NAME}_test ${NAME}_test.cpp)
  TARGET_LINK_LIBRARIES(${NAME}_test
                        ${EXTRA_LIBS} # Extra libs MUST be first.
                        CppUnitLite)
  ADD_TEST(${NAME}_test ${NAME}_test)
ENDMACRO (UNIT_TEST)

# Third parties libraries
ADD_SUBDIRECTORY(third_party)

# Orsa Core library
ADD_SUBDIRECTORY(libOrsa)

# External libraries required for the demo
ADD_SUBDIRECTORY(extras)

# Demo directory
ADD_SUBDIRECTORY(demo)

INSTALL(FILES CMakeLists.txt DESTINATION src)
INSTALL(FILES testing/testing.h DESTINATION src/testing)

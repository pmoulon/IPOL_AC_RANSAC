project(CppUnitLite)
FILE(GLOB CppUnitLite_SRC "*.h" "*.cpp")
ADD_LIBRARY(CppUnitLite ${CppUnitLite_SRC})
SET_TARGET_PROPERTIES(CppUnitLite PROPERTIES LINKER_LANGUAGE CXX)

IF(WIN32)
    ADD_DEFINITIONS(-D_CRT_SECURE_NO_WARNINGS)
ENDIF(WIN32)

INSTALL(FILES CMakeLists.txt ${CppUnitLite_SRC} DESTINATION src/third_party/CppUnitLite)

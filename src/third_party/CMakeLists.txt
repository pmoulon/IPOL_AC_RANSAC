#---------------------
#Third party libraries
#---------------------

ADD_SUBDIRECTORY(CppUnitLite)

ADD_SUBDIRECTORY(zlib)
ADD_SUBDIRECTORY(png)
ADD_SUBDIRECTORY(jpeg)

INSTALL(FILES CMakeLists.txt DESTINATION src/third_party)

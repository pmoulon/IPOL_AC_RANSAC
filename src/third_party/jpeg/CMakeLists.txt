SET(JPEG_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR} PARENT_SCOPE)
SET(JPEG_LIBRARIES jpeg PARENT_SCOPE)

# Put here the source file name for the correct system-dependent memory
# manager file.  Usually jmemnobs.c, but you may want
# to use jmemansi.c or jmemname.c if you have limited swap space.
SET(SYSDEPMEM jmemnobs.c)

# library object files common to compression and decompression
SET(COMSRC jaricom.c jcomapi.c jutils.c jerror.c jmemmgr.c ${SYSDEPMEM})
# compression library source files
SET(CLIBSRC jcapimin.c jcapistd.c jcarith.c jctrans.c jcparam.c
            jdatadst.c jcinit.c jcmaster.c jcmarker.c jcmainct.c jcprepct.c
            jccoefct.c jccolor.c jcsample.c jchuff.c jcdctmgr.c jfdctfst.c
            jfdctflt.c jfdctint.c)
# decompression library source files
SET(DLIBSRC jdapimin.c jdapistd.c jdarith.c jdtrans.c jdatasrc.c
            jdmaster.c jdinput.c jdmarker.c jdhuff.c jdmainct.c
            jdcoefct.c jdpostct.c jddctmgr.c jidctfst.c jidctflt.c
            jidctint.c jdsample.c jdcolor.c jquant1.c jquant2.c jdmerge.c)
# These objectfiles are included in libjpeg.a
SET(LIBSRC ${CLIBSRC} ${DLIBSRC} ${COMSRC})

ADD_LIBRARY(jpeg ${LIBSRC})

IF(MSVC)
    ADD_DEFINITIONS(-D_CRT_SECURE_NO_WARNINGS -D_SCL_SECURE_NO_WARNINGS)
ENDIF(MSVC)

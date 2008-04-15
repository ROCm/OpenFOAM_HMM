TEMPLATE = lib
windows-vc.net:TEMPLATE = vclib
windows-vc.net:config += release
TARGET = adf
PATHTOSRC = ../

include($$PATHTOSRC/config/ccm.pro)
include(version.pro)

HEADERS += \
	ADF.h \
	ADF_internals.h 

SOURCES += \
  ADF_fortran_2_c.c \
	ADF_interface.c \
	ADF_internals.c

windows-x86:dll:QMAKE_LFLAGS_RELEASE -= $$QMAKE_LFLAGS_RELEASE
windows-vc.net:DEFINES += WINDOWS_X86

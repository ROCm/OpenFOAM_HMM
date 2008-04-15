TEMPLATE = lib
windows-vc.net:TEMPLATE = vclib
windows-vc.net:config += release
TARGET = cgns
PATHTOSRC = ../

include($$PATHTOSRC/config/ccm.pro)
include(version.pro)

HEADERS += \
	cgns_header.h \
	cgnslib.h \
	cgnslib_f.h \
	fortran_macros.h

SOURCES += \
	adf_cond.c \
	adf_ftoc.c \
	cg_ftoc.c \
	cgns_error.c \
	cgns_internals.c \
	cgnslib.c

LIBS += -L$$CCMLIB -ladf

# Some system need to be told what fortran mangles its names.
aix64_5.1-pwr4:DEFINES += LOWERCASE
aix_5.1-pwr4:DEFINES += LOWERCASE

windows-x86:dll:QMAKE_LFLAGS_RELEASE -= $$QMAKE_LFLAGS_RELEASE
windows-vc.net:DEFINES += WINDOWS_X86

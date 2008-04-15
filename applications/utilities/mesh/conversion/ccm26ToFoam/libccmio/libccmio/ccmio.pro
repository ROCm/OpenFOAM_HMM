TEMPLATE = lib
windows-vc.net:TEMPLATE = vclib
TARGET   = ccmio
PATHTOSRC = ../

include($$PATHTOSRC/config/ccm.pro)
include(version.pro)

HEADERS	+= \
  ccmio.h \
  ccmiocore.h \
  ccmiotypes.h \
  ccmioprivate.h \
  ccmioutility.h \
  vector.h \
  ccmioversion.h

SOURCES += \
  ccmio.c \
  ccmiocore.c \
  ccmioprivate.c \
  ccmioutility.c \
  vector.c \
  ccmioversion.c

aix64_5.1-pwr4:dll:LIBS += -L../lib/$$MK_BUILDNAME/$$RELEASEMODE -ladf
aix_5.1-pwr4:dll:LIBS += -L../lib/$$MK_BUILDNAME/$$RELEASEMODE -ladf
windows-x86-gcc:dll:LIBS += -L../lib/$$MK_BUILDNAME/$$RELEASEMODE -ladf
windows-x86:dll:LIBS += ../lib/$$MK_BUILDNAME/$$RELEASEMODE/adf501000.dll
windows-x86:dll:QMAKE_LFLAGS_RELEASE -= $$QMAKE_LFLAGS_RELEASE
windows-vc.net:INCLUDEPATH += ..\libadf

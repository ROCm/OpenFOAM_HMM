all: \
  ../libccmio/Makefile \
  ../libccmio/ccmio.c \
  ../libccmio/ccmio.h \
  ../libccmio/ccmio.pro \
  ../libccmio/ccmiocore.c \
  ../libccmio/ccmiocore.h \
  ../libccmio/ccmioprivate.c \
  ../libccmio/ccmioprivate.h \
  ../libccmio/ccmioread.c \
  ../libccmio/ccmiotypes.h \
  ../libccmio/ccmioutility.c \
  ../libccmio/ccmioutility.h \
  ../libccmio/ccmioverifymesh.c \
  ../libccmio/ccmioverifymesh.h \
  ../libccmio/ccmioversion.c \
  ../libccmio/ccmioversion.h \
  ../libccmio/vector.c \
  ../libccmio/vector.h \
  compatibility.dox \
  examples.dox \
  hierarchy.dox \
  mainpage.dox \
  errors.dox \
  groups.dox \
  intropage.dox
	doxygen doxygen.conf
	cp hierarchy.png html
	touch all

copy clean qmake test: force

force:

#
# Automatic setting of emacs local variables.
# Local Variables:
# tab-width: 2
# End:

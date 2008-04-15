all copy clean qmake test:
	if [ -d libadf ];    then cd libadf;    ${MAKE} -f Makefile.qmake $@;  fi
	if [ -d libccmio ];  then cd libccmio;  ${MAKE} -f Makefile.qmake $@;  fi
	if [ -d libcgns ];   then cd libcgns;   ${MAKE} -f Makefile.qmake $@;  fi

install:
	@if [ -z "$(MK_BUILDNAME)" ]; then \
		echo "MK_BUILDNAME must be set to the computer type you are installing."; \
	else \
		set -x; \
		echo "MK_BUILDNAME=$${MK_BUILDNAME}"; \
		./Makefile.install.sc $${MK_BUILDNAME}; \
	fi

docs:
	if [ -d docs ];      then cd docs;      ${MAKE} all; fi

# Internal target suitable for batch building
batch:	all
	if [ -d libadf ];    then cd libadf;    RELEASE=1 ${MAKE} -f Makefile.qmake all;  fi
	if [ -d libccmio ];  then cd libccmio;  RELEASE=1 ${MAKE} -f Makefile.qmake all;  fi
	if [ -d libcgns ];   then cd libcgns;   RELEASE=1 ${MAKE} -f Makefile.qmake all;  fi
	if [ -d libadf ];    then cd libadf;    RELEASE=1 STATIC=1 ${MAKE} -f Makefile.qmake all;  fi
	if [ -d libccmio ];  then cd libccmio;  RELEASE=1 STATIC=1 ${MAKE} -f Makefile.qmake all;  fi
	if [ -d libcgns ];    then cd libcgns;   RELEASE=1 STATIC=1 ${MAKE} -f Makefile.qmake all;  fi
	rm -f *.a
	rm -f *.so
	-# Copy all the .a and .so from each directory in lib.  This should
	-# only be one directory, but cp lib/*/release-shared/lib*.so doesn't
	-# work on AIX and HPUX.  Also copy the static libraries first, because
	-# they are always named .a, but the shared objects may fail on AIX
	libdirs=`ls lib`; for dir in $$libdirs;  do cp lib/$${dir}/release-static/lib*.a .; cp lib/$${dir}/release-shared/lib*.s* .; done

dist:	docs
	./makedist

#
# Automatic setting of emacs local variables.
# Local Variables:
# tab-width: 2
# End:

include($${PATHTOSRC}version.pro)
include($${PATHTOSRC}config/config.pro)

aix_4.3-com:CONFIG += thread
aix64_4.3-pwr3:CONFIG += thread

# Put the current build name in CONFIG so I can branch off it.
CONFIG += $$MK_BUILDNAME

contains(RELEASEMODE, release) {
    CONFIG -= debug
    CONFIG += release
    message("Building in release mode")
}
contains(RELEASEMODE, debug) {
    CONFIG -= release
    CONFIG += debug
    message("Building in Debug mode")
}

CONFIG -= qt
equals( LIBMODE, dynamic ) {
  CONFIG -= staticlib
  CONFIG += dll
  windows-x86-gcc:QMAKE_EXTENSION_SHLIB = dll
  message( "Building dynamic libraries" )
} else {
  CONFIG -= dll
  CONFIG += staticlib
  message( "Building static libraries" )
}

dll:RELEASEMODE = $$RELEASEMODE-shared
else:RELEASEMODE = $$RELEASEMODE-static

#equals( TESTMODE, test ) {
#    DEFINES += TESTING 
#    RELEASEMODE = $${RELEASEMODE}-testing
#    message( "Building in test mode" )
#}


# This define is needed to compile cgns correctly on HP.  WRO 2005-Feb-28
hpux_11.00-pa8000: DEFINES += LOWERCASE
hpux64_11.00-pa8000: DEFINES += LOWERCASE
hpux64_11.22-itanium2: DEFINES += LOWERCASE


# Is profiling requested?
!isEmpty(PROFILEMODE) {
    RELEASEMODE = $$RELEASEMODE$$PROFILEMODE
    QMAKE_CFLAGS   += $$PROFILEMODE
    QMAKE_CXXFLAGS += $$PROFILEMODE
    QMAKE_LFLAGS   += $$PROFILEMODE
    message( "Building in profile mode: $$PROFILEMODE" )
}

# Where is qmake?
QMAKE_QMAKE = $${PATHTOSRC}config/$$MK_BUILDNAME/qmake

########################## Compiler Options #######################
# Allow user to override the default compiler options with
# environmental variables.
!isEmpty( CC ) {
  QMAKE_CC		= $$CC
}
!isEmpty( CXX ) {
  equals( QMAKE_LINK, $$QMAKE_CXX ) {
    QMAKE_LINK = $$CXX
  }
  equals( QMAKE_LINK_SHLIB, $$QMAKE_CXX ) {
    QMAKE_LINK_SHLIB = $$CXX
  }
  QMAKE_CXX		= $$CXX
}

##################################################################

CCMBASE = $$PATHTOSRC
CCMINC = $$CCMBASE/include
CCMLIB = $$CCMBASE/lib/$$MK_BUILDNAME/$$RELEASEMODE
CCMBIN = $$CCMBASE/bin/$$MK_BUILDNAME/$$RELEASEMODE
CCMOBJ = $$CCMBASE/obj/$$MK_BUILDNAME/$$RELEASEMODE

##################################################################

INCLUDEPATH += \
           $$PATHTOSRC \
           $$CCMINC

#DEPENDPATH += $$CCMINC

DESTDIR = $$CCMLIB
contains(TEMPLATE,app) {
  DESTDIR = $$CCMBIN
}

OBJECTS_DIR = $$CCMOBJ


#headers.path = $$CCMINC
#headers.files = *.h
#INSTALLS += headers

#message ( "CONFIG=$$CONFIG" )
message ( "MK_BUILDNAME=$$MK_BUILDNAME; export MK_BUILDNAME;" )
message ( "RELEASEMODE=$$RELEASEMODE; export RELEASEMODE;" )

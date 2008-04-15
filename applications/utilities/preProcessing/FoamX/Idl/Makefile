#------------------------------------------------------------------------------
# Makefile for FoamX Client-Side Idl Java Code.
#------------------------------------------------------------------------------

JAVA_IDL = idlj
JAVA_IDLFLAGS = -td ../Java -fclient
FOAMXSERVER_JAVA_DIR = ../Java/FoamXServer
FOAMXSERVER_JAVA_FILE = $(FOAMXSERVER_JAVA_DIR)/FoamXError.java

MICO_IDL    = $(MICO_ARCH_PATH)/bin/idl
MICO_INCDIR = $(MICO_ARCH_PATH)/include
MICO_IDLFLAGS = --poa --c++-suffix=C --hh-suffix=H --any --typecode
FOAMXLIB_DIR = ../C++/FoamXServer
FOAMXLIB_H_FILE = $(FOAMXLIB_DIR)/FoamX.H

IDLFILE  = FoamX.idl

all: $(FOAMXSERVER_JAVA_FILE) $(FOAMXLIB_H_FILE)

$(FOAMXSERVER_JAVA_FILE): $(IDLFILE)
	$(JAVA_IDL) $(JAVA_IDLFLAGS) $(IDLFILE)

$(FOAMXLIB_H_FILE): $(IDLFILE)
	$(MICO_IDL) $(MICO_IDLFLAGS) $(IDLFILE)
	mv FoamX.[HC] $(FOAMXLIB_DIR)

clean:
	rm -rf $(FOAMXSERVER_JAVA_DIR)
	rm -rf $(FOAMXLIB_DIR)/FoamX.[HC]

#------------------------------------------------------------------------------

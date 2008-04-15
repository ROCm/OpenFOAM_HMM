include ../../config/Makeflags.sys

OBJ_DIR        = $(BUILD_HOME)/libccmio/libccmio
TARGET_SO_NAME = ccmio
SOURCES        = $(wildcard *.cpp)

INCLUDE_PATH = \
 -I$(STAR_HOME)/base/src \
 -I$(STAR_HOME)/libccmio

# Because .cpp files #include .c files, need to
# allow dependency cache to deal with removed .c files
%.c:
	@echo '# Ignoring missing $@'

include $(STAR_HOME)/config/Makefile.shared_object
include $(STAR_HOME)/config/Makefile.objects

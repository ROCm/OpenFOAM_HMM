include ../../config/Makeflags.sys

OBJ_DIR        = $(BUILD_HOME)/libccmio
TARGET_SO_NAME = cgns
SOURCES        = $(wildcard *.c)
SOURCES        := $(filter-out au1_ftoc.c, $(SOURCES))

INCLUDE_PATH := \
 -I$(STAR_HOME)/base/src\
 -I$(STAR_HOME)/libccmio

include $(STAR_HOME)/config/Makefile.c.shared_object
include $(STAR_HOME)/config/Makefile.c.objects

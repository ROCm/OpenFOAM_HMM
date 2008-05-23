include ../../../Makefile.in

INCLUDES = $(INCDIR) -I./
DEBUGFLAGS =

ifeq ($(ddmalloc),yes)
  INCLUDES := $(INCLUDES) -I$(DMALLOCDIR)/include
  DEBUGFLAGS := $(DEBUGFLAGS) -DDMALLOC -DDEBUG
  OPTFLAGS = -g
endif

ifeq ($(dmalloc),yes)
  INCLUDES := $(INCLUDES) -I$(DMALLOCDIR)/include
  DEBUGFLAGS := $(DEBUGFLAGS) -DDMALLOC
endif

ifeq ($(debug),yes)
  DEBUGFLAGS := $(DEBUGFLAGS) -DDEBUG
  OPTFLAGS = -g
endif


CFLAGS = $(COPTIONS) $(OPTFLAGS) $(DEBUGFLAGS) $(INCLUDES)


OBJS = fused.o comm.o util.o memory.o debug.o setup.o diffuse.o grsetup.o timer.o \
       coarsen.o edge_refine.o drivers.o remap.o iidxsort.o ikeysort.o var.o

.c.o:
	$(PARCC) $(CFLAGS) -c $*.c


../libIMparmetis.a: $(OBJS)
	$(AR) $@ $(OBJS)
	$(RANLIB) $@

clean:
	rm -f *.o

realclean:
	rm -f *.o ; rm -f ../libIMparmetis.a

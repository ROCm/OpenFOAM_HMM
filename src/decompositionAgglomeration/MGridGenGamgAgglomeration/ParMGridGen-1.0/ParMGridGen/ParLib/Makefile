include ../../Makefile.in

DEBUGFLAGS =

ifeq ($(ddmalloc),yes)
  INCDIR := $(INCDIR) -I$(DMALLOCDIR)/include
  DEBUGFLAGS := $(DEBUGFLAGS) -DDMALLOC -DDEBUG
  OPTFLAGS = -g
endif

ifeq ($(dmalloc),yes)
  INCDIR := $(INCDIR)  -I$(DMALLOCDIR)/include
  DEBUGFLAGS := $(DEBUGFLAGS) -DDMALLOC  
endif

ifeq ($(debug),yes)
  DEBUGFLAGS := $(DEBUGFLAGS) -DDEBUG  
  OPTFLAGS = -g
endif

INCLUDES = -I./ -I$(IMLIBDIR) $(INCDIR)
CFLAGS = $(COPTIONS) $(OPTFLAGS) $(DEBUGFLAGS) $(INCLUDES)

OBJS = comm.o util.o memory.o debug.o setup.o parmgridgen.o grsetup.o move.o \
       ikeysort.o
IPARMETISOBJS = fused.o comm.o util.o memory.o debug.o setup.o diffuse.o grsetup.o timer.o \
        coarsen.o edge_refine.o drivers.o remap.o iidxsort.o ikeysort.o var.o
IMPARMETISOBJS = $(addprefix ../IMParMetis-2.0/ParMETISLib/, $(IPARMETISOBJS))


.c.o:
	$(PARCC) $(CFLAGS) -c $*.c


../../libparmgrid.a: $(OBJS) $(IMPARMETISOBJS)
	$(AR) $@ $(OBJS) $(IMPARMETISOBJS)
	$(RANLIB) $@

clean:
	rm -f *.o 

realclean:
	rm -f *.o ; rm -f ../../libparmgrid.a

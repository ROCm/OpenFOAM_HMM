include ../../Makefile.in

DEBUGFLAGS = 

ifeq ($(ddmalloc),yes)
  INCDIR := $(INCDIR) -I$(DMALLOCDIR)/include
  DEBUGFLAGS := $(DEBUGFLAGS) -DDMALLOC -DDEBUG
  OPTFLAGS = -g
endif

ifeq ($(dmalloc),yes)
  INCDIR := $(INCDIR) -I$(DMALLOCDIR)/include
  DEBUGFLAGS := $(DEBUGFLAGS) -DDMALLOC  
endif


ifeq ($(debug),yes)
  DEBUGFLAGS := $(DEBUGFLAGS) -DDEBUG  
  OPTFLAGS = -g
endif

INCLUDES = -I./ -I$(IMLIBDIR) $(INCDIR)
CFLAGS = $(COPTIONS) $(OPTFLAGS) $(DEBUGFLAGS) $(INCLUDES)

OBJS = coarsen.o kwayfm.o match.o merge.o mgridgen.o refine.o setup.o 
IOBJS = blas.o dfkeysort.o dkeysort.o file.o ifkeysort.o ifloatsort.o \
       iintsort.o ikeysort.o memory.o  sort.o util.o
IMOBJS = $(addprefix ../IMlib/, $(IOBJS)) 


.c.o:
	$(CC) $(CFLAGS) -c $*.c

../../libmgrid.a: $(OBJS) $(IMOBJS)
	$(AR) $@  $(OBJS) $(IMOBJS)
	$(RANLIB) $@

clean:
	rm -f *.o

realclean:
	rm -f *.o ; rm -f ../../libmgrid.a

include ../../Makefile.in

#--------------------------------------------------------------------------
DMALLOCINC = 
DMALLOCFLAGS =
DEBUGFLAGS =

ifeq ($(ddmalloc),yes)
  DMALLOCINC = -I$(DMALLOCDIR)/include
  DMALLOCFLAGS = -DDMALLOC
  DEBUGFLAGS = -DDEBUG
  OPTFLAGS = -g
endif

ifeq ($(dmalloc),yes)
  DMALLOCINC = -I$(DMALLOCDIR)/include
  DMALLOCFLAGS = -DDMALLOC
  OPTFLAGS = -g
endif

ifeq ($(debug),yes)
  DEBUGFLAGS = -DDEBUG
  OPTFLAGS = -g
endif

INCLUDES = -I./ $(INCDIR) $(DMALLOCINC)
CFLAGS = $(COPTIONS) $(OPTFLAGS) $(DMALLOCFLAGS) $(DEBUGFLAGS) $(INCLUDES)
#--------------------------------------------------------------------------

OBJS = blas.o dfkeysort.o dkeysort.o file.o ifkeysort.o ifloatsort.o \
       iintsort.o ikeysort.o memory.o  sort.o util.o

.c.o:
	$(CC) $(CFLAGS) -c $*.c

./libIMlib.a: $(OBJS)
	$(AR) $@ $(OBJS)
	$(RANLIB) $@

clean:
	rm -f *.o

realclean:
	rm -f *.o ; rm -f ./libIMlib.a

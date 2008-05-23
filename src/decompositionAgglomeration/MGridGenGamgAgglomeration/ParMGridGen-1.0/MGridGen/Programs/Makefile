include ../../Makefile.in

MGRIDLIB  = ../../libmgrid.a

ifeq ($(ddmalloc),yes)
  DEBUGFLAGS := $(DEBUGFLAGS) -DDMALLOC -DDEBUG
  INCDIR := $(INCDIR) -I$(DMALLOCDIR)/include
  LIBDIR := $(LIBDIR) -L$(DMALLOCDIR)/lib
  LIBS := $(LIBS) -ldmalloc
  OPTFLAGS = -g
endif

ifeq ($(dmalloc),yes)
  DEBUGFLAGS := $(DEBUGFLAGS) -DDMALLOC  
  INCDIR := $(INCDIR) -I$(DMALLOCDIR)/include
  LIBDIR := $(LIBDIR) -L$(DMALLOCDIR)/lib
  LIBS := $(LIBS) -ldmalloc
endif

ifeq ($(debug),yes)
  DEBUGFLAGS := $(DEBUGFLAGS) -DDEBUG  
  OPTFLAGS = -g
endif

INCLUDES = $(INCDIR) -I$(IMLIBDIR) -I../Lib
CFLAGS = $(COPTIONS) $(OPTFLAGS) $(DEBUGFLAGS) $(INCLUDES)
LD = $(CC) $(LDOPTIONS) $(LIBDIR)


MGRIDOBJS = mgridgen.o io.o
TRANSOBJS = transform.o io.o
EXTROBJS  = extract.o io.o


.c.o:
	$(CC) $(CFLAGS) -c $*.c

default: $(BINDIR)/mgridgen #$(BINDIR)/transform $(BINDIR)/extract

$(BINDIR)/mgridgen: $(MGRIDOBJS) $(MGRIDLIB) 
	$(LD)  -o $@ $(MGRIDOBJS) $(LIBS)
	chmod 744 $@

## $(BINDIR)/transform : $(TRANSOBJS) $(MGRIDLIB)
## 	$(LD) -o $@ $(TRANSOBJS) $(LIBS)

## $(BINDIR)/extract : $(EXTROBJS) $(MGRIDLIB)
##	$(LD) -o $@ $(EXTROBJS) $(LIBS)


clean:
	rm -f *.o

realclean:
	rm -f *.o ;\
        rm -f $(BINDIR)/mgridgen

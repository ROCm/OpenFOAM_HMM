include ../../Makefile.in

PARMGRIDLIB = ../../libparmgrid.a 

ifeq ($(ddmalloc),yes)
  DEBUGFLAGS := $(DEBUGFLAGS) -DDMALLOC -DDEBUG
  INCDIR := $(INCDIR)  -I$(DMALLOCDIR)/include
  LIBDIR := $(LIBDIR) -L$(DMALLOCDIR)/lib
  PARLIBS := $(PARLIBS) -ldmalloc
  OPTFLAGS = -g
endif

ifeq ($(dmalloc),yes)
  DEBUGFLAGS := $(DEBUGFLAGS) -DDMALLOC
  INCDIR := $(INCDIR)  -I$(DMALLOCDIR)/include 
  LIBDIR := $(LIBDIR) -L$(DMALLOCDIR)/lib
  PARLIBS := $(PARLIBS) -ldmalloc
endif

ifeq ($(debug),yes)
  DEBUGFLAGS := $(DEBUGFLAGS) -DDEBUG
  OPTFLAGS = -g
endif

INCLUDES = -I../ParLib -I$(IMLIBDIR) $(INCDIR)
CFLAGS = $(COPTIONS) $(OPTFLAGS) $(DEBUGFLAGS) $(INCLUDES)

OBJS = main.o parmgridgen.o io.o 

         
.c.o:
	$(PARCC) $(CFLAGS) -c $*.c


default: $(BINDIR)/parmgridgen

$(BINDIR)/parmgridgen: $(OBJS) $(PARMGRIDLIB)
	$(PARLD)  -o $@ $(OBJS) $(LIBDIR) $(PARLIBS)
	chmod 744 $@

clean:
	rm -f *.o 

realclean:
	rm -f *.o ;\
	rm -f $(BINDIR)/parmgridgen

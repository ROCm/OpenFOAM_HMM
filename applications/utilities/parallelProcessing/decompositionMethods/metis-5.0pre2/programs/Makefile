include ../Makefile.in

CURBUILDDIR = $(PRGBUILDDIR)


METISSRC       = metis.c io.c smbfactor.c 
PMETISSRC      = pmetis.c io.c cmdline_pmetis.c 
KMETISSRC      = kmetis.c io.c 
OEMETISSRC     = oemetis.c io.c smbfactor.c 
ONMETISSRC     = onmetis.c io.c smbfactor.c 
MESH2DUALSRC   = mesh2dual.c io.c 
MESH2NODALSRC  = mesh2nodal.c io.c 
PARTDMESHSRC   = partdmesh.c io.c 
PARTNMESHSRC   = partnmesh.c io.c 
GRAPHCHKSRC    = graphchk.c io.c 
KFMETISSRC     = kfmetis.c io.c cmdline_kfmetis.c smbfactor.c
CMETISSRC      = cmetis.c io.c cmdline_cmetis.c 

METISOBJS      = $(patsubst %.c, $(CURBUILDDIR)/%$(OBJEXT), $(METISSRC))
PMETISOBJS     = $(patsubst %.c, $(CURBUILDDIR)/%$(OBJEXT), $(PMETISSRC))
KMETISOBJS     = $(patsubst %.c, $(CURBUILDDIR)/%$(OBJEXT), $(KMETISSRC))
OEMETISOBJS    = $(patsubst %.c, $(CURBUILDDIR)/%$(OBJEXT), $(OEMETISSRC))
ONMETISOBJS    = $(patsubst %.c, $(CURBUILDDIR)/%$(OBJEXT), $(ONMETISSRC))
MESH2DUALOBJS  = $(patsubst %.c, $(CURBUILDDIR)/%$(OBJEXT), $(MESH2DUALSRC))
MESH2NODALOBJS = $(patsubst %.c, $(CURBUILDDIR)/%$(OBJEXT), $(MESH2NODALSRC))
PARTDMESHOBJS  = $(patsubst %.c, $(CURBUILDDIR)/%$(OBJEXT), $(PARTDMESHSRC))
PARTNMESHOBJS  = $(patsubst %.c, $(CURBUILDDIR)/%$(OBJEXT), $(PARTNMESHSRC))
GRAPHCHKOBJS   = $(patsubst %.c, $(CURBUILDDIR)/%$(OBJEXT), $(GRAPHCHKSRC))
KFMETISOBJS    = $(patsubst %.c, $(CURBUILDDIR)/%$(OBJEXT), $(KFMETISSRC))
CMETISOBJS     = $(patsubst %.c, $(CURBUILDDIR)/%$(OBJEXT), $(CMETISSRC))

ALLOBJS  = $(patsubst %.c, $(CURBUILDDIR)/%$(OBJEXT), $(wildcard *.c))
HEADERS  = $(wildcard *.h) $(TOPDIR)/include/metis.h $(wildcard $(GKLIBINCDIR)/*.h)

LIBRARIES = $(BUILDDIR)/libmetis$(LIBEXT)


TARGETS = $(BUILDDIR)/metis$(VERNUM)$(EXEEXT) \
          $(BUILDDIR)/pmetis$(VERNUM)$(EXEEXT) \
          $(BUILDDIR)/kmetis$(VERNUM)$(EXEEXT) \
          $(BUILDDIR)/oemetis$(VERNUM)$(EXEEXT) \
          $(BUILDDIR)/onmetis$(VERNUM)$(EXEEXT) \
          $(BUILDDIR)/mesh2dual$(VERNUM)$(EXEEXT) \
          $(BUILDDIR)/mesh2nodal$(VERNUM)$(EXEEXT) \
          $(BUILDDIR)/partdmesh$(VERNUM)$(EXEEXT) \
          $(BUILDDIR)/partnmesh$(VERNUM)$(EXEEXT) \
          $(BUILDDIR)/graphchk$(VERNUM)$(EXEEXT) \
          $(BUILDDIR)/kfmetis$(VERNUM)$(EXEEXT) \
          $(BUILDDIR)/cmetis$(VERNUM)$(EXEEXT) 


default: $(TARGETS)


$(BUILDDIR)/metis$(VERNUM)$(EXEEXT): $(METISOBJS) $(LIBRARIES)
	$(LD) $(LDOPTIONS) $(EXEOUTPUTFILE) $(METISOBJS) $(LIBSDIR) $(LIBS)
	chmod 744 $@
	@if [ "$(BINDIR)" ]; then cp $@* $(BINDIR); fi

$(BUILDDIR)/pmetis$(VERNUM)$(EXEEXT): $(PMETISOBJS) $(LIBRARIES)
	$(LD) $(LDOPTIONS) $(EXEOUTPUTFILE) $(PMETISOBJS) $(LIBSDIR) $(LIBS)
	chmod 744 $@
	@if [ "$(BINDIR)" ]; then cp $@* $(BINDIR); fi

$(BUILDDIR)/kmetis$(VERNUM)$(EXEEXT): $(KMETISOBJS) $(LIBRARIES)
	$(LD) $(LDOPTIONS) $(EXEOUTPUTFILE) $(KMETISOBJS) $(LIBSDIR) $(LIBS)
	chmod 744 $@
	@if [ "$(BINDIR)" ]; then cp $@* $(BINDIR); fi

$(BUILDDIR)/oemetis$(VERNUM)$(EXEEXT): $(OEMETISOBJS) $(LIBRARIES)
	$(LD) $(LDOPTIONS) $(EXEOUTPUTFILE) $(OEMETISOBJS) $(LIBSDIR) $(LIBS)
	chmod 744 $@
	@if [ "$(BINDIR)" ]; then cp $@* $(BINDIR); fi

$(BUILDDIR)/onmetis$(VERNUM)$(EXEEXT): $(ONMETISOBJS) $(LIBRARIES)
	$(LD) $(LDOPTIONS) $(EXEOUTPUTFILE) $(ONMETISOBJS) $(LIBSDIR) $(LIBS)
	chmod 744 $@
	@if [ "$(BINDIR)" ]; then cp $@* $(BINDIR); fi

$(BUILDDIR)/mesh2dual$(VERNUM)$(EXEEXT): $(MESH2DUALOBJS) $(LIBRARIES)
	$(LD) $(LDOPTIONS) $(EXEOUTPUTFILE) $(MESH2DUALOBJS) $(LIBSDIR) $(LIBS)
	chmod 744 $@
	@if [ "$(BINDIR)" ]; then cp $@* $(BINDIR); fi

$(BUILDDIR)/mesh2nodal$(VERNUM)$(EXEEXT): $(MESH2NODALOBJS) $(LIBRARIES)
	$(LD) $(LDOPTIONS) $(EXEOUTPUTFILE) $(MESH2NODALOBJS) $(LIBSDIR) $(LIBS)
	chmod 744 $@
	@if [ "$(BINDIR)" ]; then cp $@* $(BINDIR); fi

$(BUILDDIR)/partdmesh$(VERNUM)$(EXEEXT): $(PARTDMESHOBJS) $(LIBRARIES)
	$(LD) $(LDOPTIONS) $(EXEOUTPUTFILE) $(PARTDMESHOBJS) $(LIBSDIR) $(LIBS)
	chmod 744 $@
	@if [ "$(BINDIR)" ]; then cp $@* $(BINDIR); fi

$(BUILDDIR)/partnmesh$(VERNUM)$(EXEEXT): $(PARTNMESHOBJS) $(LIBRARIES)
	$(LD) $(LDOPTIONS) $(EXEOUTPUTFILE) $(PARTNMESHOBJS) $(LIBSDIR) $(LIBS)
	chmod 744 $@
	@if [ "$(BINDIR)" ]; then cp $@* $(BINDIR); fi

$(BUILDDIR)/graphchk$(VERNUM)$(EXEEXT): $(GRAPHCHKOBJS) $(LIBRARIES)
	$(LD) $(LDOPTIONS) $(EXEOUTPUTFILE) $(GRAPHCHKOBJS) $(LIBSDIR) $(LIBS)
	chmod 744 $@
	@if [ "$(BINDIR)" ]; then cp $@* $(BINDIR); fi

$(BUILDDIR)/kfmetis$(VERNUM)$(EXEEXT): $(KFMETISOBJS) $(LIBRARIES)
	$(LD) $(LDOPTIONS) $(EXEOUTPUTFILE) $(KFMETISOBJS) $(LIBSDIR) $(LIBS)
	chmod 744 $@
	@if [ "$(BINDIR)" ]; then cp $@* $(BINDIR); fi

$(BUILDDIR)/cmetis$(VERNUM)$(EXEEXT): $(CMETISOBJS) $(LIBRARIES)
	$(LD) $(LDOPTIONS) $(EXEOUTPUTFILE) $(CMETISOBJS) $(LIBSDIR) $(LIBS)
	chmod 744 $@
	@if [ "$(BINDIR)" ]; then cp $@* $(BINDIR); fi



$(BUILDDIR)/cepic$(VERNUM)$(EXEEXT): $(CEPICOBJS) $(LIBRARIES)
	$(LD)  -o $@ $(CEPICOBJ) $(LIBS)
	chmod 744 $@
	@if [ "$(BINDIR)" ]; then cp $@* $(BINDIR); fi

$(BUILDDIR)/cepic-rcb$(VERNUM)$(EXEEXT): $(CEPICRCBOBJS) $(LIBRARIES)
	$(LD)  -o $@ $(CEPICRCBOBJ) $(LIBS)
	chmod 744 $@
	@if [ "$(BINDIR)" ]; then cp $@* $(BINDIR); fi

$(BUILDDIR)/mmetis$(VERNUM)$(EXEEXT): $(MMETISOBJS) $(LIBRARIES)
	$(LD)  -o $@ $(MMETISOBJS) $(LIBS)
	chmod 744 $@
	@if [ "$(BINDIR)" ]; then cp $@* $(BINDIR); fi



clean:
	rm -f $(ALLOBJS) 

realclean:
	rm -f $(ALLOBJS) $(TARGETS)



$(ALLOBJS) : $(HEADERS) ../Makefile.in Makefile

$(CURBUILDDIR)/%$(OBJEXT) : %.c
	$(CC) $(CFLAGS) $(SOURCEFILE) $(OUTPUTFILE)



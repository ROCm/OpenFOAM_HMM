
default:
	more Makefile.help

all:
	$(MAKE) -C GKlib/trunk
	$(MAKE) -C libmetis
	$(MAKE) -C programs
	$(MAKE) -C test

clean:
	$(MAKE) -C GKlib/trunk clean
	$(MAKE) -C libmetis clean
	$(MAKE) -C programs clean
	$(MAKE) -C test clean

realclean:
	$(MAKE) -C GKlib/trunk realclean
	$(MAKE) -C libmetis realclean
	$(MAKE) -C programs realclean
	$(MAKE) -C test realclean

dist:
	$(MAKE) -C dist dist

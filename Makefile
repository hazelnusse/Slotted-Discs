package = slotted-discs
version = 0.1
tarname = $(package)
distdir = $(tarname)-$(version)

all clean render simulate:
	cd src && $(MAKE) $@

dist: $(distdir).tar.bz2

$(distdir).tar.bz2: $(distdir)
	tar cjf $@ $(distdir)
	rm -rf $(distdir)

$(distdir): FORCE
	mkdir -p $(distdir)/src
	cp Makefile $(distdir)
	cp src/Makefile $(distdir)/src
	cp src/*.cpp $(distdir)/src
	cp src/*.c $(distdir)/src
	cp src/*.h $(distdir)/src

FORCE:
	rm -rf $(distdir).tar.bz2
	rm -rf $(distdir)

distclean:
	rm -rf $(distdir).tar.bz2

.PHONY : all clean dist distclean FORCE

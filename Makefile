#
# Makefile.local must define the following variables
#  BINDIR
#  LIBDIR
#  INCDIR
#  PYDIR
#  CPP
#  MPICPP
#
# See Makefile.local.example for an example (in fact, you may just
# be able to do 'cp Makefile.local.example Makefile.local)
#
include Makefile.local

EXE_INSTALL=frb-compare
EXE_NOINSTALL=test-rng write-pulse
SCRIPT_INSTALL=frb-dump.py frb-compare-postprocess.py
PY_INSTALL=frb_olympics.py

all: libfrb_olympics.so frb_olympics_c.so $(EXE_INSTALL) $(EXE_NOINSTALL)

%.o: %.cpp frb_olympics.hpp
	$(CPP) -c -o $@ $<

# use suffix .mo for MPI object files
%.mo: %.cpp frb_olympics.hpp
	$(MPICPP) -c -o $@ $<

libfrb_olympics.so: frb_misc.o frb_pulse.o frb_rng.o frb_search_algorithm_base.o frb_search_params.o frb_downsample.o frb_rechunk.o frb_simple_direct.o frb_sloth.o frb_simple_tree.o frb_bonsai.o
	$(CPP) -o $@ -shared $^ -ljstree -lfftw3

frb_olympics_c.cpp: frb_olympics_c.pyx _frb_olympics_c.pxd frb_olympics.hpp
	cython --cplus $<

frb_olympics_c.so: frb_olympics_c.cpp libfrb_olympics.so
	$(CPP) -shared -o $@ $< -L. -lfrb_olympics

test-rng: test-rng.o libfrb_olympics.so
	$(CPP) -o $@ $^

write-pulse: write-pulse.o libfrb_olympics.so
	$(CPP) -o $@ $^

frb-compare: frb-compare.mo libfrb_olympics.so
	$(MPICPP) -o $@ $^

install: libfrb_olympics.so frb_olympics_c.so $(EXE_INSTALL)
	cp -f frb_olympics.hpp $(INCDIR)/frb_olympics.hpp
	cp -f libfrb_olympics.so $(LIBDIR)/libfrb_olympics.so
	cp -f frb_olympics_c.so $(PYDIR)/frb_olympics_c.so
	cp -f $(EXE_INSTALL) $(BINDIR)/
	cp -f $(SCRIPT_INSTALL) $(BINDIR)/
	cp -f $(PY_INSTALL) $(PYDIR)/

clean:
	rm -f *~ *.o *.mo *.so frb_olympics_c.cpp $(EXE_INSTALL) $(EXE_NOINSTALL)

uninstall:
	for f in $(EXE_INSTALL) $(SCRIPT_INSTALL); do rm -f $(BINDIR)/$$f; done
	for f in $(PY_INSTALL); do rm -f $(PYDIR)/$$f; done
	rm -f $(INCDIR)/frb_olympics.hpp $(LIBDIR)/libfrb_olympics.so $(PYDIR)/frb_olympics_c.so


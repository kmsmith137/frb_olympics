include Makefile.local

ifndef PYDIR
$(error Fatal: Makefile.local must define PYDIR variable)
endif

all:
	echo "make: nothing to compile here, the only targets are 'make install', 'make uninstall', and 'make clean'"

install: 
	cp -f frb_olympics.py $(PYDIR)/

clean:
	rm -f *~ *.pyc

uninstall:
	rm -f $(PYDIR)/frb_olympics.py

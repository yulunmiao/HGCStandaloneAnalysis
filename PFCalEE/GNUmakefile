# $Id: GNUmakefile,v 1.2 2000-10-19 12:22:10 stanaka Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := PFCalEE
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

CPPFLAGS += $(shell $(ROOTSYS)/bin/root-config --cflags) 
EXTRALIBS += $(shell $(ROOTSYS)/bin/root-config --glibs)

.PHONY: all
all: lib bin

dictionary:
	rootcint -v -f src/dict.cc -c include/HGCSSSimHit.hh include/LinkDef.h
	sed "s/include\/HGCSSSimHit/HGCSSSimHit/" src/dict.h > include/dict.h
	rm src/dict.h

include $(G4INSTALL)/config/binmake.gmk

visclean:
	rm -f g4*.prim g4*.eps g4*.wrl
	rm -f .DAWN_*


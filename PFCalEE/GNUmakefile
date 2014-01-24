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

CPPFLAGS += -Iuserlib/include/ $(shell $(ROOTSYS)/bin/root-config --cflags)
EXTRALIBS += $(shell $(ROOTSYS)/bin/root-config --glibs) -Luserlib/lib -lPFCalEEuserlib

SUBDIRS = userlib

.PHONY: $(SUBDIRS) all
all: $(SUBDIRS) lib bin

include $(G4INSTALL)/config/binmake.gmk

$(SUBDIRS):
	$(MAKE) -C $@

visclean:
	rm -f g4*.prim g4*.eps g4*.wrl
	rm -f .DAWN_*

userclean:
	$(MAKE) -C userlib/ clean

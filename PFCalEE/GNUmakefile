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

.PHONY: $(SUBDIRS) all
all: $(SUBDIRS) lib bin

include $(G4INSTALL)/config/binmake.gmk

INCFLAGS += -I$(HEPMC_DIR)/include/

ifdef G4LIB_USE_PYTHIA
      LDLIBS1 += -L$(HEPMC_DIR)/lib -lHepMC -lHepMCfio -L$(G4TMPDIR) -lPythia6 -lg2c
else
      LDLIBS1 += -L$(HEPMC_DIR)/lib -lHepMC -lHepMCfio #$(G4TMPDIR)/HEPEvtcom.o
endif


SUBDIRS = userlib

$(SUBDIRS):
	$(MAKE) -C $@

visclean:
	rm -f g4*.prim g4*.eps g4*.wrl
	rm -f .DAWN_*

userclean:
	$(MAKE) -C userlib/ clean

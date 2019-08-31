name := MAGNETOCOSMICS
G4TARGET := $(name)
G4EXLIB := true
TSYLIBDIR :=./fortran/lib
#USE_UNILIB :=1


VERSIONG4 :=$(notdir $(G4INSTALL))
VERSIONG4 :=$(subst geant4.,,$(VERSIONG4))
VERSIONG4 :=$(subst ., ,$(VERSIONG4))
SUBVERSIONG4 :=$(word 2,$(VERSIONG4))
VERSIONG4 :=$(word 1,$(VERSIONG4))
ifneq ($(findstring $(VERSIONG4), 1 2 3 4 5 6 ),)
  CPPFLAGS += -DBEFORE_V7
endif 

ifdef USE_OLD_GPS
CPPFLAGS += -DUSE_OLD_GPS
endif
ifdef USE_UNILIB
CPPFLAGS += -DUSE_UNILIB
endif

ifdef BUILD_STATIC
CPPFLAGS += -static
endif

ifdef USE_TSY04
CPPFLAGS += -DUSE_TSY04
endif
ifdef USE_PALEO
CPPFLAGS += -DUSE_PALEO
endif

CPPFLAGS += -Wno-deprecated -L$(TSYLIBDIR)
 
EXTRALIBS := -lFortranMagField  -lg2c
ifdef USE_UNILIB
CPPFLAGS += -I./unilib30/include  -L./unilib30/lib
EXTRALIBS := -lUnilib30 -lFortranMagField -lg2c 
endif
ifdef USE_PALEO
CPPFLAGS +=  -L./paleocode/lib
EXTRALIBS := -lFortranPaleoField -lFortranMagField  -lg2c
endif
 
  #fortran dir
FORTRAN_SUBDIRS :=  fortran/  

ifdef USE_PALEO
FORTRAN_SUBDIRS +=$(FORTRAN_SUBDIRS)paleocode/
endif

.PHONY: all
all: lib_fortran lib bin link bin_magcos 
clean_all: clean_fortran clean

lib_fortran:
	@for dir in $(FORTRAN_SUBDIRS); do (cd $$dir && $(MAKE)); done

clean_fortran:
	@for dir in $(FORTRAN_SUBDIRS); do (cd $$dir && $(MAKE) clean ); done



link: bin
ifdef BUILD_STATIC
include ./binmake_mod.gmk
else
include $(G4INSTALL)/config/binmake.gmk
endif


bin_magcos:
	@if [ ! -d bin ]; then mkdir bin; fi
ifdef BUILD_STATIC
	@cp $(G4BINDIR)/MAGNETOCOSMICS bin/STATIC_MAGNETOCOSMICS
else
	@if [ -f bin/MAGNETOCOSMICS ]; then rm bin/MAGNETOCOSMICS; fi
	@ln -s $(G4BINDIR)/MAGNETOCOSMICS bin/MAGNETOCOSMICS
endif

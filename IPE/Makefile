# GNU Makefile template for user ESMF application

################################################################################
################################################################################
## This Makefile must be able to find the "esmf.mk" Makefile fragment in the  ##
## 'include' line below. Following the ESMF User's Guide, a complete ESMF     ##
## installation should ensure that a single environment variable "ESMFMKFILE" ##
## is made available on the system. This variable should point to the         ##
## "esmf.mk" file.                                                            ##
##                                                                            ##
## This example Makefile uses the "ESMFMKFILE" environment variable.          ##
##                                                                            ##
## If you notice that this Makefile cannot find variable ESMFMKFILE then      ##
## please contact the person responsible for the ESMF installation on your    ##
## system.                                                                    ##
## As a work-around you can simply hardcode the path to "esmf.mk" in the      ##
## include line below. However, doing so will render this Makefile a lot less ##
## flexible and non-portable.                                                 ##
################################################################################

ifneq ($(origin ESMFMKFILE), environment)
$(error Environment variable ESMFMKFILE was not set.)
endif

include $(ESMFMKFILE)
IPE_INCFLAGS= -I$(IPE)/include 
IPE_LDFLAGS= -L $(IPE)/lib -lipe
SMS_INCFLAGS= -I$(SMS)/include
SMS_LDFLAGS= -L$(SMS)/lib -lsms
GPTL_INCFLAGS = -I/contrib/gptl/gptl-v5.4.4_impi_noomp/include
GPTL_LDFLAGS  = -L/contrib/gptl/gptl-v5.4.4_impi_noomp/lib -lgptl

################################################################################
################################################################################

.SUFFIXES: .f90 .F90 .f .c .C

%.o : %.f
	$(ESMF_F90COMPILER) -c $(ESMF_F90COMPILEOPTS) $(ESMF_F90COMPILEPATHS) $(ESMF_F90COMPILEFREENOCPP)   $(IPE_INCFLAGS) $(SMS_INCFLAGS) $(GPTL_INCFLAGS) $<
%.o : %.f90
	$(ESMF_F90COMPILER) -c $(ESMF_F90COMPILEOPTS) $(ESMF_F90COMPILEPATHS) $(ESMF_F90COMPILEFREENOCPP)   $(IPE_INCFLAGS) $(SMS_INCFLAGS) $(GPTL_INCFLAGS) $<
%.o : %.F90
	$(ESMF_F90COMPILER) -c $(ESMF_F90COMPILEOPTS) $(ESMF_F90COMPILEPATHS) $(ESMF_F90COMPILEFREECPP) $(ESMF_F90COMPILECPPFLAGS) -DESMF_VERSION_MAJOR=$(ESMF_VERSION_MAJOR)  $(IPE_INCFLAGS) $(SMS_INCFLAGS) $(GPTL_INCFLAGS) $<
        
%.o : %.c
	$(ESMF_CXXCOMPILER) -c $(ESMF_CXXCOMPILEOPTS) $(ESMF_CXXCOMPILEPATHSLOCAL) $(ESMF_CXXCOMPILEPATHS) $(ESMF_CXXCOMPILECPPFLAGS)  $(IPE_INCFLAGS) $(SMS_INCFLAGS) $(GPTL_INCFLAGS) $<

%.o : %.C
	$(ESMF_CXXCOMPILER) -c $(ESMF_CXXCOMPILEOPTS) $(ESMF_CXXCOMPILEPATHSLOCAL) $(ESMF_CXXCOMPILEPATHS) $(ESMF_CXXCOMPILECPPFLAGS)  $(IPE_INCFLAGS) $(SMS_INCFLAGS) $(GPTL_INCFLAGS) $<


# -----------------------------------------------------------------------------
#include Objects
# module dependencies:
# include Depends
mainApp.o: driver_v7.o
driver_v7.o: wamCap_v7.o ipeCap.o
wamCap_v7.o: module_sub_myWAM_Init.o module_update_WAM.o module_finalize_WAM.o
ipeCap.o: module_myIPE_Init.o  module_sub_myIPE_Init.o module_sub_myIPE_Init.o module_update_IPE.o module_finalize_IPE.o
module_sub_myWAM_Init.o: module_myWAM_Init.o
module_sub_myIPE_Init.o: module_myIPE_Init.o


# --------------------------------------------------------------------------
# -----------------------------------------------------------------------------
.PHONY: dust clean distclean info edit
dust:
	rm -f PET*.ESMF_LogFile *.nc
ipeclean:
	rm -f mainApp *.o *.mod
distclean: dust clean

info:
	@echo ==================================================================
	@echo ESMFMKFILE=$(ESMFMKFILE)
	@echo ==================================================================
	@cat $(ESMFMKFILE)
	@echo ==================================================================

edit:
	nedit mainApp.F90 driver_v7.F90 wamCap_v7.F90 ipeCap.F90 \
module_myIPE_Init.F90 module_sub_myIPE_Init.F90 module_update_IPE.F90 module_finalize_IPE.F90 \
module_myWAM_Init.F90 module_sub_myWAM_Init.F90 module_update_WAM.F90 module_finalize_WAM.F90 &
# DO NOT DELETE
#
# ------------------------------------------------------------------------------
#  # --- NUOPC additions ----------------------------------------------------------

#%.o : %.F90
#	$(ESMF_F90COMPILER) -c $(DEP_FRONTS) $(DEP_INCS) $(ESMF_F90COMPILEOPTS) $(ESMF_F90COMPILEPATHS) $(ESMF_F90COMPILEFREECPP) $(ESMF_F90COMPILECPPFLAGS) $<

.PRECIOUS: %.o

ipe.mk : ipeCap.o
	@echo "# ESMF self-describing build dependency makefile fragment" > $@
	@echo >> $@
	@echo "ESMF_DEP_FRONT     = ipeCap"                 >> $@
	@echo "ESMF_DEP_INCPATH   = `pwd`"                  >> $@
	@echo "ESMF_DEP_CMPL_OBJS = `pwd`/"$<               >> $@
	@echo "ESMF_DEP_LINK_OBJS = `pwd`/libipe_nuopc.a `pwd`/libipe.a"   >> $@
	@echo "ESMF_DEP_SHRD_PATH = $(SMS)/lib"         >> $@
	@echo "ESMF_DEP_SHRD_LIBS = sms"                >> $@
	@echo
	@echo "Finished generating ESMF self-describing build dependency makefile fragment:" $@
	@echo

nuopc: ipe.mk libipe_nuopc.a

libipe_nuopc.a: ipeCap.o
	ar cr $@ *.o

.PHONY: clean nuopcinstall
clean:
	rm -f *.o *.mod *.a *.mk

PWDIR := `pwd`

ifndef DESTDIR
DESTDIR := $(PWDIR)
endif

INSTDATE := $(shell date '+%Y-%m-%d-%H-%M-%S')
ifndef INSTDIR
INSTDIR  := IPE_$(INSTDATE)
endif

nuopcinstall:
	@gmake nuopc
	@mkdir -p $(DESTDIR)/$(INSTDIR)
	@cp libipe_nuopc.a ipecap.mod $(DESTDIR)/$(INSTDIR)
	@cp $(IPE)/lib/libipe.a $(DESTDIR)/$(INSTDIR)
	@sed -e 's;'$(PWDIR)';'$(DESTDIR)/$(INSTDIR)';g' ipe.mk > $(DESTDIR)/$(INSTDIR)/ipe.mk
	@touch VERSION
	@if [ -d .svn ]; then \
	  echo "SVN Repository" > VERSION; \
	  svn info . | grep URL >> VERSION; \
	  svn info . | grep "Last Changed Rev" >> VERSION; \
	fi
	@cp VERSION $(DESTDIR)/$(INSTDIR)/
	@echo Installation into \"$(DESTDIR)/$(INSTDIR)\" complete!
	@echo
# ------------------------------------------------------------------------------
#

clean:
	$(RM) libflip.a $(LIBOBJS) *.mod

OPT=-O0 -g -traceback

LIBOBJS = CMINOR.o \
          CTIPE-int.o \
          ELECXS.o \
          FLIP_GRID.o \
          INIT-PROFILES.o \
          KEMPRN.o \
          MINORA.o \
          Neut_Heating.o \
          Photoel-Freqs.o \
          RSDENA_EXB.o \
          RSJACA.o \
          RSLPSD.o \
          RSLPST.o \
          RSPE2B.o \
          RSPRIM.o \
          RSTEMC_EXB.o \
          Rates.o

%.o: %.f
	$(FC) $(OPT) -c $<

libflip.a: $(LIBOBJS)
	$(AR) $(ARFLAGS) $@ $(LIBOBJS)


CMINOR.o : CMINOR.f CTIPE-int.o  
CTIPE-int.o : CTIPE-int.f  
ELECXS.o : ELECXS.f  
FLIP_GRID.o : FLIP_GRID.f  
INIT-PROFILES.o : INIT-PROFILES.f CTIPE-int.o  
KEMPRN.o : KEMPRN.f CTIPE-int.o  
MINORA.o : MINORA.f CTIPE-int.o  
Neut_Heating.o : Neut_Heating.f CTIPE-int.o  
Photoel-Freqs.o : Photoel-Freqs.f  
RSDENA_EXB.o : RSDENA_EXB.f CTIPE-int.o  
RSJACA.o : RSJACA.f  
RSLPSD.o : RSLPSD.f CTIPE-int.o  
RSLPST.o : RSLPST.f CTIPE-int.o  
RSPE2B.o : RSPE2B.f CTIPE-int.o  
RSPRIM.o : RSPRIM.f CTIPE-int.o  
RSTEMC_EXB.o : RSTEMC_EXB.f CTIPE-int.o  
Rates.o : Rates.f  
#initialize_module_parameters.o : initialize_module_parameters.f CTIPE-int.o  
#module_tirosdata_ipe.o : module_tirosdata_ipe.f  
#sub_ionize_ipe.o : sub_ionize_ipe.f  
#sub_tiros_init_ipe.o : sub_tiros_init_ipe.f module_tirosdata_ipe.o  
#sub_tiros_ionize_ipe.o : sub_tiros_ionize_ipe.f module_tirosdata_ipe.o  


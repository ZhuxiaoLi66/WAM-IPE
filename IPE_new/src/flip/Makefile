clean:
	$(RM) libflip.a $(LIBOBJS)

LIBOBJS = CTIPE-int.o INIT-PROFILES.o Neut_Heating.o Rates.o

%.o: %.f
	$(FC) $(OPT) -c $<

libflip.a: $(LIBOBJS)
	$(AR) $(ARFLAGS) $@ $(LIBOBJS)


CTIPE-int.o: CTIPE-int.f Rates.o
INIT-PROFILES.o: INIT-PROFILES.f CTIPE-int.o
Neut_Heating.o: Neut_Heating.f CTIPE-int.o
Rates.o: Rates.f

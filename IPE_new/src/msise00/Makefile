clean:
	$(RM) libnrlmsis.a $(LIBOBJS)

LIBOBJS = physics_constants.o utils_spline.o physics_msis.o utils_constants.o hwm93.o

%.o: %.f90
	$(FC) $(OPT) -c $<

%.o: %.F90
	$(FC) $(CPPFLAGS) $(OPT) -c $<

libnrlmsis.a: $(LIBOBJS)
	$(AR) $(ARFLAGS) $@ $(LIBOBJS)

utils_constants.o: utils_constants.f90
physics_constants.o: physics_constants.f90 utils_constants.o
utils_spline.o: utils_spline.f90 utils_constants.o
physics_msis.o: physics_msis.f90 physics_constants.o utils_constants.o
hwm93.o: hwm93.f90

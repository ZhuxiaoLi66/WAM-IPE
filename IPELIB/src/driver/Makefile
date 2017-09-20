# driver Makefile

TOP = ..

PSRCS = driver_ipe.f90

OBJS   = $(wildcard $(addprefix ../,$(addsuffix /*.o,$(filter-out driver,$(DIRS)))))

.PHONY: exe

exe: $(EXE)

include $(TOP)/Makefile.common

$(EXE): driver_ipe.o $(OBJS)
ifeq ($(PARALLELISM),parallel)
	$(MODCMD) && $(FCP) $(FFLAGS) $(INCS) -o $@ $^ $(SMS_INCFLAGS) $(SMS_LDFLAGS)
else
	$(MODCMD) && $(FCP) $(FFLAGS) $(INCS) -o $@ $^ $(GPTL_LDFLAGS)
endif

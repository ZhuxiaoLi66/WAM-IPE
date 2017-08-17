# driver Makefile

TOP = ..

PSRCS = driver_ipe.f90

MODS   = $(wildcard $(addprefix ../,$(addsuffix /*.mod,$(filter-out driver,$(DIRS)))))
OBJS   = $(wildcard $(addprefix ../,$(addsuffix /*.o,  $(filter-out driver,$(DIRS)))))

.PHONY: all includes

all: $(EXE) includes

include $(TOP)/Makefile.common

$(EXE): driver_ipe.o $(LIB)
ifeq ($(PARALLELISM),parallel)
	$(MODCMD) && $(FCP) $(FFLAGS) $(INCS) -o $@ $^ -L$(dir $(LIB)) -lipe $(SMS_INCFLAGS) $(SMS_LDFLAGS)
else
	$(MODCMD) && $(FCS) $(FFLAGS) $(INCS) -o $@ $^ -L$(dir $(LIB)) -lipe 
endif

$(LIB): $(OBJS)
	ar ruv $@ $(OBJS)

includes: $(MODS)
	cp -uv $^ $(INCDIR)

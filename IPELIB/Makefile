TARGETS = clean

.PHONY: all clean

all:
	$(info Valid targets are:)
	$(foreach x,$(TARGETS),$(info - $(x)))
	$(error Invalid target specified)

clean:
	$(RM) -v bin/*
	$(RM) -v lib/*
	make -C src clean
	(cd tst && ./ipets clean)

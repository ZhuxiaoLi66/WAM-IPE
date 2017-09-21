# dummygptl Makefile

TOP = ..

ARCV = libgptl.a
SRCS = $(wildcard *.F90)
OBJS = $(addsuffix .o, $(basename $(SRCS)))

all: $(ARCV)

$(ARCV): $(OBJS)
	ar rc $@ $^

%.o: %.F90
	$(FCS) -c $<

clean:
	$(RM) *.o *.a

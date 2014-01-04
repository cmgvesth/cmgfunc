CCOPT = -xO4

all : pseg pmerge

pseg : pseg.c lnfac.h pgenwin.h pgenwin.o
	cc $(CCOPT) -o pseg pseg.c pgenwin.o -lm

pmerge : pmerge.c pgenwin.h pgenwin.o
	cc $(CCOPT) -o pmerge pmerge.c pgenwin.o -lm

pgenwin.o : pgenwin.c pgenwin.h
	cc $(CCOPT) -c pgenwin.c

clean:
	rm -f pseg pseg.o pgenwin.o

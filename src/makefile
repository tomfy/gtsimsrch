IDIR =../include
CC=gcc
CFLAGS=-I$(IDIR)

ODIR=obj
LDIR =../lib

LIBS=-lm -lpthread

_DEPS = vect.h gtset.h pedigree.h various.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ1 = vect.o duplicatesearch.o gtset.o pedigree.o various.o
OBJ1 = $(patsubst %,$(ODIR)/%,$(_OBJ1))

_OBJ2 = vect.o pedigree_test.o gtset.o pedigree.o
OBJ2 = $(patsubst %,$(ODIR)/%,$(_OBJ2))

_OBJ3 = vect.o  gtset.o pedigree.o find_parents.o
OBJ3 = $(patsubst %,$(ODIR)/%,$(_OBJ3))

_OBJ4 = vect.o  vcf_to_gts.o various.o
OBJ4 = $(patsubst %,$(ODIR)/%,$(_OBJ4))

_OBJ5 = plnkout_to_dsout.o 
OBJ5 = $(patsubst %,$(ODIR)/%,$(_OBJ5))


$(ODIR)/%.o: %.c $(DEPS)
	$(CC) -c -O3 -o $@ $< $(CFLAGS)

duplicatesearch: $(OBJ1)
	$(CC) -O3 -o $@ $^ $(CFLAGS) $(LIBS)

pedigree_test: $(OBJ2)
	$(CC) -O3 -o $@ $^ $(CFLAGS) $(LIBS)

find_parents: $(OBJ3)
	$(CC) -O3 -o $@ $^ $(CFLAGS) $(LIBS)

vcf_to_gts: $(OBJ4)
	$(CC) -O3 -o $@ $^ $(CFLAGS) $(LIBS)

plnkout_to_dsout: $(OBJ5)
	$(CC) -O3 -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: all

all: duplicatesearch pedigree_test find_parents vcf_to_gts plnkout_to_dsout

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~ 

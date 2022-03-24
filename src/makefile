IDIR =../include
CC=gcc
CFLAGS=-I$(IDIR)

ODIR=obj
LDIR =../lib

LIBS=-lm

_DEPS = vect.h gtset.h pedigree.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ1 = vect.o simsearch.o gtset.o pedigree.o various.o
OBJ1 = $(patsubst %,$(ODIR)/%,$(_OBJ1))

_OBJ2 = vect.o pedigree_test.o gtset.o pedigree.o
OBJ2 = $(patsubst %,$(ODIR)/%,$(_OBJ2))


$(ODIR)/%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

simsearch: $(OBJ1)
	$(CC) -O3 -o $@ $^ $(CFLAGS) $(LIBS)

pedigree_test: $(OBJ2)
	$(CC) -O3 -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: all

all: simsearch pedigree_test

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~ 

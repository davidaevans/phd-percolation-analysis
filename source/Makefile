PROG =	sc

SRCS = load.c dump.c geometry.c keywords.c main.c overlap.c percolate.c spanning.c read.c recipes.c simulate.c calculate.c rdf.c

OBJS = load.o dump.o geometry.o keywords.o main.o overlap.o percolate.o spanning.o read.o recipes.o simulate.o calculate.o rdf.o

LIBS = -lm
CC = cc
CFLAGS = -O2 -Wall
LDFLAGS = -s

all: $(PROG)

$(PROG): $(OBJS)
	$(CC) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

clean:
	rm -f $(PROG) $(OBJS)

load.o: global.h prototypes.h
dump.o: global.h prototypes.h
geometry.o: global.h prototypes.h
main.o: global.h prototypes.h
overlap.o: global.h prototypes.h
percolate.o: global.h prototypes.h
spanning.o: global.h prototypes.h
read.o: global.h prototypes.h
recipes.o: global.h prototypes.h
simulate.o: global.h prototypes.h
calculate.o: global.h prototypes.h
rdf.o: global.h prototypes.h
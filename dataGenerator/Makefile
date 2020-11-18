ifndef USEJEMALLOC
JEMALLOC =
else
JEMALLOCLD = $(shell jemalloc-config --libdir)
JEMALLOC = -DUSEJEMALLOC -L$(JEMALLOCLD) -ljemalloc
endif

#CFLAGS = -mcx16 -O3 -std=c++14 -march=native -Wall -Wextra
CFLAGS = -mcx16 -O3 -std=c++14 -march=native

OMPFLAGS = -DOPENMP -fopenmp
CILKFLAGS = -DCILK -fcilkplus
HGFLAGS = -DHOMEGROWN -pthread

ifdef CLANG
CC = clang++-5.0
PFLAGS = $(CILKFLAGS)
else ifdef OPENCILK
CC = ~/clang
PFLAGS = -DUSEMALLOC -DOPENCILK -fopencilk -lstdc++ -lm
#PFLAGS = -DOPENCILK -fopencilk -lstdc++ -lm
else ifdef CILK
CC = g++
PFLAGS = $(CILKFLAGS)
else ifdef OPENMP
CC = g++
PFLAGS = $(OMPFLAGS)
else ifdef HOMEGROWN
CC = g++
PFLAGS = $(HGFLAGS)
else ifdef SERIAL
CC = g++
PFLAGS =
else # default is sequential
CC = g++
PFLAGS =
endif

#########

GENERATORS = uniform

.PHONY: all clean
all: $(GENERATORS)

$(COMMON) :
	ln -s ../../common/$@ .

%.o : %.C $(COMMON)
	$(CC) $(PFLAGS) -c $< -o $@

uniform : uniform.C geometryData.h $(COMMON)
	$(CC) $(PFLAGS) -o $@ $@.C
clean :
	rm -f *.o $(GENERATORS)

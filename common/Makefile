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

OBJS = bench.o
BNCHMRK = bench

TIME = $(BNCHMRK)Time

# Make the benchmark, and make the output check code in the common directory
all : $(BNCHMRK)

# Make all implementation objects
%.o : %.C
	$(CC) $(CFLAGS) $(PFLAGS) $(JEMALLOC) -c $< -o $@

# Make timing code object file
$(TIME).o : $(TIME).C
	$(CC) $(CFLAGS) $(PFLAGS) $(JEMALLOC) -c $< -o $@

# Make benchmark
$(BNCHMRK) : $(TIME).o $(OBJS)
	$(CC) $(CFLAGS) $(PFLAGS) -o $@ $(TIME).o $(OBJS) $(LIBS) $(JEMALLOC)

clean :
	rm -f $(BNCHMRK) *.o

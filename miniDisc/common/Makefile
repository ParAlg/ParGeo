ifndef USEJEMALLOC
JEMALLOC =
else
JEMALLOCLD = $(shell jemalloc-config --libdir)
JEMALLOC = -DUSEJEMALLOC -L$(JEMALLOCLD) -ljemalloc
endif

#CFLAGS = -mcx16 -O3 -std=c++14 -march=native -Wall -Wextra
CFLAGS = -mcx16 -O3 -std=c++17 -march=native

OMPFLAGS = -DOPENMP -fopenmp
CILKFLAGS = -DCILK -fcilkplus
HGFLAGS = -DHOMEGROWN -pthread

ifdef SILENT
OTHERFLAGS = -DSILENT
else
OTHERFLAGS =
endif

# ifdef CLANG
# CC = clang++-5.0
# PFLAGS = $(CILKFLAGS)
# else
ifdef OPENCILK
CC = ~/clang
PFLAGS = -DUSEMALLOC -DOPENCILK -fopencilk -lstdc++ -lm  $(OTHERFLAGS)
#PFLAGS = -DOPENCILK -fopencilk -lstdc++ -lm
else ifdef CILK
CC = g++
PFLAGS = $(CILKFLAGS) $(OTHERFLAGS)
# else ifdef OPENMP
# CC = g++
# PFLAGS = $(OMPFLAGS) $(OTHERFLAGS)
else ifdef HOMEGROWN
CC = g++
PFLAGS = $(HGFLAGS)  $(OTHERFLAGS)
# else ifdef SERIAL
# CC = g++
# PFLAGS = $(OTHERFLAGS)
else # default is sequential
CC = g++
PFLAGS = $(OTHERFLAGS)
endif

#########

OBJS = miniDisc.o
BNCHMRK = miniDisc

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

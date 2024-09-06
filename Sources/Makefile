# Compiler
CC = mpicc

# Flags for release and debug
OPTIMIZATION = -O3
DEBUG_FLAGS = -g
CFLAGS_RELEASE = $(OPTIMIZATION)
CFLAGS_DEBUG = $(DEBUG_FLAGS)

# Source files
SRC_SIM 	=       helper 		\
			parallel 	\
			mpifft		\
			calc		\
			sort		\
			prime		\
			rsa		\
		    main

# Default target (release)
all: release

# Release target
release: CFLAGS = $(CFLAGS_RELEASE)
release: rsa

# Debug target
debug: CFLAGS = $(CFLAGS_DEBUG)
debug: rsa

# Building the final executable
rsa: $(SRC_SIM:%=%.o)
	$(CC) $(CFLAGS) -o rsa $(SRC_SIM:%=%.o) -lm

# Rule for compiling object files
%.o : %.c
	$(CC) -c $(CFLAGS) $*.c -o $*.o

# Clean up the build
clean:
	/bin/rm -f $(SRC_SIM:%=%.o) rsa

# Dependencies
sort.o		: sort.h
helper.o     	: helper.h 
mpifft.o	: mpifft.h helper.h 
calc.o		: calc.h helper.h prime.h
prime.o		: prime.h helper.h sort.h
parallel.o 	: parallel.h
rsa.o		: rsa.h helper.h calc.h prime.h
main.o       	: helper.h parallel.h calc.h mpifft.h rsa.h prime.h

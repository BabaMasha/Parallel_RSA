CC = mpicc
CFLAGS = -O3
SRC_SIM 	=       helper 		\
			parallel 	\
			mpifft		\
			calc		\
			sort		\
			prime		\
			rsa		\
		    main


all:
	@make rsa

rsa: $(SRC_SIM:%=%.o)
	$(CC) $(CFLAGS) -o rsa $(SRC_SIM:%=%.o)  -lm

%.o : %.c
	$(CC) -c $(CFLAGS) $*.c -o $*.o

clean:
	/bin/rm -f $(SRC_SIM:%=%.o) rsa 

sort.o		: sort.h
helper.o     	: helper.h 
mpifft.o	: mpifft.h helper.h 
calc.o		: calc.h helper.h prime.h
prime.o		: prime.h helper.h sort.h
parallel.o 	: parallel.h
rsa.o		: rsa.h helper.h calc.h prime.h


#$mpifft_test.o	: helper.h mpifft.h 
$main.o       	: helper.h parallel.h calc.h mpifft.h rsa.h prime.h

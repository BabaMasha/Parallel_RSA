Parallel RSA with C/MPI by Maryna Babayeva
Utrecht University
Course in Scientific Computing
2011

#calc.c, calc.h: 
basic mathematical operations in parallel (add, mult, sub, div,...).

#helper.c, helper.h: 
helper functions as file IO, or printing to the terminal of the high-precision numbers. To store a high-precision number in a file a special notation will be used:
1.23*10^2 will be stored in a file as
--------------
2
0123
--------------

or
--------------
<exp>
<sign><significant>
--------------

#prime.c, prime.h: 
contains prime sieve, prime generator with primality test

#rsa.c, rsa.h: 
contains routines for RSA encryption/decryption and key generation, uses smallPrime.txt and data in KEYS folder


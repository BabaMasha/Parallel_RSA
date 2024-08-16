#include "helper.h"
#include "prime.h"
#include "mpi.h"
#include "sort.h"

//#define TEST FALSE
//#define TESTIT FALSE

//check if index is local to the processor
unsigned int isLocal(unsigned int globalIndex, int width, int p, int s)
{
	int ret = TRUE;
	if(p > 1)
		ret = (globalIndex/width)%p == s;
	return ret;
}
//return local index
unsigned int local(unsigned int globalIndex, int width, int p, int s)
{
	unsigned int ret = globalIndex;
	if(p > 1)
		ret = (globalIndex-width*s)/p+globalIndex%width;
	return ret;
}

//return global index
unsigned int global(unsigned int localIndex, int width, int p, int s)
{
	unsigned int ret = localIndex;
	if(p > 1)
		ret = (localIndex/width)*width*p+width*s+localIndex%width;
	return ret;
}

unsigned int getBlockSize(int width, int p, int s, unsigned int n)
{
	unsigned long add = 0;
	unsigned long rest = n%(width*p);
	if(p>1)
	{
		if (rest>width*s)
		{
			if (width*(s+1)<rest)
				add = width;
			else
				add = rest-width*s;
		}
		else
			add=0;
		
		return ceil(n+1-rest)/p+add;
	}else{
		return n;
	}
	
}

void initArray(unsigned int* primeArray, unsigned int blockSize)
{	
	unsigned int i, counter;
	counter=0;
	int procID;
	MPI_Comm_rank(MPI_COMM_WORLD, &procID);
	for(i = 0;  i<=blockSize; i=i+2)
	{
		counter++;
		primeArray[i] = FALSE;
		primeArray[i+1] = TRUE;
	}
#ifdef TESTIT
	printf("init iterations %d\n", counter); fflush(stdout);
#endif
	if(procID==0) primeArray[1] = FALSE;
	
}

unsigned int mark(unsigned int* localPrimeArray, unsigned int key, unsigned int maxNum, int numProcs, int procID, int width)
{
	unsigned int counter = 0;
	unsigned int i, k;
	
	i = key*key;
	
	for(i; i <= maxNum; i = i+2*key)
	{
		
		counter++;
		if(isLocal(i, width, numProcs, procID))
		{
			k = local(i, width, numProcs, procID);	
			break;
		}
		
	}
	
	unsigned int blockSize = getBlockSize(width, numProcs, procID, maxNum);
	
	//mark all multiples
	for(k; k <= blockSize; k = k+2*key)
	{
		counter++;
		localPrimeArray[k] = FALSE;
	}
	
	return counter;
}

void makeKeyList(unsigned int* keyArray, unsigned int maxNum)
{
	
    unsigned int i, j, k, t, stopAt, counter;	
	unsigned int *keys;
	int procID;
	MPI_Comm_rank(MPI_COMM_WORLD, &procID);
	keys = (unsigned int *)calloc(maxNum+1, sizeof(unsigned int));
	
	stopAt = ceil(sqrt(maxNum));
	counter=0;
	
	for(i = 0; i <= maxNum; i = i + 2)
	{
		keys[i]=FALSE;
		keys[i+1]=TRUE;
		counter++;
	}
	
	// 1 is not a prime
	keys[1] = FALSE;
	
	for (k = 3; k <= stopAt; k++) 
	{
        if (keys[k]==TRUE) 
		{
			j=0;
			for(i = k*(k+j), j=2; i <= maxNum; i=k*(k+j), j=2*j)
			{
				keys[i] = FALSE; 
				counter++;
			}
			
		}
	}
	
	t=0;
	for(i=3; i <= maxNum; i++)
	{

		if(keys[i])
		{
			keyArray[t] = i;
#ifdef TEST
			if(procID == 0)
				printf("keyArray[%d]: %d\n", t, keyArray[t]); fflush(stdout);
#endif
			t++;
		}
	}
#ifdef TESTIT
	printf("key iterations %d\n", counter); fflush(stdout);
#endif
	free(keys);
}

void makePrimeList(unsigned int maxNum, 
				   int numProcs,
				   int procID,
				   unsigned int* globalPrimeArray, 
				   unsigned int width)
{
    unsigned int i, j, k, t, stopAt, counter;
	unsigned int *localPrimeArray, *keys;
	unsigned int *offset, *blockSize2;
	double runtime;
	unsigned int blockSize;
	
	stopAt = ceil(sqrt(maxNum));
	
	//the sieve algorithm
	blockSize = getBlockSize(width, numProcs, procID, maxNum);
	
	//allocate memory
	localPrimeArray = (unsigned int *)calloc(blockSize, sizeof(unsigned int));
	keys = (unsigned int *)calloc(stopAt, sizeof(unsigned int));
	offset = (unsigned  int *)calloc(numProcs, sizeof(unsigned  int));

	MPI_Barrier(MPI_COMM_WORLD);
	runtime = -MPI_Wtime();
	
	initArray(localPrimeArray, blockSize);
	//generate a list with all prime key factors up to sqrt(n)
	makeKeyList(keys, stopAt); 
	counter = 0;
	
	//Filter primes
    t=0;	
	do
	{	
		//filter non-primes
		counter += mark(localPrimeArray, keys[t], maxNum, numProcs, procID, width);
		t++;
		
		//stop when sqrt(n) is reached
	}while((keys[t] <= stopAt) && (keys[t] != 0));
	
	//set 2 as prime
	if(isLocal(2, width, numProcs, procID)) 
		localPrimeArray[local(2, width, numProcs, procID)] = TRUE; 
		
	MPI_Barrier(MPI_COMM_WORLD);
	runtime += MPI_Wtime();
	
#ifdef TESTIT
	printf("id %d: make iterations %d\n", procID, counter); fflush(stdout);
#endif

	unsigned int *numPrimes = (unsigned  int *)calloc(numProcs, sizeof(unsigned  int));
	unsigned int Primes = 0;
	for(i=0,j=0;i<blockSize; i++)
	{
		if(localPrimeArray[i]>0)
		{
			localPrimeArray[j] = global(i, width, numProcs, procID);
			numPrimes[procID]++;
			j++;
			if(i>1)
			   localPrimeArray[i]=0;
			
		}
	}

	MPI_Allreduce(&numPrimes[procID], &Primes,1,MPI_UNSIGNED,MPI_SUM,MPI_COMM_WORLD);
	MPI_Allgather(&numPrimes[procID], 1, MPI_UNSIGNED, numPrimes, 1, MPI_UNSIGNED, MPI_COMM_WORLD );
	
	//combine all primes in global array
	offset[0] = 0;
	for(i = 1; i < numProcs; i++){
			offset[i] = offset[i-1]+numPrimes[i-1];
	}
	
	MPI_Allgatherv(localPrimeArray, numPrimes[procID], MPI_UNSIGNED, 
				   globalPrimeArray, (int*)numPrimes, (int*)offset, 
				   MPI_UNSIGNED, MPI_COMM_WORLD);

	int m = log_2(numProcs);
	para_QuickSort((int *)globalPrimeArray,0,Primes-1,m,0,procID);


	MPI_Bcast(globalPrimeArray, Primes, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
	
	if(procID == 0) 
	{
		printf("Primes=%d\n", Primes);
	}
	
#ifdef TEST	
	if(procID == 0) 
	{
		printf("In %f seconds we found %d primes less than or equal to %d.\n",
			   runtime, numPrimes[procID-1], maxNum-1);
	}
#endif
	
	free(localPrimeArray);
	free(keys);
	free(offset);
	free(numPrimes);
}

// create a random prime number of desired length
void paraGenPrime(IEEE_754_FloatNum *prime, 
				  UINT32 *smallPrimeList,
				  UINT32 numProcs, 
				  UINT32 procID,
				  UINT32 length)
{
	UINT32 width = BLOCK_SIZE(procID, numProcs, prime->precision);
	
	IEEE_754_FloatNum smallPrime;
	smallPrime.base	= BASE;
	smallPrime.exp	= 0;
	smallPrime.sign	= FALSE;
	smallPrime.significand = (UINT32 *)calloc(width, sizeof(UINT32));
	smallPrime.precision = PRECISION;
	memset(smallPrime.significand, 0, width*sizeof(UINT32));
	
	IEEE_754_FloatNum one;
	one.base	= BASE;
	one.exp	= 0;
	one.sign	= FALSE;
	one.significand = (UINT32 *)calloc(width, sizeof(UINT32));
	one.precision = PRECISION;
	memset(one.significand, 0, width*sizeof(UINT32));
	
	IEEE_754_FloatNum two;
	two.base	= BASE;
	two.exp	= 0;
	two.sign	= FALSE;
	two.significand = (UINT32 *)calloc(width, sizeof(UINT32));
	two.precision = PRECISION;
	memset(two.significand, 0, width*sizeof(UINT32));
	
	IEEE_754_FloatNum res;
	res.base	= BASE;
	res.exp	= 0;
	res.sign	= FALSE;
	res.significand = (UINT32 *)calloc(width, sizeof(UINT32));
	res.precision = PRECISION;
	memset(res.significand, 0, width*sizeof(UINT32));
	
	IEEE_754_FloatNum three;
	three.base	= BASE;
	three.exp	= 0;
	three.sign	= FALSE;
	three.significand = (UINT32 *)calloc(width, sizeof(UINT32));
	three.precision = PRECISION;
	memset(three.significand, 0, width*sizeof(UINT32));
		
	//three = 3;
	if(procID == 0)
	{
		one.significand[0]=1;
		two.significand[0]=2;
		three.significand[0]=3;
	}

	//get random number
	paraRand(prime, numProcs, procID, length);
	
	//make number odd if neccessary
	if(isODD(prime,numProcs, procID) != TRUE)
		paraAdd(prime, &three, prime, numProcs, procID, ADD_CALL);
		
	//get a small prime
	char *fileName;
	if(KEYLENGTH>32)
	{
		fileName = makeFileName("P", KEYLENGTH/2);
		read_IEEE_754_FloatNum(&smallPrime, numProcs, procID, fileName);
	}
	else
	{
		read_IEEE_754_FloatNum(&smallPrime, numProcs, procID, "smallPrime.txt");
	}

	print_IEEE_754_FloatNum(&smallPrime, "smallPrime", numProcs, procID);
	
	//check gcd(q, Π ) = 1
	paraGCD(prime, &smallPrime, &res, numProcs, procID);
	int j = 0;
	while(equal_IEEE(&res, &one, numProcs, procID) != TRUE)
	{
		j++;
		//printf("j=%d\n", j);
		//satisfying gcd(q, Π ) = 1? otherwise prime = prime + 2
		paraGCD(prime, &smallPrime, &res, numProcs, procID);
		print_IEEE_754_FloatNum(&res, "gcd", numProcs, procID);
		paraAdd(prime, &smallPrime, prime, numProcs, procID, ADD_CALL);
		paraAdd(prime, &two, prime, numProcs, procID, ADD_CALL);
	}

	print_IEEE_754_FloatNum(prime, "1. PRIME", numProcs, procID);
	int i = 0;
	while(isPrime(prime, smallPrimeList, numProcs, procID, length) != TRUE)
	{	
		i++;
		paraAdd(prime, &two, prime, numProcs, procID, ADD_CALL);
	}
	if(procID == 0)
		printf("nach %d versuchen\n", i); 
	
	free(one.significand);
	free(two.significand);
	free(three.significand);
	free(res.significand);
	free(smallPrime.significand);

}

BOOLEAN isPrime(IEEE_754_FloatNum *prime,
				UINT32 *smallPrimeList,
				UINT32 numProcs, 
				UINT32 procID,
				UINT32 length)
{
	UINT32 k = 0;
	BOOLEAN ret = TRUE;
	UINT32 width = BLOCK_SIZE(procID, numProcs, prime->precision);

	IEEE_754_FloatNum null;
	null.base	= BASE;
	null.exp	= 0;
	null.sign	= FALSE;
	null.significand = (UINT32 *)calloc(width, sizeof(UINT32));
	null.precision = PRECISION;
	memset(null.significand, 0, width*sizeof(UINT32));
	
	IEEE_754_FloatNum smallPrime;
	smallPrime.base	= BASE;
	smallPrime.exp	= 0;
	smallPrime.sign	= FALSE;
	smallPrime.significand = (UINT32 *)calloc(width, sizeof(UINT32));
	smallPrime.precision = PRECISION;
	memset(smallPrime.significand, 0, width*sizeof(UINT32));
	
	IEEE_754_FloatNum a;
	a.base	= BASE;
	a.exp	= 0;
	a.sign	= FALSE;
	a.significand = (UINT32 *)calloc(width, sizeof(UINT32));
	a.precision = PRECISION;
	
	IEEE_754_FloatNum n_minus_two;
	n_minus_two.base	= BASE;
	n_minus_two.exp	= 0;
	n_minus_two.sign	= FALSE;
	n_minus_two.significand = (UINT32 *)calloc(width, sizeof(UINT32));
	n_minus_two.precision = PRECISION;
	
	IEEE_754_FloatNum two;
	two.base	= BASE;
	two.exp	= 0;
	two.sign	= FALSE;
	two.significand = (UINT32 *)calloc(width, sizeof(UINT32));
	two.precision = PRECISION;
	memset(two.significand, 0, width*sizeof(UINT32));
	
	//two = 2
	if(procID == 0)
	{
		two.significand[0]=2;
	}
	
	//check whether multiple of known primes
	k = 0;
	while(smallPrimeList[k] > 0 && k<167)
	{
		toIEEE_754(smallPrimeList[k],&smallPrime,numProcs,procID);
		paraMod(prime, &smallPrime, &a, numProcs, procID);
		if(equal_IEEE(&a, &null, numProcs, procID) == TRUE)
		{
			ret = FALSE;
			goto exit;
		}
		k++;

	}

	paraSub(prime, &two, &n_minus_two, numProcs, procID);
		
	//repeat k times Miller-Rabin test
	k=0;
	while(k < 5)
	{
		k++;
		//pick a random integer a in the range [2, n − 2]
		do
		{
			paraRand(&a, numProcs, procID, length);
			//printf("JA");
		}while(smaller_IEEE(&a, &n_minus_two, numProcs, procID) != TRUE &&
			   smaller_IEEE(&a, &two, numProcs, procID) == TRUE) ;
	
		if( millerRabinTest(prime, &a, numProcs,procID) != TRUE)
		{
			ret = FALSE;
			goto exit;
		}
					
	}//end WHILE loop

exit:
	free(null.significand);
	free(a.significand);
	free(n_minus_two.significand);
	free(two.significand);
	free(smallPrime.significand);

	return ret;
	
}

//TRUE if prime
BOOLEAN millerRabinTest(IEEE_754_FloatNum *prime,
						IEEE_754_FloatNum *a,
						UINT32 numProcs, 
						UINT32 procID)
{
	UINT32 width = BLOCK_SIZE(procID, numProcs, prime->precision);
	BOOLEAN ret = FALSE;
	
	IEEE_754_FloatNum x;
	x.base	= BASE;
	x.exp	= 0;
	x.sign	= FALSE;
	x.significand = (UINT32 *)calloc(width, sizeof(UINT32));
	x.precision = PRECISION;
	
	IEEE_754_FloatNum d;
	d.base	= BASE;
	d.exp	= 0;
	d.sign	= FALSE;
	d.significand = (UINT32 *)calloc(width, sizeof(UINT32));
	d.precision = PRECISION;
	
	IEEE_754_FloatNum n_minus_one;
	n_minus_one.base	= BASE;
	n_minus_one.exp		= 0;
	n_minus_one.sign	= FALSE;
	n_minus_one.significand = (UINT32 *)calloc(width, sizeof(UINT32));
	n_minus_one.precision = PRECISION;
	
	IEEE_754_FloatNum one;
	one.base	= BASE;
	one.exp		= 0;
	one.sign	= FALSE;
	one.significand = (UINT32 *)calloc(width, sizeof(UINT32));
	one.precision = PRECISION;
	memset(one.significand, 0, width*sizeof(UINT32));
	
	//one = 1;
	if(procID == 0)
	{
		one.significand[0]=1;
	}
	
	UINT32 s, r;
	
	//1. Get s and d, s.t. n − 1 = 2^s·d with d odd
	paraSub(prime, &one, &n_minus_one, numProcs, procID);
	memCPY(&d, &n_minus_one, width);
	
	s=0;

	while(isODD(&d, numProcs, procID) != TRUE)
	{
		paraDivBy2(&d, &d, numProcs, procID);
		s++;
	}	
	
	//2. x = a^d mod n
	paraModExp(a, &d, prime, &x, numProcs, procID);
	
	//if x = 1 OR x = n-1
	if( equal_IEEE(&x, &one, numProcs,procID) == TRUE || equal_IEEE(&x, &n_minus_one, numProcs,procID) == TRUE)
	{
		ret = TRUE;
		goto exit;
	}
	else
	{
		for(r = 0; r < s; r++)
		{
			//3. x = x^2 mod n
			paraPow2(&x, &d, numProcs,procID);
			paraMod(&d, prime, &x, numProcs, procID);

			//if x = 1
			if (equal_IEEE(&x, &one, numProcs,procID) == TRUE)
			{
				ret = FALSE;
				goto exit;
			}
			
			//if x = n-1
			if (equal_IEEE(&x, &n_minus_one, numProcs,procID) == TRUE)
			{
				ret = TRUE;
				goto exit;
			}
			
		}//end FOR loop
		
		ret = FALSE;

	}//end ELSE
	

	
exit:	
	free(n_minus_one.significand);
	free(one.significand);
	free(d.significand);
	free(x.significand);
	return ret;
}



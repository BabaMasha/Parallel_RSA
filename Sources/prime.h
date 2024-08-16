
unsigned int isLocal(unsigned int globalIndex, int width, int p, int s);

unsigned int local(unsigned int globalIndex, int width, int p, int s);

unsigned int global(unsigned int localIndex, int width, int p, int s);

unsigned int getBlockSize(int width, int p, int s,unsigned int n);

void initArray(unsigned int* primeArray, unsigned int blockSize);

unsigned int mark(unsigned int* localPrimeArray, unsigned int key, unsigned int maxNum, int numProcs, int procID, int width);

void makeKeyList(unsigned int* keyArray, unsigned int maxNum);

void makePrimeList(unsigned int maxNum, 
				   int numProcs,
				   int procID,
				   unsigned int* globalPrimeArray, 
				   unsigned int width);
void paraGenPrime(IEEE_754_FloatNum *prime,
				  UINT32 *smallPrimeList,
				  UINT32 numProcs, 
				  UINT32 procID,
				  UINT32 length);
BOOLEAN isPrime(IEEE_754_FloatNum *prime,
				UINT32 *smallPrimeList,
				UINT32 numProcs, 
				UINT32 procID,
				UINT32 length);
BOOLEAN millerRabinTest(IEEE_754_FloatNum *prime,
						IEEE_754_FloatNum *a,
						UINT32 numProcs, 
						UINT32 procID);
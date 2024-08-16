#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include "helper.h"
#include "parallel.h"
#include "calc.h"
#include "prime.h"
#include "rsa.h"

void RSA_genKey(UINT32 procID, UINT32 numProcs, UINT32 keylength)
{
	UINT32 width = BLOCK_SIZE(procID, numProcs, PRECISION);

	IEEE_754_FloatNum E;
	E.base	= BASE;
	E.exp	= 0;
	E.sign	= FALSE;
	E.significand = (UINT32 *)calloc(width, sizeof(UINT32));
	E.precision = PRECISION;
	memset(E.significand, 0, width*sizeof(UINT32));
	
	IEEE_754_FloatNum P;
	P.base	= BASE;
	P.exp	= 0;
	P.sign	= FALSE;
	P.significand = (UINT32 *)calloc(width, sizeof(UINT32));
	P.precision = PRECISION;
	memset(P.significand, 0, width*sizeof(UINT32));
	
	IEEE_754_FloatNum Q;
	Q.base	= BASE;
	Q.exp	= 0;
	Q.sign	= FALSE;
	Q.significand = (UINT32 *)calloc(width, sizeof(UINT32));
	Q.precision = PRECISION;
	memset(Q.significand, 0, width*sizeof(UINT32));
	
	IEEE_754_FloatNum P1;
	P1.base	= BASE;
	P1.exp	= 0;
	P1.sign	= FALSE;
	P1.significand = (UINT32 *)calloc(width, sizeof(UINT32));
	P1.precision = PRECISION;
	memset(P1.significand, 0, width*sizeof(UINT32));
	
	IEEE_754_FloatNum Q1;
	Q1.base	= BASE;
	Q1.exp	= 0;
	Q1.sign	= FALSE;
	Q1.significand = (UINT32 *)calloc(width, sizeof(UINT32));
	Q1.precision = PRECISION;
	memset(Q1.significand, 0, width*sizeof(UINT32));
	
	IEEE_754_FloatNum H;
	H.base	= BASE;
	H.exp	= 0;
	H.sign	= FALSE;
	H.significand = (UINT32 *)calloc(width, sizeof(UINT32));
	H.precision = PRECISION;
	memset(H.significand, 0, width*sizeof(UINT32));
	
	IEEE_754_FloatNum G;
	G.base	= BASE;
	G.exp	= 0;
	G.sign	= FALSE;
	G.significand = (UINT32 *)calloc(width, sizeof(UINT32));
	G.precision = PRECISION;
	memset(G.significand, 0, width*sizeof(UINT32));
	
	IEEE_754_FloatNum D;
	D.base	= BASE;
	D.exp	= 0;
	D.sign	= FALSE;
	D.significand = (UINT32 *)calloc(width, sizeof(UINT32));
	D.precision = PRECISION;
	memset(D.significand, 0, width*sizeof(UINT32));
	
	IEEE_754_FloatNum DP;
	DP.base	= BASE;
	DP.exp	= 0;
	DP.sign	= FALSE;
	DP.significand = (UINT32 *)calloc(width, sizeof(UINT32));
	DP.precision = PRECISION;
	memset(DP.significand, 0, width*sizeof(UINT32));
	
	IEEE_754_FloatNum DQ;
	DQ.base	= BASE;
	DQ.exp	= 0;
	DQ.sign	= FALSE;
	DQ.significand = (UINT32 *)calloc(width, sizeof(UINT32));
	DQ.precision = PRECISION;
	memset(DQ.significand, 0, width*sizeof(UINT32));
	
	IEEE_754_FloatNum QP;
	QP.base	= BASE;
	QP.exp	= 0;
	QP.sign	= FALSE;
	QP.significand = (UINT32 *)calloc(width, sizeof(UINT32));
	QP.precision = PRECISION;
	memset(QP.significand, 0, width*sizeof(UINT32));
	
	IEEE_754_FloatNum one;
	one.base	= BASE;
	one.exp		= 0;
	one.sign	= FALSE;
	one.significand = (UINT32 *)calloc(width, sizeof(UINT32));
	one.precision = PRECISION;
	memset(one.significand, 0, width*sizeof(UINT32));
	
	IEEE_754_FloatNum N;
	N.base	= BASE;
	N.exp		= 0;
	N.sign	= FALSE;
	N.significand = (UINT32 *)calloc(width, sizeof(UINT32));
	N.precision = PRECISION;
	memset(N.significand, 0, width*sizeof(UINT32));
	
	if(procID == 0)
	{
		one.significand[0]=1;
	}
	
	char *fileName;
	
	//make a list of small primes
	UINT32 maxNum = (UINT32)pow(10,3);
	unsigned int *globalPrimeArray = (unsigned int *)calloc(maxNum+1, sizeof(unsigned int));
	makePrimeList(maxNum, numProcs, procID, globalPrimeArray, 2);
	
	read_IEEE_754_FloatNum(&E, numProcs, procID, "KEYS/E.txt");
	print_IEEE_754_FloatNum(&E, "E", numProcs, procID);

	int i = 1;
	do
	{
		if(procID==0)
			printf("%d try...\n", i);
		
		//generate random primes
//		paraGenPrime(&P, globalPrimeArray, numProcs, procID, KEYLENGTH);
//		print_IEEE_754_FloatNum(&P, "P", numProcs, procID);
//
//		paraGenPrime(&Q, globalPrimeArray, numProcs, procID, KEYLENGTH);
//		print_IEEE_754_FloatNum(&Q, "Q", numProcs, procID);

		//read alredy pregenerated prime from file
		fileName = makeFileName("P", KEYLENGTH);
		read_IEEE_754_FloatNum(&P, numProcs, procID, fileName);
		print_IEEE_754_FloatNum(&P, ">>P", numProcs, procID); 
		
		fileName = makeFileName("Q", KEYLENGTH);
		read_IEEE_754_FloatNum(&Q, numProcs, procID, fileName);
		print_IEEE_754_FloatNum(&Q, ">>Q", numProcs, procID); 
		
		
		paraSub(&Q, &one, &Q1, numProcs, procID);
		print_IEEE_754_FloatNum(&Q1, "Q1", numProcs, procID);

		paraSub(&P, &one, &P1, numProcs, procID);
		print_IEEE_754_FloatNum(&P1, "P1", numProcs, procID);

		
		paraMult(&Q1, &P1, &H, numProcs, procID);
		print_IEEE_754_FloatNum(&H, "H", numProcs, procID);

		
		paraGCD(&E, &H, &G, numProcs, procID);
		print_IEEE_754_FloatNum(&G, "G", numProcs, procID);


		i++;
	}
	while(equal_IEEE(&G, &one, numProcs, procID) != TRUE);
	if(procID == 0)
		printf("gen key versuche: %d\n", i-1);

	if(procID == 0)
		printf("Now calc the others...\n");
	
	paraInvMod(&E, &H, &D, numProcs, procID);//D  = E^-1 mod ((P-1)*(Q-1))
	paraMod(&D, &P1, &DP, numProcs, procID);//DP = D mod (P - 1)
	paraMod(&D, &Q1, &DQ, numProcs, procID);//DQ = D mod (Q - 1)
	paraInvMod(&Q, &P, &QP, numProcs, procID);//QP = Q^-1 mod P
	paraMult(&Q, &P, &N, numProcs, procID);
	
	//store primes to file
	fileName = makeFileName("P", KEYLENGTH);
	write_IEEE_754_FloatNum(&P, numProcs, procID, fileName);
	fileName = makeFileName("Q", KEYLENGTH);
	write_IEEE_754_FloatNum(&Q, numProcs, procID, fileName);
	
	fileName = makeFileName("D", KEYLENGTH);
	write_IEEE_754_FloatNum(&D, numProcs, procID, fileName);
	fileName = makeFileName("N", KEYLENGTH);
	write_IEEE_754_FloatNum(&N, numProcs, procID, fileName);
	
	fileName = makeFileName("DP", KEYLENGTH);
	write_IEEE_754_FloatNum(&DP, numProcs, procID, fileName);
	fileName = makeFileName("DQ", KEYLENGTH);
	write_IEEE_754_FloatNum(&DQ, numProcs, procID, fileName);
	fileName = makeFileName("QP", KEYLENGTH);
	write_IEEE_754_FloatNum(&QP, numProcs, procID, fileName);

	print_IEEE_754_FloatNum(&D, "D", numProcs, procID);
	print_IEEE_754_FloatNum(&DP, "DP", numProcs, procID);
	print_IEEE_754_FloatNum(&DQ, "DQ", numProcs, procID);
	print_IEEE_754_FloatNum(&QP, "QP", numProcs, procID);
	print_IEEE_754_FloatNum(&N, "N", numProcs, procID);

	free(fileName);
	free(N.significand);
	free(Q1.significand);
	free(P1.significand);
	free(Q.significand);
	free(P.significand);
	free(one.significand);
	free(H.significand);
	free(E.significand);
	free(D.significand);
	free(DP.significand);
	free(DQ.significand);
	free(QP.significand);
	free(G.significand);
	free(globalPrimeArray);	
}

//m=c^D (mod N) with CRT
void RSA_decrypt(IEEE_754_FloatNum *c, 
				 IEEE_754_FloatNum *m,  
				 IEEE_754_FloatNum *Q,
				 IEEE_754_FloatNum *P,
				 IEEE_754_FloatNum *QP,
				 IEEE_754_FloatNum *DP,
				 IEEE_754_FloatNum *DQ,
				 UINT32 procID, 
				 UINT32 numProcs )
{
	UINT32 width = BLOCK_SIZE(procID, numProcs, c->precision);

	IEEE_754_FloatNum T;
	T.base	= BASE;
	T.exp		= 0;
	T.sign	= FALSE;
	T.significand = (UINT32 *)calloc(width, sizeof(UINT32));
	T.precision = c->precision;
	memset(T.significand, 0, width*sizeof(UINT32));
	
	IEEE_754_FloatNum T1;
	T1.base	= BASE;
	T1.exp		= 0;
	T1.sign	= FALSE;
	T1.significand = (UINT32 *)calloc(width, sizeof(UINT32));
	T1.precision = c->precision;
	memset(T1.significand, 0, width*sizeof(UINT32));
	
	IEEE_754_FloatNum T2;
	T2.base	= BASE;
	T2.exp		= 0;
	T2.sign	= FALSE;
	T2.significand = (UINT32 *)calloc(width, sizeof(UINT32));
	T2.precision = c->precision;
	memset(T2.significand, 0, width*sizeof(UINT32));
	
	//faster decryption using the CRT
	//T1 = input ^ DP mod P
	// T2 = input ^ DQ mod Q
	paraModExp(c, DP, P, &T1, numProcs, procID);
	paraModExp(c, DQ, Q, &T2, numProcs, procID);

	// T = (T1 - T2) * (Q^-1 mod P) mod P
	paraSub(&T1, &T2, &T, numProcs, procID);
	paraMult(&T, QP, &T1, numProcs, procID);
	paraMod(&T1, P, &T, numProcs, procID);

    //output = T2 + T * Q
	paraMult(&T, Q, &T1, numProcs, procID);
	paraAdd(&T1, &T2, m, numProcs, procID, ADD_CALL);
	
	free(T.significand);
	free(T1.significand);
	free(T2.significand);
}

//c=m^E (mod N)
void RSA_encrypt(IEEE_754_FloatNum *m, 
				 IEEE_754_FloatNum *c,  
				 IEEE_754_FloatNum* N, 
				 IEEE_754_FloatNum* E,  
				 UINT32 procID, 
				 UINT32 numProcs )
{
	paraModExp(m, E, N, c, numProcs, procID);
}

void read_RSA_context(IEEE_754_FloatNum* E,
					  IEEE_754_FloatNum* P,
					  IEEE_754_FloatNum* Q,
					  IEEE_754_FloatNum* N,
					  IEEE_754_FloatNum* D,
					  IEEE_754_FloatNum* DP,
					  IEEE_754_FloatNum* DQ,
					  IEEE_754_FloatNum* QP,
					  UINT32 procID, 
					  UINT32 numProcs)
{
	char *fileName;
	read_IEEE_754_FloatNum(E, numProcs, procID, "KEYS/E.txt");
	E->exp = 4;
	print_IEEE_754_FloatNum(E, ">>E", numProcs, procID); 
	
	fileName = makeFileName("P", KEYLENGTH);
	read_IEEE_754_FloatNum(P, numProcs, procID, fileName);
	print_IEEE_754_FloatNum(P, ">>P", numProcs, procID); 
	
	fileName = makeFileName("Q", KEYLENGTH);
	read_IEEE_754_FloatNum(Q, numProcs, procID, fileName);
	print_IEEE_754_FloatNum(Q, ">>Q", numProcs, procID); 
	
	fileName = makeFileName("N", KEYLENGTH);
	read_IEEE_754_FloatNum(N, numProcs, procID, fileName);
	print_IEEE_754_FloatNum(N, ">>N", numProcs, procID);
	
	fileName = makeFileName("D", KEYLENGTH);
	read_IEEE_754_FloatNum(D, numProcs, procID, fileName);
	print_IEEE_754_FloatNum(D, ">>D", numProcs, procID);
	
	fileName = makeFileName("DP", KEYLENGTH);
	read_IEEE_754_FloatNum(DP, numProcs, procID, fileName);
	print_IEEE_754_FloatNum(DP, ">>DP", numProcs, procID);
	
	fileName = makeFileName("DQ", KEYLENGTH);
	read_IEEE_754_FloatNum(DQ, numProcs, procID, fileName);
	print_IEEE_754_FloatNum(DQ, ">>DQ", numProcs, procID);
	
	fileName = makeFileName("QP", KEYLENGTH);
	read_IEEE_754_FloatNum(QP, numProcs, procID, fileName);
	print_IEEE_754_FloatNum(QP, ">>QP", numProcs, procID);
	free(fileName);
	
}
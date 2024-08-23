#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include "helper.h"
#include "parallel.h"
#include "calc.h"
#include "mpifft.h"
#include "prime.h"
#include "rsa.h"


int main(int argc, char** argv)
{
	double exectime, start;
	int ierr = MPI_Init (&argc, &argv);
	
	if (ierr != MPI_SUCCESS){
		printf("Cannot initialize MPI!\n");
		MPI_Finalize();
		exit(0);
	}
	
	int numProcs;
	MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
	int procID;
	MPI_Comm_rank(MPI_COMM_WORLD, &procID);

	//-----------------------------------------
	//2. Test struct
	UINT32 width = BLOCK_SIZE(procID, numProcs, PRECISION);
	IEEE_754_FloatNum m, c, DP, DQ, QP, D, N, res, E, P, Q, x, y;
	
	x.base	= BASE;
	x.exp	= 0;
	x.sign	= FALSE;
	x.significand = (UINT32 *)calloc(width, sizeof(UINT32));
	x.precision = PRECISION;
	
	y.base	= BASE;
	y.exp	= 0;
	y.sign	= FALSE;
	y.significand = (UINT32 *)calloc(width, sizeof(UINT32));
	y.precision = PRECISION;
	
	P.base	= BASE;
	P.exp	= 0;
	P.sign	= FALSE;
	P.significand = (UINT32 *)calloc(width, sizeof(UINT32));
	P.precision = PRECISION;
	
	Q.base	= BASE;
	Q.exp	= 0;
	Q.sign	= FALSE;
	Q.significand = (UINT32 *)calloc(width, sizeof(UINT32));
	Q.precision = PRECISION;
	
	N.base	= BASE;
	N.exp	= 0;
	N.sign	= FALSE;
	N.significand = (UINT32 *)calloc(width, sizeof(UINT32));
	N.precision = PRECISION;
	
	D.base	= BASE;
	D.exp	= 0;
	D.sign	= FALSE;
	D.significand = (UINT32 *)calloc(width, sizeof(UINT32));
	D.precision = PRECISION;
	
	E.base	= BASE;
	E.exp	= 0;
	E.sign	= FALSE;
	E.significand = (UINT32 *)calloc(width, sizeof(UINT32));
	E.precision = PRECISION;
	
	m.base	= BASE;
	m.exp	= 0;
	m.sign	= FALSE;
	m.significand = (UINT32 *)calloc(width, sizeof(UINT32));
	m.precision = PRECISION;
	
	QP.base	= BASE;
	QP.exp	= 0;
	QP.sign	= FALSE;
	QP.significand = (UINT32 *)calloc(width, sizeof(UINT32));
	QP.precision = PRECISION;
	
	DP.base	= BASE;
	DP.exp	= 0;
	DP.sign	= FALSE;
	DP.significand = (UINT32 *)calloc(width, sizeof(UINT32));
	DP.precision = PRECISION;
	
	DQ.base	= BASE;
	DQ.exp	= 0;
	DQ.sign	= FALSE;
	DQ.significand = (UINT32 *)calloc(width, sizeof(UINT32));
	DQ.precision = PRECISION;
	
	res.base	= BASE;
	res.exp		= 0;
	res.sign	= FALSE;
	res.significand = (UINT32 *)calloc(width, sizeof(UINT32));
	res.precision = PRECISION;
	memset(res.significand, 0, width*sizeof(UINT32));
	
	c.base	= BASE;
	c.exp		= 0;
	c.sign	= FALSE;
	c.significand = (UINT32 *)calloc(width, sizeof(UINT32));
	c.precision = PRECISION;
	memset(c.significand, 0, width*sizeof(UINT32));

	IEEE_754_FloatNum one;
	one.base	= BASE;
	one.exp	= 0;
	one.sign	= FALSE;
	one.significand = (UINT32 *)calloc(width, sizeof(UINT32));
	one.precision = PRECISION;
	memset(one.significand, 0, width*sizeof(UINT32));
	
	if(procID==0)
	{
		one.significand[0] = 1;
	}
	
	if(procID == 0)
    {
		printf("********* SETTINGS *********\n");
		printf("ITERS: %d\n", ITERS);
		printf("PRECISION: %d\n", PRECISION);
		printf("KEYLENGTH: %d\n", KEYLENGTH);
		printf("BASE: %d\n", BASE);
		printf("numProcs: %d\n", numProcs);
		printf("digits of pi: %d\n", (int)log10((double)BASE)*PRECISION);
	}
/*	MPI_Barrier(MPI_COMM_WORLD);
	start = MPI_Wtime();

	MPI_Barrier(MPI_COMM_WORLD);
	if(procID == 0)
    {
		exectime = MPI_Wtime() - start;
		printf( "computed pi time: %f\n", exectime );
		
    }/**/
	char *fileName;
//	fileName = makeFileName("E", PRECISION/2);
//	read_IEEE_754_FloatNum(&E, numProcs, procID, fileName);
//	print_IEEE_754_FloatNum(&E, ">>E", numProcs, procID);

//	x.precision = y.precision = P.precision = res.precision = PRECISION/2;
//	
//	fileName = makeFileName("c", 0);
//	read_IEEE_754_FloatNum(&x, numProcs, procID, fileName);
//	print_IEEE_754_FloatNum(&x, ">>c", numProcs, procID);
//	
//	fileName = makeFileName("DP", KEYLENGTH);
//	read_IEEE_754_FloatNum(&y, numProcs, procID, fileName);
//	print_IEEE_754_FloatNum(&y, ">>DP", numProcs, procID);
//	
//	fileName = makeFileName("P", KEYLENGTH);
//	read_IEEE_754_FloatNum(&P, numProcs, procID, fileName);
//	print_IEEE_754_FloatNum(&P, ">>P", numProcs, procID);
	

	
//	UINT32 maxNum = (UINT32)pow(10,3);
//	unsigned int *globalPrimeArray = (unsigned int *)calloc(maxNum+1, sizeof(unsigned int));
//	
//	makePrimeList(maxNum, numProcs, procID, globalPrimeArray, 2);
	
	
//	RSA_genKey		(procID, numProcs, KEYLENGTH);
//	paraInvMod		(&E, &y, &res, numProcs, procID);
//	paraFloor		(&a, &res, numProcs, procID);
//	paraExtEuclRecursive(&E, &y, &res, &x, numProcs,procID);
//	paraExtendedEuclidean(&x, &y, &res, &m, numProcs,procID);

//	paraExtendedEuclidean2(&x, &y, &res, &m, numProcs,procID);
//	paraGenPrime	(&res, globalPrimeArray, numProcs, procID, KEYLENGTH);
//	paraPow2		(&E, &res,numProcs,procID);
	
//	paraModExpMonty	(&y, &E, &x, &res, numProcs, procID);
//	paraModExp(&x, &y, &P, &res, numProcs, procID);
//	print_IEEE_754_FloatNum(&res, "RES", numProcs, procID);

//	paraRand(&E, numProcs, procID, PRECISION);
//	paraRand(&x, numProcs, procID, PRECISION);
//	paraRand(&y, numProcs, procID, PRECISION);
//	paraDiv(&x, &y, &res, numProcs, procID);
//	printf("Smaller: %d\n", smaller_IEEE(&res, &one, numProcs, procID));
//	printf("Smaller: %d\n", smallerOne_IEEE(&res, numProcs, procID));
//	memCPY(&res, &x, width);
//	paraMod(&x, &y, &res, numProcs, procID);
//	
//	if(procID == 0)
//    {
//		exectime = MPI_Wtime() - start;
//		printf( "computed time TEST1: %f\n", exectime );
//		
//    }
//	print_IEEE_754_FloatNum(&res, "RES1", numProcs, procID);
//	
//	MPI_Barrier(MPI_COMM_WORLD);
//	start = MPI_Wtime();
//	
//	memset(res.significand, 0, width*sizeof(UINT32));
//	paraMod3(&x, &y, &res, numProcs, procID);
//	
//	if(procID == 0)
//    {
//		exectime = MPI_Wtime() - start;
//		printf( "computed time TEST3: %f\n", exectime );
//		
//    }
//	print_IEEE_754_FloatNum(&res, "RES3", numProcs, procID);

//	paraGCD2(&x, &y, &res, numProcs, procID);
//	paraNorm(&y, numProcs, procID);
//	makeSameExp(&y, &x, numProcs, procID);
//	paraAdd(&x, &y, &res, numProcs, procID, ADD_CALL);
//	paraMult(&x, &y, &res, numProcs, procID);
//	paraSub(&P, &m, &res, numProcs, procID);

//	paraDiv(&x, &E, &res, numProcs, procID);
//	print_IEEE_754_FloatNum(&res, "nach DIV", numProcs, procID);
//	paraFloor(&res, &b, numProcs, procID);
//	print_IEEE_754_FloatNum(&b, "nach FLOOR", numProcs, procID);
//	paraMult(&E, &E, &res, numProcs, procID);
//	print_IEEE_754_FloatNum(&a, "MULT", numProcs, procID);
//	paraSub(&y, &a, &res, numProcs, procID);
	
//	paraSqrt(&y, &res, numProcs, procID, 1);
//	print_IEEE_754_FloatNum(&res, "res SQRT", numProcs, procID);
	
	read_RSA_context(&E,&P,&Q,&N,&D,&DP,&DQ,&QP,procID,numProcs);
	read_IEEE_754_FloatNum(&m, numProcs, procID, "KEYS/m.txt");
	print_IEEE_754_FloatNum(&m, ">>m", numProcs, procID);
	
	MPI_Barrier(MPI_COMM_WORLD);
	start = MPI_Wtime();
	
	RSA_encrypt(&m, &c, &N, &E, procID, numProcs );
	
	if(procID == 0)
    {
		exectime = MPI_Wtime() - start;
		printf( "computed time ENC: %f\n", exectime );
		
    }
	
	fileName = makeFileName("c", 0);
	write_IEEE_754_FloatNum(&c, numProcs, procID, fileName);
	
	print_IEEE_754_FloatNum(&c, "c decryption", numProcs, procID);
  
	memset(m.significand, 0, width*sizeof(UINT32));
	m.exp = 0;
	print_IEEE_754_FloatNum(&m, "m vorher", numProcs, procID);
	
	MPI_Barrier(MPI_COMM_WORLD);
	start = MPI_Wtime();
	paraModExp(&c, &D, &N, &m, numProcs, procID);
	if(procID == 0)
	{
        exectime = MPI_Wtime() - start;
        printf( "computed time decryption: %f\n", exectime );
	}
	print_IEEE_754_FloatNum(&m, "m decrypted", numProcs, procID);
	

//	memCPY(&Q, &one, width);
//	memCPY(&P, &one, width);
//	memCPY(&DP, &one, width);
//	memCPY(&DQ, &one, width);
//	memCPY(&QP, &one, width);
//	memCPY(&N, &one, width);
//	memCPY(&D, &one, width);
//	memCPY(&c, &one, width);
//	
//	m.precision = c.precision = PRECISION/2;
//	
//	DP.precision = DQ.precision = QP.precision = Q.precision = P.precision = PRECISION/2;
//	read_RSA_context(&E,&P,&Q,&N,&D,&DP,&DQ,&QP,procID,numProcs);
//	
//	fileName = makeFileName("c", 0);
//	read_IEEE_754_FloatNum(&c, numProcs, procID, fileName);
//	print_IEEE_754_FloatNum(&c, "c CRT", numProcs, procID);

	
	MPI_Barrier(MPI_COMM_WORLD);
	start = MPI_Wtime();
	RSA_decrypt(&c, &m, &Q, &P, &QP, &DP, &DQ, procID, numProcs );
	if(procID == 0)
    {
		exectime = MPI_Wtime() - start;
		printf( "computed time CRT: %f\n", exectime );
		
    }
	print_IEEE_754_FloatNum(&m, "m", numProcs, procID);


	
	
//	print_IEEE_754_FloatNum(&res, "RES", numProcs, procID);
//	print_IEEE_754_FloatNum(&m, "m", numProcs, procID);

//	print_IEEE_754_FloatNum(&E, "E", numProcs, procID);
//	print_IEEE_754_FloatNum(&x, "X", numProcs, procID);
//	print_IEEE_754_FloatNum(&y, "Y", numProcs, procID);
//
//	fileName = makeFileName("E", PRECISION);
//	write_IEEE_754_FloatNum(&E, numProcs, procID, fileName);
//	
//	fileName = makeFileName("Y", PRECISION);
//	write_IEEE_754_FloatNum(&y, numProcs, procID, fileName);
//	
//	fileName = makeFileName("X", PRECISION);
//	write_IEEE_754_FloatNum(&x, numProcs, procID, fileName);

//free memory
	free(DP.significand);
	free(DQ.significand);
	free(QP.significand);
	free(D.significand);
	free(N.significand);
	free(E.significand);
	free(P.significand);
	free(Q.significand);
	free(res.significand);
	free(x.significand);
	free(y.significand);


	
	
	
	MPI_Finalize();
		
	return 0;
	
}

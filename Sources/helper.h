#ifndef __HELPER_H__
#define __HELPER_H__

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

typedef signed char         INT8;          /*        -128 .. +127            */
typedef unsigned char       UINT8;         /*           0 .. 255             */
typedef signed short        INT16;         /*      -32768 .. +32767          */
typedef unsigned short      UINT16;        /*           0 .. 65535           */ 
typedef signed int         INT32;         /* -2147483648 .. +2147483647     */
typedef unsigned int       UINT32;        /*           0 .. 4294967295      */
typedef UINT8               BOOLEAN;       /* for use with TRUE/FALSE        */


#ifndef TRUE                               /* conditional check */
#define TRUE      ((BOOLEAN) 1)
#endif

#ifndef FALSE                              /* conditional check */
#define FALSE     ((BOOLEAN) 0)
#endif

#define ZERO 0
#ifndef NULL_PTR
#define NULL_PTR  (0)
#endif

//#define INTEGER TRUE
#define KEYLENGTH 2048 //bin length of m=pq
#define ITERS 20 //define iterations for div, sqrt, pi 
#define CARRY_TAG 2
#define BASE 10
#define PRECISION KEYLENGTH //define precision = 4 x KEYLENGTH
//#define PRECISION 2048//2^11
//#define PRECISION 4096//2^12
//#define PRECISION 8192 //2^13
//#define PRECISION 16384	//2^14
//#define PRECISION 32768	//2^15
//#define PRECISION 65536	//2^16
//#define PRECISION 262144	//2^18
//#define PRECISION 1048576 //2^20
//#define PRECISION 16777216 //2^24
#define SUB_CALL TRUE
#define ADD_CALL FALSE
#define LEFT FALSE
#define RIGHT TRUE

#define LOW_INDEX(s, p, n)	((s)*(n)/(p))
#define HIGH_INDEX(s, p, n) (LOW_INDEX((s)+1, p, n)-1)
#define BLOCK_SIZE(s, p, n)	(HIGH_INDEX(s, p, n) - LOW_INDEX(s, p, n)+1)
#define BLOCK_OWNER(j, p, n) ((p)*((j)+1)-1)/(n))

#define SZDBL (sizeof(double))
#define SZINT (sizeof(int))
#define MAX(a,b) ((a)>(b) ? (a) : (b))
#define MIN(a,b) ((a)<(b) ? (a) : (b))

typedef enum
	{
		ERROR_INEXACT, //set if the rounded (and returned) value is different from the mathematically exact result of the operation.
		ERROR_UNDERFLOW, //set if the rounded value is tiny (as specified in IEEE 754) and inexact (or maybe limited to if it has denormalisation loss, as per the 1984 version of IEEE 754), returning a subnormal value including the zeros.
		ERROR_OVERFLOW, //set if the absolute value of the rounded value is too large to be represented. An infinity or maximal finite value is returned, depending on which rounding is used.
		ERROR_DEVIDE_BY_ZERO, //set if the result is infinite given finite operands, returning an infinity, either +inf or −inf.
		ERROR_INVALID, //set if a real-valued result cannot be returned e.g. sqrt(−1) or 0/0, returning a quiet NaN.
		NO_ERROR	
	}IEEE_754_ErrorType;

//X = s*B^e*sum(x_i*B^i,i = 0..n-1)
typedef struct
	{
		UINT32 base; //B
		INT32 exp; //e
		BOOLEAN sign; //s
		UINT32 *significand; //x[n], x = [0,9]
		UINT32 precision; //n, how many digits total for memory allocation
	}IEEE_754_FloatNum;

double *vecallocd(int n);
int *vecalloci(int n);
double **matallocd(int m, int n);
void vecfreed(double *pd);
void vecfreei(int *pi);
void matfreed(double **ppd);
void increasePrecision(IEEE_754_FloatNum *oldPrecNum,
					   IEEE_754_FloatNum *newPrecNum,
					   UINT32 newPRecision,
					   UINT32 numProcs, 
					   UINT32 procID);
void reducePrecision(IEEE_754_FloatNum *oldPrecNum,
					 IEEE_754_FloatNum *newPrecNum,
					 UINT32 newPRecision,
					 UINT32 numProcs, 
					 UINT32 procID);

void toIEEE_754(UINT32 intNumber,IEEE_754_FloatNum *x,UINT32 numProcs,UINT32 procID);
char* makeFileName(char *input, int keyLength);
void comparePI(IEEE_754_FloatNum *number, IEEE_754_FloatNum *number2, UINT32 numProcs, UINT32 procID);
void read_PI_1000000(IEEE_754_FloatNum *number, UINT32 numProcs, UINT32 procID);
void read_IEEE_754_FloatNum(IEEE_754_FloatNum *number, UINT32 numProcs, UINT32 procID, char *fileName);
void write_IEEE_754_FloatNum(IEEE_754_FloatNum *number, UINT32 numProcs, UINT32 procID, char *fileName);
void print_IEEE_754_FloatNum(IEEE_754_FloatNum *number, char *name, UINT32 numProcs, UINT32 procID);
int min( int a, int b);
int max( int a, int b);
BOOLEAN isODD(IEEE_754_FloatNum *x,UINT32 numProcs,UINT32 procID);
BOOLEAN equal_IEEE(IEEE_754_FloatNum *a, IEEE_754_FloatNum *b,UINT32 numProcs, UINT32 procID);
BOOLEAN smallerOne_IEEE(IEEE_754_FloatNum *a, UINT32 numProcs, UINT32 procID);
BOOLEAN smaller_IEEE(IEEE_754_FloatNum *a, IEEE_754_FloatNum *b,UINT32 numProcs, UINT32 procID);
BOOLEAN minMax_IEEE(IEEE_754_FloatNum *min, IEEE_754_FloatNum *max,IEEE_754_FloatNum *a, IEEE_754_FloatNum *b, BOOLEAN sign,  UINT32 numProcs, UINT32 procID);
void toCyclic(UINT32 *x, double *compX,UINT32 lengthOfx,UINT32 numProcs,UINT32 procID);
void toBlock(double *a,double *b,UINT32 p, UINT32 s,UINT32 n);
IEEE_754_FloatNum create_IEEE_NULL(UINT32 precision, UINT32 numProcs, UINT32 procID);
void kill_IEEE_NULL(IEEE_754_FloatNum x);
void memCPY(IEEE_754_FloatNum *a, IEEE_754_FloatNum *b,UINT32 width);
#endif     

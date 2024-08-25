#include "calc.h"
#include "prime.h"
#include <mpi.h>
#include "mpifft.h"

//x*2 = res = x+x;
void paraMultBy2(IEEE_754_FloatNum *x,
				  IEEE_754_FloatNum *res,
				  UINT32 numProcs, 
				  UINT32 procID)
{
	paraAdd(x, x, res, numProcs,procID, ADD_CALL);
}

//integer div by two with remainder as return value
UINT32 paraDivBy2(IEEE_754_FloatNum *x,
			  IEEE_754_FloatNum *res,
			  UINT32 numProcs, 
			  UINT32 procID)
{
	UINT32 width = BLOCK_SIZE(procID, numProcs, x->precision);
	UINT32 rest = 0;
	UINT32 digs = x->exp;

	UINT32 i, j, k;
	for(i=0, j=procID*width; i<width && j<digs; i++, j++)
	{
		;
	}
	//if x uneven set rest to 1 and subtract 1
	if(i < width && j<=digs)
	{
		if(x->significand[i]%2 != 0)
		{
			rest = 1;
			x->significand[i] -= 1;
		}
	}
	
	MPI_Allreduce(&rest,&rest,1,MPI_UNSIGNED,MPI_SUM,MPI_COMM_WORLD);

	
	IEEE_754_FloatNum oneHalf;
	oneHalf.base	= BASE;
	oneHalf.exp		= -1;
	oneHalf.sign	= FALSE;
	oneHalf.significand = (UINT32 *)calloc(width, sizeof(UINT32));
	oneHalf.precision = x->precision;
	memset(oneHalf.significand, 0, width*sizeof(UINT32));
	
	if(procID == 0)
		oneHalf.significand[0] = 5;
	
	paraMult(x, &oneHalf, res, numProcs, procID);
	
	free(oneHalf.significand);
	return rest;
}

//sqaurring of a number with modified multiplication
void paraPow2(IEEE_754_FloatNum *x,
			  IEEE_754_FloatNum *res,
			  UINT32 numProcs, 
			  UINT32 procID)
{
	UINT32 width = BLOCK_SIZE(procID, numProcs, res->precision);
	
	IEEE_754_FloatNum null;
	null.base	= BASE;
	null.exp	= 0;
	null.sign	= FALSE;
	null.significand = (UINT32 *)calloc(width, sizeof(UINT32));
	null.precision = x->precision;
	memset(null.significand, 0, width*sizeof(UINT32));
		
	if(equal_IEEE(x, &null, numProcs, procID) == TRUE)
	{
		memCPY(res, &null, width);
		free(null.significand);
		return;
	}
	
	free(null.significand);
	
	IEEE_754_FloatNum one;
	one.base	= BASE;
	one.exp	= 0;
	one.sign	= FALSE;
	one.significand = (UINT32 *)calloc(width, sizeof(UINT32));
	one.precision = x->precision;
	memset(one.significand, 0, width*sizeof(UINT32));	
	
	if(procID == 0)
	{
		one.significand[0]=1;
	}
	
	if(equal_IEEE(x, &one, numProcs, procID) == TRUE)
	{
		memCPY(res, x, width);
		free(one.significand);
		return;
	}
	
	free(one.significand);
	
	UINT32 i, finExp;
	int j, jglob, p, s;
	p = numProcs;
	s = procID;
	double *compX = (double*)calloc(width*4, sizeof(double));
	double *compY = (double*)calloc(width*4, sizeof(double));
	
	finExp = x->exp + x->exp;
	
	//block to cyclic redistribution
	toCyclic(x->significand, compX, width, numProcs, procID);
	
	for(i=1;i<=width;i++)
	{
		compX[2*(width-i)]=compX[width-i];
		compX[2 * (width-i) + 1] = 0; 
		compX[2 * (i-1 + width)] = 0; 
		compX[2 * (i-1 + width) + 1] = 0;
	}
	
	//reset res
	res->exp = 0;
	
	//perform FFT on numbers
	MPI_FFT(compX, x->precision*2, procID, numProcs, 1);

	double reX, imX, reY, imY;
	for(i = 0; i < width*2; i++)
	{
		reX = compX[2*i+0];
		imX = compX[2*i+1];
		reY = compX[2*i+0];
		imY = compX[2*i+1];
		
		compX[2*i+0] = reX*reY-imX*imY;
		compX[2*i+1] = reX*imY+reY*imX;
	}
	
	//perform an inverse FFT
	
	//remove imaginary part
	for(j=0, i=0; j < width*4; j=j+2, i++)
	{
        //jglob= j*numProcs+procID;
		compX[i] = compX[j];
	}
	
	
	//cyclic to block redistribution compX->compY, but just for the half of the array
	toBlock(compX, compY, numProcs, procID, 2*width);

	IEEE_754_FloatNum temp;
	temp.base	= BASE;
	temp.exp	= 0;
	temp.sign	= FALSE;
	temp.significand = (UINT32 *)calloc(2*width, sizeof(UINT32));
	temp.precision = 2*x->precision;
	memset(temp.significand, 0, width*sizeof(UINT32));
	
	IEEE_754_FloatNum IEEE_Null = create_IEEE_NULL(2*x->precision, numProcs, procID);
	IEEE_754_FloatNum IEEE_Temp = create_IEEE_NULL(2*x->precision, numProcs, procID);

	for(i=0; i < 2*width; i++){
		IEEE_Temp.significand[i]=(int)(compY[i]+0.5);
	}
		
	IEEE_Temp.exp = 0;
	//carry add on the result
	paraAdd(&IEEE_Temp, &IEEE_Null, &temp, numProcs, procID, ADD_CALL);

	reducePrecision(&temp, res, x->precision, numProcs, procID);
	
	res->exp = finExp + res->exp;
	res->sign= x->sign | x->sign;
	
	//normalize res
	paraNorm(res, numProcs, procID);	
	
	free(compX);
	free(compY);
	
	kill_IEEE_NULL(IEEE_Temp);
	kill_IEEE_NULL(IEEE_Null);
	kill_IEEE_NULL(temp);
}

// res = base^exp (mod modulus) with binary right-to-left method
void paraModExp(IEEE_754_FloatNum *base, 
				IEEE_754_FloatNum *exp,
				IEEE_754_FloatNum *modulus,
				IEEE_754_FloatNum *res,
				UINT32 numProcs, 
				UINT32 procID)
{
	UINT32 width = BLOCK_SIZE(procID, numProcs, base->precision);
	
	IEEE_754_FloatNum tempE;
	tempE.significand = (UINT32 *)calloc(width, sizeof(UINT32));
	memCPY(&tempE, exp, width);
	
	IEEE_754_FloatNum tempBase;
	tempBase.significand = (UINT32 *)calloc(width, sizeof(UINT32));
	memCPY(&tempBase, base, width);
	
	IEEE_754_FloatNum tempRes;
	tempRes.significand = (UINT32 *)calloc(width, sizeof(UINT32));
	
	IEEE_754_FloatNum one;
	one.base	= BASE;
	one.exp		= 0;
	one.sign	= FALSE;
	one.significand = (UINT32 *)calloc(width, sizeof(UINT32));
	one.precision = base->precision;
	memset(one.significand, 0, width*sizeof(UINT32));
	memset(res->significand, 0, width*sizeof(UINT32));
	
	IEEE_754_FloatNum temp;
	temp.base	= BASE;
	temp.exp		= 0;
	temp.sign	= FALSE;
	temp.significand = (UINT32 *)calloc(width, sizeof(UINT32));
	temp.precision = base->precision;
	memset(temp.significand, 0, width*sizeof(UINT32));
	
	IEEE_754_FloatNum null;
	null.base	= BASE;
	null.exp	= 0;
	null.sign	= FALSE;
	null.significand = (UINT32 *)calloc(width, sizeof(UINT32));
	null.precision = base->precision;
	memset(null.significand, 0, width*sizeof(UINT32));
	memset(res->significand, 0, width*sizeof(UINT32));
	
	UINT32 rest = 0;
	
	//res = 1;
	if(procID == 0)
	{
		res->significand[0] = 1;
		one.significand[0]=1;
	}
	memCPY(&tempRes, res, width);

	int i = 0;
	do
	{
		
		i++;
		
		//Multiplication with 0.5
		rest = paraDivBy2(&tempE, &tempE, numProcs, procID);
		
		if(rest == 1)
		{
			//result = (result * base) mod modulus
			paraMult(res, &tempBase, &tempRes, numProcs, procID);
			paraMod(&tempRes, modulus, res, numProcs, procID);
		}
				
		//base = (base * base) mod modulus
		paraPow2(&tempBase, &tempRes, numProcs,procID);
		paraMod(&tempRes, modulus, &tempBase, numProcs, procID);
				
	}while(equal_IEEE(&tempE, &null, numProcs,procID) != TRUE);
	
	free(temp.significand);
	free(tempE.significand);
	free(tempRes.significand);
	free(tempBase.significand);
	free(one.significand);
	free(null.significand);
	
	
}

//res = x^-1 mod modulus
void paraInvMod(IEEE_754_FloatNum *x, 
				IEEE_754_FloatNum *modulus,
				IEEE_754_FloatNum *res,
				UINT32 numProcs, 
				UINT32 procID)
{
	UINT32 width = BLOCK_SIZE(procID, numProcs, x->precision);

	IEEE_754_FloatNum aux;
	aux.base	= BASE;
	aux.exp	= 0;
	aux.sign	= FALSE;
	aux.significand = (UINT32 *)calloc(width, sizeof(UINT32));
	aux.precision = x->precision;
	memset(aux.significand, 0, width*sizeof(UINT32));
	
	paraExtendedEuclidean(x, modulus, res, &aux, numProcs,procID);

	free(aux.significand);
}

void paraRand(IEEE_754_FloatNum *x, 
			 UINT32 numProcs, 
			 UINT32 procID,
			 UINT32 length)
{
	UINT32 width = BLOCK_SIZE(procID, numProcs, x->precision);
	
	IEEE_754_FloatNum tempRes;
	tempRes.base	= BASE;
	tempRes.exp	= 0;
	tempRes.sign	= FALSE;
	tempRes.significand = (UINT32 *)calloc(width, sizeof(UINT32));
	tempRes.precision = x->precision;
	
	
	int i;
	
	for(i=0; i<width; i++)
	{
		x->significand[i] = rand() % BASE;		
	}
	
	x->exp = length-1;
	x->sign = FALSE;

	//avoid leading zero
	if(procID == 0 && x->significand[0]==0)
		x->significand[0] = 3;
	
	paraFloor(x, &tempRes, numProcs, procID);
	memCPY(x, &tempRes, width);
		
}

// res = x (mod y) w/ paraDiv
void paraMod1(IEEE_754_FloatNum *x, 
			 IEEE_754_FloatNum *y,
			 IEEE_754_FloatNum *res,
			 UINT32 numProcs, 
			 UINT32 procID)
{
	UINT32 width = BLOCK_SIZE(procID, numProcs, x->precision);
	
	IEEE_754_FloatNum null;
	null.base	= BASE;
	null.exp	= 0;
	null.sign	= FALSE;
	null.significand = (UINT32 *)calloc(width, sizeof(UINT32));
	null.precision = x->precision;
	memset(null.significand, 0, width*sizeof(UINT32));
	
	IEEE_754_FloatNum one;
	one.base	= BASE;
	one.exp	= 0;
	one.sign	= FALSE;
	one.significand = (UINT32 *)calloc(width, sizeof(UINT32));
	one.precision = x->precision;
	memset(one.significand, 0, width*sizeof(UINT32));
	
	if(procID==0)
	{
		one.significand[0] = 1;
	}

	
	IEEE_754_FloatNum temp;
	temp.base	= BASE;
	temp.exp	= 0;
	temp.sign	= FALSE;
	temp.significand = (UINT32 *)calloc(width, sizeof(UINT32));
	temp.precision = x->precision;
	memset(temp.significand, 0, width*sizeof(UINT32));


	if(equal_IEEE(x, y, numProcs, procID) == TRUE)
	{
		memCPY(res, &null, width);
		goto exit;
	}
	
	if(x->sign == TRUE && y->sign == FALSE)
	{
		paraDiv(x, y, res, numProcs, procID);
		
		res->sign = FALSE;
	
		//res < |1|
		if(smallerOne_IEEE(res, numProcs, procID) == TRUE )
		{
			paraAdd(y, x, res, numProcs, procID, ADD_CALL);
			goto exit;
		}
		

		paraFloor(res, &temp, numProcs, procID);
		
		paraAdd(&temp, &one, res, numProcs, procID, ADD_CALL);
		res->sign = TRUE;

	}
	
	if((x->sign == TRUE && y->sign == TRUE) || (x->sign == FALSE && y->sign == FALSE))
	{
		paraDiv(x, y, res, numProcs, procID);

		//res < 1
		if(smallerOne_IEEE(res, numProcs, procID) == TRUE )
		{
			memCPY(res, x, width);
			goto exit;
		}
		
		paraFloor(res, &temp, numProcs, procID);
		memCPY(res, &temp, width);
	}
	

	paraMult(res, y, &temp, numProcs, procID);
	paraSub(x, &temp, res, numProcs, procID);


exit:	
	free(null.significand);
	free(one.significand);
	free(temp.significand);

	return;
}


// res = x (mod y) w/o paraDiv
void paraMod(IEEE_754_FloatNum *x, 
			 IEEE_754_FloatNum *y,
			 IEEE_754_FloatNum *res,
			 UINT32 numProcs, 
			 UINT32 procID)
{
	UINT32 width = BLOCK_SIZE(procID, numProcs, x->precision);
	
	IEEE_754_FloatNum null;
	null.base	= BASE;
	null.exp	= 0;
	null.sign	= FALSE;
	null.significand = (UINT32 *)calloc(width, sizeof(UINT32));
	null.precision = x->precision;
	memset(null.significand, 0, width*sizeof(UINT32));
	
	IEEE_754_FloatNum one;
	one.base	= BASE;
	one.exp	= 0;
	one.sign	= FALSE;
	one.significand = (UINT32 *)calloc(width, sizeof(UINT32));
	one.precision = x->precision;
	memset(one.significand, 0, width*sizeof(UINT32));
	
	if(procID==0)
	{
		one.significand[0] = 1;
	}
	
	
	IEEE_754_FloatNum temp;
	temp.significand = (UINT32 *)calloc(width, sizeof(UINT32));
	memCPY(&temp, x, width);
	
	if(equal_IEEE(x, y, numProcs, procID) == TRUE)
	{
		memCPY(res, &null, width);
		goto exit;
	}
	int i=0;
	
	if(x->sign == TRUE && y->sign == FALSE)
	{
		paraDiv(x, y, res, numProcs, procID);
		
		res->sign = FALSE;
		
		//res < |1|
		if(smallerOne_IEEE(res, numProcs, procID) == TRUE )
		{
			paraAdd(y, x, res, numProcs, procID, ADD_CALL);
			goto exit;
		}
		
		
		paraFloor(res, &temp, numProcs, procID);
		
		paraAdd(&temp, &one, res, numProcs, procID, ADD_CALL);
		res->sign = TRUE;
		
		paraMult(res, y, &temp, numProcs, procID);
		paraSub(x, &temp, res, numProcs, procID);		
		
	}
	
	INT32 yExp = y->exp;
	if((x->sign == TRUE && y->sign == TRUE) || (x->sign == FALSE && y->sign == FALSE))
	{
		//if x < y then res = x
		if(smaller_IEEE(x, y, numProcs, procID) == TRUE )
		{
			memCPY(res, x, width);
			goto exit;
		}

		while(smaller_IEEE(&temp, y, numProcs, procID) != TRUE)// && i<10)
		{
			i++;
			y->exp = temp.exp;
			if(smaller_IEEE(&temp, y, numProcs, procID) != TRUE)
			{
				paraSub(&temp, y, res, numProcs, procID);
				memCPY(&temp, res, width);
			}
			else
			{
				y->exp = temp.exp-1;
				paraSub(&temp, y, res, numProcs, procID);
				memCPY(&temp, res, width);				
			}
			memCPY(&temp, res, width);
			y->exp = yExp;
		}
		
		y->exp = yExp;
	}
		
exit:	
	free(null.significand);
	free(one.significand);
	free(temp.significand);
	
	return;
}

void paraFloor(IEEE_754_FloatNum *x,
			   IEEE_754_FloatNum *res,
			   UINT32 numProcs, 
			   UINT32 procID)
{
	UINT32 width = BLOCK_SIZE(procID, numProcs, x->precision);
	UINT32 digs = x->exp+1;
	memset(res->significand, 0, width*sizeof(UINT32));
	UINT32 i, j;

	if(x->exp < 0)
	{
		res->exp = 0;
		res->sign = FALSE;
		return;
	}
	
	for(i=0, j=procID*width; i<width && j<digs; i++, j++)
	{
		res->significand[i] = x->significand[i];
	}
	
	res->exp = x->exp;
	res->sign = x->sign;
}

//extended Euclidian algorithm
void paraExtendedEuclidean(IEEE_754_FloatNum *a,
						   IEEE_754_FloatNum *b,
						   IEEE_754_FloatNum *x, 
						   IEEE_754_FloatNum *y,
						   UINT32 numProcs, 
						   UINT32 procID)
{
	UINT32 width = BLOCK_SIZE(procID, numProcs, b->precision);
	
	IEEE_754_FloatNum aTemp;
	IEEE_754_FloatNum bTemp;
	
	aTemp.significand = (UINT32 *)calloc(width, sizeof(UINT32));
	bTemp.significand = (UINT32 *)calloc(width, sizeof(UINT32));
	
	UINT32 newWidth = BLOCK_SIZE(procID, numProcs, 2*b->precision);
	
	IEEE_754_FloatNum aTempNew;
	aTempNew.base	= BASE;
	aTempNew.exp	= 0;
	aTempNew.sign	= FALSE;
	aTempNew.significand = (UINT32 *)calloc(newWidth, sizeof(UINT32));
	aTempNew.precision = b->precision;
	memset(aTempNew.significand, 0, width*sizeof(UINT32));
		
	IEEE_754_FloatNum bTempNew;
	bTempNew.base	= BASE;
	bTempNew.exp	= 0;
	bTempNew.sign	= FALSE;
	bTempNew.significand = (UINT32 *)calloc(newWidth, sizeof(UINT32));
	bTempNew.precision = b->precision;
	memset(bTempNew.significand, 0, width*sizeof(UINT32));
	
	IEEE_754_FloatNum newRes;
	newRes.base	= BASE;
	newRes.exp	= 0;
	newRes.sign	= FALSE;
	newRes.significand = (UINT32 *)calloc(newWidth, sizeof(UINT32));
	newRes.precision = b->precision;
	memset(newRes.significand, 0, width*sizeof(UINT32));
	
	memCPY(&aTemp, a, width);
	memCPY(&bTemp, b, width);
	
	increasePrecision(&aTemp, &aTempNew, 2*b->precision,numProcs,procID);
	increasePrecision(&bTemp, &bTempNew, 2*b->precision,numProcs,procID);
	
	IEEE_754_FloatNum null;
	null.base	= BASE;
	null.exp	= 0;
	null.sign	= FALSE;
	null.significand = (UINT32 *)calloc(width, sizeof(UINT32));
	null.precision = b->precision;
	memset(null.significand, 0, width*sizeof(UINT32));
	
	IEEE_754_FloatNum aux;
	aux.base	= BASE;
	aux.exp	= 0;
	aux.sign	= FALSE;
	aux.significand = (UINT32 *)calloc(width, sizeof(UINT32));
	aux.precision = b->precision;
	memset(aux.significand, 0, width*sizeof(UINT32));
	
	IEEE_754_FloatNum one;
	one.base	= BASE;
	one.exp	= 0;
	one.sign	= FALSE;
	one.significand = (UINT32 *)calloc(width, sizeof(UINT32));
	one.precision = b->precision;
	memset(one.significand, 0, width*sizeof(UINT32));
	
	IEEE_754_FloatNum lastX;
	lastX.base	= BASE;
	lastX.exp	= 0;
	lastX.sign	= FALSE;
	lastX.significand = (UINT32 *)calloc(width, sizeof(UINT32));
	lastX.precision = b->precision;
	memset(lastX.significand, 0, width*sizeof(UINT32));
	
	IEEE_754_FloatNum lastY;
	lastY.base	= BASE;
	lastY.exp	= 0;
	lastY.sign	= FALSE;
	lastY.significand = (UINT32 *)calloc(width, sizeof(UINT32));
	lastY.precision = b->precision;
	memset(lastY.significand, 0, width*sizeof(UINT32));
	
	IEEE_754_FloatNum quotient;
	quotient.base	= BASE;
	quotient.exp	= 0;
	quotient.sign	= FALSE;
	quotient.significand = (UINT32 *)calloc(width, sizeof(UINT32));
	quotient.precision = b->precision;
	memset(quotient.significand, 0, width*sizeof(UINT32));
	
	IEEE_754_FloatNum temp;
	temp.base	= BASE;
	temp.exp	= 0;
	temp.sign	= FALSE;
	temp.significand = (UINT32 *)calloc(width, sizeof(UINT32));
	temp.precision = b->precision;
	memset(temp.significand, 0, width*sizeof(UINT32));
	
	if(procID == 0)
	{
		one.significand[0] = 1;
	}
	
	//initialize
	memCPY(x, &null, width);
	memCPY(y, &one, width);
	memCPY(&lastX, &one, width);
	memCPY(&lastY, &null, width);

	
	int i=1;
	
	//while b != 0
	while(equal_IEEE(&bTemp, &null, numProcs, procID) != TRUE)
	{

		//temp = b;	 
		memCPY(&temp, &bTemp, width);

		//quotient = a / b;	
		increasePrecision(&aTemp, &aTempNew, 2*b->precision,numProcs,procID);
		increasePrecision(&bTemp, &bTempNew, 2*b->precision,numProcs,procID);
		paraDiv(&aTempNew, &bTempNew, &newRes, numProcs, procID);

		reducePrecision(&newRes, &aux, b->precision, numProcs, procID);
		paraFloor(&aux, &quotient, numProcs, procID);

		//b = a % b;
		paraMod(&aTemp, &bTemp, &aux, numProcs, procID);
		memCPY(&bTemp, &aux, width);

		//a = temp;
		memCPY(&aTemp, &temp, width);
		
		//temp = x;
		memCPY(&temp, x, width);		
		
		//x = lastx-quotient*x;
		paraMult(&quotient, x, &aux, numProcs, procID);
		paraSub(&lastX, &aux, x, numProcs, procID);
		
		//lastx = temp;
		memCPY(&lastX, &temp, width);
		
		//temp = y;
		memCPY(&temp, y, width);
		
		//y = lasty-quotient*y;
		paraMult(&quotient, y, &aux, numProcs, procID);
		paraSub(&lastY, &aux, y, numProcs, procID);
		
		//lasty = temp;
		memCPY(&lastY, &temp, width);
		
		i++;
			
	}
	
	memCPY(x, &lastX, width);
	memCPY(y, &lastY, width);
	
	free(null.significand);
	free(one.significand);
	free(lastX.significand);
	free(lastY.significand);
	free(quotient.significand);
	free(temp.significand);
	free(aux.significand);
	free(aTempNew.significand);
	free(bTempNew.significand);
	free(newRes.significand);


}

//Euclidian gcd algorithm
void paraGCD(IEEE_754_FloatNum *x, 
			 IEEE_754_FloatNum *y,
			 IEEE_754_FloatNum *res,
			 UINT32 numProcs, 
			 UINT32 procID)
{
	UINT32 width = BLOCK_SIZE(procID, numProcs, x->precision);
	IEEE_754_FloatNum max;
	IEEE_754_FloatNum min;
	
	IEEE_754_FloatNum temp;
	temp.base	= BASE;
	temp.exp	= 0;
	temp.sign	= FALSE;
	temp.significand = (UINT32 *)calloc(width, sizeof(UINT32));
	temp.precision = x->precision;
	memset(temp.significand, 0, width*sizeof(UINT32));
	
	IEEE_754_FloatNum temp1;
	temp1.base	= BASE;
	temp1.exp	= 0;
	temp1.sign	= FALSE;
	temp1.significand = (UINT32 *)calloc(width, sizeof(UINT32));
	temp1.precision = x->precision;
	memset(temp1.significand, 0, width*sizeof(UINT32));
	
	IEEE_754_FloatNum null;
	null.base	= BASE;
	null.exp	= 0;
	null.sign	= FALSE;
	null.significand = (UINT32 *)calloc(width, sizeof(UINT32));
	null.precision = x->precision;
	memset(null.significand, 0, width*sizeof(UINT32));
	
	min.significand = (UINT32 *)calloc(width, sizeof(UINT32));
	max.significand = (UINT32 *)calloc(width, sizeof(UINT32));
	
	
	int i=0;

	//get the smaller number
	minMax_IEEE(&min, &max, x, y, FALSE, numProcs, procID);

	while(equal_IEEE(&min, &null, numProcs, procID) != TRUE)
	{
		
		//t = b
		memCPY(&temp, &min, width);

		//b = a mod b
		paraMod(&max, &min, &temp1, numProcs, procID);
		memCPY(&min, &temp1, width);
		
		//a = t
		memCPY(&max, &temp, width);
		
		i++;
				
	}
	memCPY(res, &max, width);
	
	free(temp.significand);
	free(temp1.significand);
	free(min.significand);
	free(max.significand);
	free(null.significand);

}

UINT32 addWithCarry(UINT32 *x, 
					 UINT32 *y, 
					 UINT32 *res,
					 UINT32 base,
					 UINT32 width,
					 UINT32 offset){
	UINT32 carry = FALSE;
	UINT32 temp = 0; 
	UINT32 i = width+offset;
	UINT32 j;
	
	for(j=0; j < width; j++)
	{
		i--;
		
		temp = x[i]+y[i]+carry;
		carry = temp/base;
		res[i] = temp%base;
	}
	
	return carry;
}

void paraAdd(IEEE_754_FloatNum *x, 
			   IEEE_754_FloatNum *y,
			   IEEE_754_FloatNum *res,
			   UINT32 numProcs, 
			   UINT32 procID,
			   BOOLEAN callID)
{
	UINT32 width = BLOCK_SIZE(procID, numProcs, x->precision);
	
	int j, k=0;
	UINT32 buf, maxCarry=1;
	
	UINT32 *carry;
	UINT32 *localCarry;
		
	IEEE_754_FloatNum max;
	IEEE_754_FloatNum min;
	
	if(callID != SUB_CALL && x->sign == FALSE && y->sign == TRUE)
	{
		x->sign = y->sign = FALSE;
		paraSub(x, y, res, numProcs, procID);
		x->sign = FALSE;
		y->sign = TRUE;
		return;
	}
	
	if(callID != SUB_CALL && x->sign == TRUE && y->sign == FALSE)
	{
		x->sign = y->sign = FALSE;
		paraSub(y, x, res, numProcs, procID);
		x->sign = TRUE;
		y->sign = FALSE;
		return;		
	}
	
	carry = (UINT32 *)calloc(numProcs, sizeof(UINT32));
	localCarry = (UINT32 *)calloc(width, sizeof(UINT32));
	memset(localCarry, 0, width*sizeof(UINT32));

	min.significand = (UINT32 *)calloc(width, sizeof(UINT32));
	max.significand = (UINT32 *)calloc(width, sizeof(UINT32));

	//adjust exponent if necessary
	double exectime, start;

	
	if(x->exp != y->exp)
	{
		
		//get the smaller number
		minMax_IEEE(&min, &max, x, y, FALSE, numProcs, procID);
		makeSameExp(&min, &max, numProcs, procID);

	}
	else
	{
		memCPY(&min, x, width);
		memCPY(&max, y, width);
	}
	

	res->exp = max.exp;
	res->sign = x->sign;
	//1. add local parts of number	
	carry[procID] = addWithCarry(min.significand, 
								 max.significand, 
								 res->significand,
								 BASE,
								 width,
								 0);
	
	free(min.significand);
	free(max.significand);
	
	if(numProcs>1)
	{
		//2. exchange localCarry
		while(maxCarry > 0)
		{
			k++;
			maxCarry=0;
			MPI_Allgather (&carry[procID], 1, MPI_UNSIGNED, carry, 1, MPI_UNSIGNED, MPI_COMM_WORLD );
			
			for(j=1;j<numProcs;j++)
			{
				if(carry[j] > maxCarry)
					maxCarry = carry[j];
			}
			if(procID != 0)
				carry[procID] = FALSE;

			if(procID != (numProcs-1))//if proc does not owe last piece of number
			{
				if(carry[procID+1] != FALSE)
				{
					localCarry[width-1] = carry[procID+1];	
					carry[procID+1] = FALSE;
					
					//add the carries to the local part of the result
					carry[procID] += addWithCarry(res->significand, 
												  localCarry, 
												  res->significand,
												  BASE,
												  width,
												  0);
					
				}
			}		
			else
			{
				//set all carries to zero as they had been already added 
				//to the result, also from the last proc
				carry[numProcs-1] = FALSE;
			}
		}
		
		MPI_Allgather (&carry[procID], 1, MPI_UNSIGNED, carry, 1, MPI_UNSIGNED, MPI_COMM_WORLD );
		k++;
		//3. shift res to the right if neccessary
		if(carry[0] != FALSE && callID != SUB_CALL) //make the shift to the right
		{
			paraShiftRight(res, numProcs, procID);
			//add carry
			if(procID==0)
				res->significand[0] = carry[procID];
			k++;
		}
	}
	else //if only 1 processor is used!
	{
		//shift and add the carry
		if(carry[procID] != FALSE && callID != SUB_CALL)
		{			
			buf = res->significand[LOW_INDEX(procID, numProcs, res->precision)]; 
			for (j = width-1; j > 0; j--) 
			{
				res->significand[j] = res->significand[j-1];
			}
			
			res->significand[0] = carry[procID];
			res->exp++ ;
		}
	}
		
	//free memory	
	free(carry);
	free(localCarry);
	
}

//shift number one to the right, increase exponent, and add 1 to most significant digit.
void paraShiftRight(IEEE_754_FloatNum* x,
					UINT32 numProcs, 
					UINT32 procID){
	int lnbr, rnbr, i, j;
	UINT32 buf, buf2=0;
	MPI_Status status;
	MPI_Request	send_request, recv_request;
	
	UINT32 width = BLOCK_SIZE(procID, numProcs, x->precision);
	
	if (procID == 0) {
		lnbr = numProcs - 1;
		rnbr = procID + 1;
	} else if (procID == numProcs - 1) {//last processor
		lnbr = procID - 1;
		rnbr = 0; // loop the processors
	} else {
		lnbr = procID - 1;
		rnbr = procID + 1;
	}
	
	//set last number to zero to shift it in from the other side
	if(procID == (numProcs-1))
		x->significand[width-1] = 0;
	
	buf = x->significand[width-1]; //1st entry
	
	MPI_Isend(&buf,  1, MPI_UNSIGNED, rnbr, 10, MPI_COMM_WORLD, &send_request);
	MPI_Irecv(&buf2, 1, MPI_UNSIGNED, lnbr, 10, MPI_COMM_WORLD, &recv_request);
	
	i=width;
	for (j = 1; j < width; j++) {
		i--;
		x->significand[i] = x->significand[i-1];
	}
	
	MPI_Wait(&send_request, &status);
	MPI_Wait(&recv_request, &status);
	
	x->significand[0] = buf2; //last entry
		
	x->exp++; //increase the exp due to the shift
}

//shift number one to the left and subtract 1 from exponent.
void paraShiftLeft(IEEE_754_FloatNum* x,
					UINT32 numProcs, 
					UINT32 procID){
	int lnbr, rnbr, i, j;
	UINT32 buf, buf2=0;
	MPI_Status status;
	MPI_Request	send_request, recv_request;
	
	UINT32 width = BLOCK_SIZE(procID, numProcs, x->precision);
	
	if (procID == 0) {
		lnbr = numProcs - 1;
		rnbr = procID + 1;
	} else if (procID == numProcs - 1) {//last processor
		lnbr = procID - 1;
		rnbr = 0; // loop the processors
	} else {
		lnbr = procID - 1;
		rnbr = procID + 1;
	}
	
	buf = x->significand[0]; //1st entry
	
	MPI_Isend(&buf,  1, MPI_UNSIGNED, lnbr, 11, MPI_COMM_WORLD, &send_request);
	MPI_Irecv(&buf2, 1, MPI_UNSIGNED, rnbr, 11, MPI_COMM_WORLD, &recv_request);
	
	for (j = 0; j < width; j++) {
		x->significand[j] = x->significand[j+1];
	}
	MPI_Wait(&send_request, &status);
	MPI_Wait(&recv_request, &status);
	
	x->significand[width-1] = buf2; //last entry
	
	x->exp--; //decrease the exp due to the shift
}

//x- y = res
void paraSub(IEEE_754_FloatNum *x, 
			 IEEE_754_FloatNum *y,
			 IEEE_754_FloatNum *res,
			 UINT32 numProcs, 
			 UINT32 procID){
	BOOLEAN same=FALSE;
	UINT32 width = BLOCK_SIZE(procID, numProcs, x->precision);

	IEEE_754_FloatNum null;
	null.base	= BASE;
	null.exp	= 0;
	null.sign	= FALSE;
	null.significand = (UINT32 *)calloc(width, sizeof(UINT32));
	null.precision = x->precision;
	memset(null.significand, 0, width*sizeof(UINT32));
	
	if(equal_IEEE(x, &null, numProcs, procID) == TRUE &&
	   equal_IEEE(y, &null, numProcs, procID) == TRUE)
	{
		memCPY(res, &null, width);
		free(null.significand);
		return;
	}
	
	if(equal_IEEE(x, &null, numProcs, procID) == TRUE)
	{
		memCPY(res, y, width);
		res->sign = TRUE;
		free(null.significand);
		return;
	}
	
	if(equal_IEEE(y, &null, numProcs, procID) == TRUE)
	{
		memCPY(res, x, width);
		free(null.significand);
		return;
	}
	
	free(null.significand);

	//-x-(-y) = -x+y
	if(x->sign == TRUE &&
	   y->sign == TRUE)
	{
		x->sign = TRUE;
		y->sign = FALSE;
		paraAdd(x, y, res, numProcs, procID, ADD_CALL);
		x->sign = TRUE;
		y->sign = TRUE;
		return;
	}
	
	//x - (-y) = x+y
	if(x->sign == FALSE &&
	   y->sign == TRUE)
	{
		x->sign = FALSE;
		y->sign = FALSE;
		paraAdd(x, y, res, numProcs, procID, ADD_CALL);
		res->sign = FALSE;
		x->sign = FALSE;
		y->sign = TRUE;
		return;
	}
	
	//- x - y = -(x+y)
	if(x->sign == TRUE &&
	   y->sign == FALSE)
	{
		x->sign = FALSE;
		y->sign = FALSE;
		paraAdd(x, y, res, numProcs, procID, ADD_CALL);
		res->sign = TRUE;
		x->sign = TRUE;
		y->sign = FALSE;
		return;
	}
	
	IEEE_754_FloatNum max;
	max.significand = (UINT32 *)calloc(width, sizeof(UINT32));
	IEEE_754_FloatNum min;
	min.significand = (UINT32 *)calloc(width, sizeof(UINT32));
	
	//1. get the smaller number and obtain the sign of the result
	same = minMax_IEEE(&min, &max, x, y, TRUE, numProcs, procID); 

	if(same != TRUE)
	{
		//0. adjust if necessary
		if(x->exp != y->exp)
			makeSameExp(&min, &max, numProcs, procID);
			
		//2. make complement of res -> res.comp
		paraCompl(&min, numProcs, procID);

		//3. x+y.comp with paraAdd()
		paraAdd(&min, &max, res, numProcs, procID, SUB_CALL);

		//4. normalize result
		paraNorm(res, numProcs, procID);

		res->sign = min.sign;
	}
	else
	{
		memset(res->significand, 0, width*sizeof(UINT32));
		res->exp = 0;
		res->sign = FALSE;
	}
	
	free(min.significand);
	free(max.significand);

}

void paraCompl(IEEE_754_FloatNum *x,
			   UINT32 numProcs, 
			   UINT32 procID){
	UINT32 width = BLOCK_SIZE(procID, numProcs, x->precision);
	int i;
	for(i=0; i < width; i++)
	{
		x->significand[i] = (x->base - 1) - x->significand[i];		
	}
	if(procID==(numProcs-1))
		x->significand[width-1] += 1;
}

void paraMult(IEEE_754_FloatNum *x,
			  IEEE_754_FloatNum *y,
			  IEEE_754_FloatNum *res,
			  UINT32 numProcs, 
			  UINT32 procID)
{

	UINT32 width = BLOCK_SIZE(procID, numProcs, x->precision);
	
	IEEE_754_FloatNum null;
	null.base	= BASE;
	null.exp	= 0;
	null.sign	= FALSE;
	null.significand = (UINT32 *)calloc(width, sizeof(UINT32));
	null.precision = x->precision;
	memset(null.significand, 0, width*sizeof(UINT32));

	if(equal_IEEE(y, &null, numProcs, procID) == TRUE ||
	   equal_IEEE(x, &null, numProcs, procID) == TRUE)
	{
		memCPY(res, &null, width);
		free(null.significand);
		return;
	}
	
	free(null.significand);

	IEEE_754_FloatNum one;
	one.base	= BASE;
	one.exp	= 0;
	one.sign	= FALSE;
	one.significand = (UINT32 *)calloc(width, sizeof(UINT32));
	one.precision = x->precision;
	memset(one.significand, 0, width*sizeof(UINT32));	

	if(procID == 0)
	{
		one.significand[0]=1;
	}
	
	if(equal_IEEE(x, &one, numProcs, procID) == TRUE)
	{
		memCPY(res, y, width);
		free(one.significand);
		return;
	}
	
	if(equal_IEEE(y, &one, numProcs, procID) == TRUE)
	{
		memCPY(res, x, width);
		free(one.significand);
		return;
	}
	
	free(one.significand);
	
	UINT32 i, finExp;
	int j, jglob, p, s;
	p = numProcs;
	s = procID;
	double *compX = (double*)calloc(width*4, sizeof(double));
	double *compY = (double*)calloc(width*4, sizeof(double));

	finExp = x->exp + y->exp;
	
	//block to cyclic redistribution
	toCyclic(x->significand, compX, width, numProcs, procID);
	toCyclic(y->significand, compY, width, numProcs, procID);

	for(i=1;i<=width;i++)
	{
		compX[2*(width-i)]=compX[width-i];
		compX[2 * (width-i) + 1] = 0; 
		compX[2 * (i-1 + width)] = 0; 
		compX[2 * (i-1 + width) + 1] = 0;
		
		compY[2*(width-i)]=compY[width-i];
		compY[2 * (width-i) + 1] = 0; 
		compY[2 * (i-1 + width)] = 0; 
		compY[2 * (i-1 + width) + 1] = 0;
	}
	
	//reset res
	res->exp = 0;

	//perform FFT on numbers
	MPI_FFT(compX, x->precision*2, procID, numProcs, 1);
	MPI_FFT(compY, x->precision*2, procID, numProcs, 1);
	
	double reX, imX, reY, imY;
	for(i = 0; i < width*2; i++)
	{
		reX = compX[2*i+0];
		imX = compX[2*i+1];
		reY = compY[2*i+0];
		imY = compY[2*i+1];
		
		compX[2*i+0] = reX*reY-imX*imY;
		compX[2*i+1] = reX*imY+reY*imX;
	}
	
	//perform an inverse FFT
	MPI_FFT(compX, x->precision*2, procID, numProcs, -1);
	
	//remove imaginary part
	for(j=0, i=0; j < width*4; j=j+2, i++)
	{
        //jglob= j*numProcs+procID;
		compX[i] = compX[j];
	}


	//cyclic to block redistribution compX->compY, but just for the half of the array
	toBlock(compX, compY, numProcs, procID, 2*width);
	//compY now has the the right numbers and also reduced precision!
	IEEE_754_FloatNum temp;
	temp.base	= BASE;
	temp.exp	= 0;
	temp.sign	= FALSE;
	temp.significand = (UINT32 *)calloc(2*width, sizeof(UINT32));
	temp.precision = 2*x->precision;
	memset(temp.significand, 0, width*sizeof(UINT32));
	
	IEEE_754_FloatNum IEEE_Null = create_IEEE_NULL(2*x->precision, numProcs, procID);
	IEEE_754_FloatNum IEEE_Temp = create_IEEE_NULL(2*x->precision, numProcs, procID);

	for(i=0; i < 2*width; i++){
		IEEE_Temp.significand[i]=(int)(compY[i]+0.5);
	}
	
	IEEE_Temp.exp = 0;
	//carry add on the result
	paraAdd(&IEEE_Temp, &IEEE_Null, &temp, numProcs, procID, ADD_CALL);

	reducePrecision(&temp, res, x->precision, numProcs, procID);
	
	res->exp = finExp + res->exp;
	res->sign= x->sign | y->sign;
	
	//normalize res
	paraNorm(res, numProcs, procID);
   
	free(compX);
	free(compY);
	
	kill_IEEE_NULL(IEEE_Temp);
	kill_IEEE_NULL(IEEE_Null);
	kill_IEEE_NULL(temp);

}

//remove leading zeros
void paraNorm(IEEE_754_FloatNum* x,
			  UINT32 numProcs, 
			  UINT32 procID){
	
	//printf("prec in NORM=%d\n", x->precision);
	UINT32 width = BLOCK_SIZE(procID, numProcs, x->precision);
	
	UINT32 *n = (UINT32 *)calloc(numProcs, sizeof(UINT32));
	UINT32 *l = (UINT32 *)calloc(numProcs, sizeof(UINT32));
	memset(l, 0, numProcs*sizeof(UINT32));
	UINT32 zeros=0;
	UINT32 j, i = 0;
	//determine amout of leading zeros for each part of the number
	while(x->significand[i] == 0 && i < width)
	{
		n[procID]++;
		i++;
	}
	
	//broadcast amout of zeros
	MPI_Allgather(&n[procID], 1, MPI_UNSIGNED, n, 1, MPI_UNSIGNED, MPI_COMM_WORLD );
	
	if(n[0] > 0)
	{
		i = 0;
		while(n[i] == width && i < numProcs-1)
		{
			i++;
		}
		
		zeros = n[i] + i*width; //total amount of leading zeros
		
		//reset n[j] fuer j >= i
		if(i != numProcs-1 && numProcs>1)
		{
			i++;
			for(i; i<numProcs; i++)
				n[i]=0;
		}
		
		UINT32 arrayLength = (x->precision-zeros);
		UINT32 *array = (UINT32 *)calloc(arrayLength, sizeof(UINT32));
		memset(array, 0, arrayLength*sizeof(UINT32));
		
		//prepare array
		for(j=0;j<numProcs-1;j++)
		{
			//find out a place inside the array
			l[j+1]=width-n[j]+l[j];
		}
		
		for(j=n[procID], i=l[procID];  j<width && i<arrayLength; j++, i++)
		{
			array[i] = x->significand[j];
		}
		
		//all reduce using SUM
		MPI_Allreduce(array,array,arrayLength,MPI_UNSIGNED,MPI_SUM,MPI_COMM_WORLD);
		
		//clean what was there before
		memset(x->significand, 0, width*sizeof(UINT32));
		
		//redisribute the array
		for(i=0, j=procID*width; i<width && j<arrayLength; i++, j++)
		{
			x->significand[i] = array[j];
		}
		
		free(array);
	}
	
	free(n);
	free(l);
	x->exp -= zeros; //decrease the exp due to the shift
}

 
void paraDiv(IEEE_754_FloatNum *x, 
			  IEEE_754_FloatNum *y,
			  IEEE_754_FloatNum *res,
			  UINT32 numProcs, 
			  UINT32 procID)
{
	//division of large integers using Newton
	//1 / y = res
	UINT32 width = BLOCK_SIZE(procID, numProcs, x->precision);
	
	IEEE_754_FloatNum null;
	null.base	= BASE;
	null.exp	= 0;
	null.sign	= FALSE;
	null.significand = (UINT32 *)calloc(width, sizeof(UINT32));
	null.precision = x->precision;
	memset(null.significand, 0, width*sizeof(UINT32));
	

	if(equal_IEEE(x, &null, numProcs, procID) == TRUE)
	{
		memCPY(res, &null, width);
		free(null.significand);
		return;
	}
	
	int finExp = y->exp;
	res->exp = 0;
	
	IEEE_754_FloatNum temp1;
	temp1.base	= BASE;
	temp1.exp	= -1;
	temp1.sign	= FALSE;
	temp1.significand = (UINT32 *)calloc(width, sizeof(UINT32));
	temp1.precision = x->precision;
	memset(temp1.significand, 0, width*sizeof(UINT32));
	
	IEEE_754_FloatNum temp2;
	temp2.base	= BASE;
	temp2.exp	= 0;
	temp2.sign	= FALSE;
	temp2.significand = (UINT32 *)calloc(width, sizeof(UINT32));
	temp2.precision = x->precision;
	memset(temp2.significand, 0, width*sizeof(UINT32));

	
	IEEE_754_FloatNum temp3;
	temp3.base	= BASE;
	temp3.exp	= 0;
	temp3.sign	= FALSE;
	temp3.significand = (UINT32 *)calloc(width, sizeof(UINT32));
	temp3.precision = x->precision;
	memset(temp3.significand, 0, width*sizeof(UINT32));

	
	UINT32 iter = (int)log2((y->precision+1)/log2(17))*7;
	y->exp=0;
	UINT32 init = 0;
	double u = 0;
	double u_max, u_min;
	double min=0;
	double max=0;
	
	if(procID == 0)
	{
		if(y->significand[0]>1)
			min=1;
		if(y->significand[0]<(BASE-1))
			max=1;
		u_min = (double)y->significand[0]+(double)y->significand[1]/(double)BASE-min;
		u_max = (double)y->significand[0]+(double)y->significand[1]/(double)BASE+max;
		//u = (0.5*(1.0/u_max+1.0/u_min)*(double)BASE);
		u = 1.0/((double)y->significand[0]+(double)y->significand[1]/(double)BASE);
		init = (int)u;
		temp1.significand[0] = init; 	//put init in number
		temp1.significand[1] = (int)((u - (double)init)*(double)BASE); 
		temp2.significand[0] = 1;
	}
	int i, j=0;
	
	for(i=1; i<=ITERS*4; i++)
	{
		paraPow2(&temp1, &temp3, numProcs,procID); //xk^2 = temp3
		paraMult(&temp3, y, res, numProcs, procID); //b*xk^2 = res
		paraSub(&temp1, res, &temp3, numProcs, procID); //xk-b*xk^2 = temp3
		paraAdd(&temp3, &temp1, &temp1, numProcs, procID, ADD_CALL);//

		//validation
		paraMult(&temp1, y, &temp3, numProcs, procID); 
		paraSub(&temp3, &temp2, &temp3, numProcs, procID);

		if((temp3.exp <= (int)(-1.0*((double)x->precision*0.995-0.5))) || 
		   (i > 0 && equal_IEEE(&temp3, &null, numProcs, procID) == TRUE))
		{
			break;
		}
		
		//if divergent start over with another init value
		if(equal_IEEE(&temp3, &null, numProcs, procID) != TRUE && j == 0 && i > 4 && temp3.exp > -1)
		{
			i=0;
			j++;
			memset(temp1.significand, 0, width*sizeof(UINT32));
			
			if(procID == 0)
			{
				temp1.significand[0] = 1;
			}

			temp1.exp  = -16; //machine precision
			y->exp = finExp;
			temp1.sign = FALSE;	
		}

	}

	temp1.exp=0;
	int xTempExp = x->exp;
	x->exp = y->exp = 0;
	BOOLEAN xSign = x->sign;
	BOOLEAN ySign = y->sign;

	x->sign = y->sign = FALSE;
	if(smaller_IEEE(x, y, numProcs, procID) == TRUE)
	{	
		temp1.exp = -1;
	}
	
	y->exp = finExp;
	x->exp = xTempExp;
	finExp = x->exp-finExp+temp1.exp;

	res->exp=0;
	temp1.exp=0;
	memset(res->significand, 0, width*sizeof(UINT32));

	//correct by error
	paraSub(&temp1, &temp3, &temp1, numProcs, procID);

	//finally mult x*res = res
	paraMult(x, &temp1, res, numProcs, procID); //D*x^2 = temp3
	
	res->exp = finExp;
	x->sign = xSign;
	y->sign = ySign;
	res->sign = x->sign | y->sign;
	
	free(temp1.significand);
	free(temp2.significand);
	free(temp3.significand);
	free(null.significand);

}

void makeSameExp(IEEE_754_FloatNum *min, 
				 IEEE_754_FloatNum *max,
				 UINT32 numProcs, 
				 UINT32 procID)
{

	INT32 shift = abs(abs(min->exp)-abs(max->exp));
	UINT32 width = BLOCK_SIZE(procID, numProcs, min->precision);

	int j, k, i;
	
	if(shift >= min->precision)
	{
		memset(min->significand, 0, width*sizeof(UINT32));
		min->exp = 0;
	}
	else
	{
		UINT32 *n = (UINT32 *)calloc(numProcs, sizeof(UINT32));
		UINT32 *l = (UINT32 *)calloc(numProcs, sizeof(UINT32));
		UINT32 arrayLength = (min->precision-shift);
		UINT32 *array = (UINT32 *)calloc(arrayLength, sizeof(UINT32));
		memset(array, 0, arrayLength*sizeof(UINT32));
		
		//prepare array
		k=0;
		for(i=0; i<numProcs && k<shift; i++)
		{
			if(k<(shift-shift%width))
			{
				n[i]=width;
				k+=width;
			}
			else
			{
				n[i]=shift%width;
				k+=shift%width;
			}
			
		}
		
		for(i=procID*width, j=0; i < arrayLength && j < width; i++, j++)
		{
			array[i] = min->significand[j];
		}
		
		MPI_Allreduce(array,array,arrayLength,MPI_UNSIGNED,MPI_SUM,MPI_COMM_WORLD);

		//clean what was there before
		memset(min->significand, 0, width*sizeof(UINT32));
		
		//redisribute the array
		//prepare array
		for(j=0;j<numProcs-1;j++)
		{
			//find out a place inside the array
			l[j+1]=width-n[j]+l[j];
		}
		for(i=n[procID], j=l[procID]; i<width && j<arrayLength; i++, j++)
		{
			min->significand[i] = array[j];
		}
		
		min->exp += shift;
		free(array);
		free(n);
		free(l);
	}
		
}

void paraSqrt(IEEE_754_FloatNum *x,
			 IEEE_754_FloatNum *res,
			 UINT32 numProcs, 
			 UINT32 procID,
			 int n)
{
	//sqrt of large integers using Newton
	//1 / sqrt(y) = res
	UINT32 width = BLOCK_SIZE(procID, numProcs, x->precision);
	if(n<=1)
		res->exp = 0;
	
	IEEE_754_FloatNum temp1;
	temp1.base	= BASE;
	temp1.exp	= -1;
	temp1.sign	= FALSE;
	temp1.significand = (UINT32 *)calloc(width, sizeof(UINT32));
	temp1.precision = x->precision;
	memset(temp1.significand, 0, width*sizeof(UINT32));
	
	IEEE_754_FloatNum temp2;
	temp2.base	= BASE;
	temp2.exp	= 0;
	temp2.sign	= FALSE;
	temp2.significand = (UINT32 *)calloc(width, sizeof(UINT32));
	temp2.precision = x->precision;
	memset(temp2.significand, 0, width*sizeof(UINT32));
	
	IEEE_754_FloatNum temp3;
	temp3.base	= BASE;
	temp3.exp	= 0;
	temp3.sign	= FALSE;
	temp3.significand = (UINT32 *)calloc(width, sizeof(UINT32));
	temp3.precision = x->precision;
	memset(temp3.significand, 0, width*sizeof(UINT32));
	
	IEEE_754_FloatNum oneHalf;
	oneHalf.base	= BASE;
	oneHalf.exp		= -1;
	oneHalf.sign	= FALSE;
	oneHalf.significand = (UINT32 *)calloc(width, sizeof(UINT32));
	oneHalf.precision = x->precision;
	memset(oneHalf.significand, 0, width*sizeof(UINT32));
	
	UINT32 iter = (int)log2((x->precision+1)/log2(17))*7;
	
	int expSign = (x->exp >= 0) ?  1 :  -1;
	
	int finExp = (x->exp-(x->exp%2)*expSign)/2;
	int xExp = x->exp;
	x->exp=(x->exp%2);

	UINT32 init = 0;
	double u = 0;
	double low=0;
	double high=0;

	if(procID == 0)
	{
		if(x->significand[0]>1)
			low=1;
		if(x->significand[0]<(BASE-1))
			high=1;
		u = (0.5*(sqrt((double)x->significand[0]-low)+sqrt((double)x->significand[0]+high))*(double)1);
		init = (int)u;
		
		if(n==1)
		temp1.significand[0] = init*BASE/10; 	//put init in number

		oneHalf.significand[0]=5*BASE/10;
		temp2.significand[0] = 1;
	}
	if(n>1)
		memCPY(&temp1, res, width);

	int i,j=0;
	
	for(i=0; i<=ITERS*2; i++)
	{
		
		paraMult(&temp1, &temp1, &temp3, numProcs, procID); //x^2 = temp3
		paraMult(x, &temp3, res, numProcs, procID); //D*x^2 = res		
		paraSub(&temp2, res, &temp3, numProcs, procID); //1-D*x^2 = temp3		
		paraMult(&temp1, &oneHalf, res, numProcs, procID); //x/2 = res
		paraMult(res, &temp3, &temp3, numProcs, procID); //x/2*() = temp3		
		paraAdd(&temp1, &temp3, &temp1, numProcs, procID, ADD_CALL);//x+x/2*() = temp1		
		//------validation
		paraMult(x, &temp1, res, numProcs, procID);
		paraMult(res, res, &temp3, numProcs, procID);	
		paraSub(&temp3, x, &temp3, numProcs, procID); 
		if(temp3.exp < (int)(-1.0*((double)x->precision*0.95-0.5)))
		{
			break;
		}

		
	}
#ifdef TEST
	if(procID==0)
	printf("SQRT iters: i=%d\n", i);
	print_IEEE_754_FloatNum(&temp3, "ERROR SQRT", numProcs, procID); 
#endif	
	x->exp = xExp;
	res->exp = finExp;


	free(temp1.significand);
	free(temp2.significand);
	free(temp3.significand);
	free(oneHalf.significand);
}

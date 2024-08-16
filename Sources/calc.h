#ifndef __CALC_H_
#define __CALC_H_

#include "helper.h"

void paraMultBy2(IEEE_754_FloatNum *x,
				 IEEE_754_FloatNum *res,
				 UINT32 numProcs, 
				 UINT32 procID);

void paraInvMod(IEEE_754_FloatNum *x, 
				IEEE_754_FloatNum *modulus,
				IEEE_754_FloatNum *res,
				UINT32 numProcs, 
				UINT32 procID);

void paraExtendedEuclidean(IEEE_754_FloatNum *a,
						   IEEE_754_FloatNum *b,
						   IEEE_754_FloatNum *x, 
						   IEEE_754_FloatNum *y,
						   UINT32 numProcs, 
						   UINT32 procID);


void paraFloor(IEEE_754_FloatNum *x,
			   IEEE_754_FloatNum *res,
			   UINT32 numProcs, 
			   UINT32 procID);
void paraPow2(IEEE_754_FloatNum *x,
			  IEEE_754_FloatNum *res,
			  UINT32 numProcs, 
			  UINT32 procID);
UINT32 paraDivBy2(IEEE_754_FloatNum *x,
				  IEEE_754_FloatNum *res,
				  UINT32 numProcs, 
				  UINT32 procID);
void paraRand(IEEE_754_FloatNum *x, 
			  UINT32 numProcs, 
			  UINT32 procID,
			  UINT32 length);
void paraMod(IEEE_754_FloatNum *x, 
			 IEEE_754_FloatNum *y,
			 IEEE_754_FloatNum *res,
			 UINT32 numProcs, 
			 UINT32 procID);
void paraMod2(IEEE_754_FloatNum *x, 
			 IEEE_754_FloatNum *y,
			 IEEE_754_FloatNum *res,
			 UINT32 numProcs, 
			 UINT32 procID);
void paraModExp(IEEE_754_FloatNum *base, 
				IEEE_754_FloatNum *exp,
				IEEE_754_FloatNum *modulus,
				IEEE_754_FloatNum *res,
				UINT32 numProcs, 
				UINT32 procID);

void paraSub(IEEE_754_FloatNum *x, 
			 IEEE_754_FloatNum *y,
			 IEEE_754_FloatNum *res,
			 UINT32 numProcs, 
			 UINT32 procID);

void paraAdd(IEEE_754_FloatNum *x, 
				IEEE_754_FloatNum *y,
				IEEE_754_FloatNum* res,
				UINT32 numProcs, 
				UINT32 procID,
				BOOLEAN callID);


UINT32 addWithCarry(UINT32 *x, 
					 UINT32 *y, 
					 UINT32 *res,
					 UINT32 base,
					 UINT32 width,
					 UINT32 offset);

void paraMult(IEEE_754_FloatNum *x, 
			  IEEE_754_FloatNum *y,
			  IEEE_754_FloatNum *res,
			  UINT32 numProcs, 
			  UINT32 procID);

void paraDiv(IEEE_754_FloatNum *x, 
			 IEEE_754_FloatNum *y,
			 IEEE_754_FloatNum *res,
			 UINT32 numProcs, 
			 UINT32 procID);

void paraSqrt(IEEE_754_FloatNum *x,
			  IEEE_754_FloatNum *res,
			  UINT32 numProcs, 
			  UINT32 procID,
			  int n);

void paraGCD(IEEE_754_FloatNum *x, 
			 IEEE_754_FloatNum *y,
			 IEEE_754_FloatNum *res,
			 UINT32 numProcs, 
			 UINT32 procID);


void powOf2(UINT32 p,
			  IEEE_754_FloatNum *res,
			  UINT32 numProcs, 
			  UINT32 procID);

void paraShiftRight(IEEE_754_FloatNum* x,
					UINT32 numProcs, 
					UINT32 procID);

void paraShiftLeft(IEEE_754_FloatNum* x,
				   UINT32 numProcs, 
				   UINT32 procID);

void paraCompl(IEEE_754_FloatNum *res,
			   UINT32 numProcs, 
			   UINT32 procID);

void paraNorm(IEEE_754_FloatNum* x,
			  UINT32 numProcs, 
			  UINT32 procID);

void makeSameExp(IEEE_754_FloatNum *min, 
				 IEEE_754_FloatNum *max,
				 UINT32 numProcs, 
				 UINT32 procID);

#endif

#ifndef __RSA_H__
#define __RSA_H__

void RSA_genKey(UINT32 procID, UINT32 numProcs, UINT32 keylength);
void RSA_encrypt(IEEE_754_FloatNum *m, 
				 IEEE_754_FloatNum *c,  
				 IEEE_754_FloatNum* N, 
				 IEEE_754_FloatNum* E,  
				 UINT32 procID, 
				 UINT32 numProcs );

void RSA_decrypt(IEEE_754_FloatNum *c, 
				 IEEE_754_FloatNum *m,  
				 IEEE_754_FloatNum *Q,
				 IEEE_754_FloatNum *P,
				 IEEE_754_FloatNum *QP,
				 IEEE_754_FloatNum *DP,
				 IEEE_754_FloatNum *DQ,
				 UINT32 procID, 
				 UINT32 numProcs );

void read_RSA_context(IEEE_754_FloatNum* E,
					  IEEE_754_FloatNum* P,
					  IEEE_754_FloatNum* Q,
					  IEEE_754_FloatNum* N,
					  IEEE_754_FloatNum* D,
					  IEEE_754_FloatNum* DP,
					  IEEE_754_FloatNum* DQ,
					  IEEE_754_FloatNum* QP,
					  UINT32 procID, 
					  UINT32 numProcs);
#endif     

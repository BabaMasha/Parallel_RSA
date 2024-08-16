#include <stdio.h>
#include <stdlib.h>

#include "helper.h"
#include "mpi.h"
void increasePrecision(IEEE_754_FloatNum *oldPrecNum,
					   IEEE_754_FloatNum *newPrecNum,
					   UINT32 newPRecision,
					   UINT32 numProcs, 
					   UINT32 procID)
{
	UINT32 width = BLOCK_SIZE(procID, numProcs, oldPrecNum->precision);
	UINT32 *array = (UINT32 *)calloc(newPrecNum->precision, sizeof(UINT32));
	memset(array, 0, oldPrecNum->precision);
	int i;
	MPI_Allgather(oldPrecNum->significand, width, MPI_UNSIGNED, array, width, MPI_UNSIGNED, MPI_COMM_WORLD);
	
	width = BLOCK_SIZE(procID, numProcs, newPrecNum->precision);
	MPI_Scatter( array, width,  MPI_UNSIGNED, newPrecNum->significand, width,  MPI_UNSIGNED, 0, MPI_COMM_WORLD); 
	newPrecNum->sign = oldPrecNum->sign;
	newPrecNum->exp = oldPrecNum->exp;
	newPrecNum->precision = newPRecision;
	free(array);	
}

void reducePrecision(IEEE_754_FloatNum *oldPrecNum,
					 IEEE_754_FloatNum *newPrecNum,
					 UINT32 newPRecision,
					 UINT32 numProcs, 
					 UINT32 procID)
{
	UINT32 width = BLOCK_SIZE(procID, numProcs, oldPrecNum->precision);
	UINT32 *array = (UINT32 *)calloc(oldPrecNum->precision, sizeof(UINT32));
	memset(array, 0, oldPrecNum->precision);
	int i;
	MPI_Allgather(oldPrecNum->significand, width, MPI_UNSIGNED, array, width, MPI_UNSIGNED, MPI_COMM_WORLD);
	
	width = BLOCK_SIZE(procID, numProcs, newPrecNum->precision);
	MPI_Scatter( array, width,  MPI_UNSIGNED, newPrecNum->significand, width,  MPI_UNSIGNED, 0, MPI_COMM_WORLD); 
	newPrecNum->sign = oldPrecNum->sign;
	newPrecNum->exp = oldPrecNum->exp;
	newPrecNum->precision = newPRecision;
	free(array);

}
void toIEEE_754(UINT32 intNumber,
				IEEE_754_FloatNum *x,
				UINT32 numProcs, 
				UINT32 procID)
{
	UINT32 width = BLOCK_SIZE(procID, numProcs, x->precision);
	UINT32 *array = (UINT32 *)calloc(x->precision, sizeof(UINT32));
	memset(array, 0, x->precision);
	
	unsigned int i, div, res, digs;
	digs = 0;
	div = BASE;
	//determine # of digits
	res = intNumber;
	while(res>0)
	{
		res = (int)intNumber/div;
		div *= BASE;
		digs++;
	}
	
	MPI_Barrier(MPI_COMM_WORLD);

	x->exp = digs-1;
	x->sign = FALSE;

	//put number in array
	int digsIter=digs-1;
	for(i=0; i<digs; i++)
	{
		array[i] =  (UINT32)(intNumber/(int)pow(BASE,digsIter));
		intNumber -= array[i]*(int)pow(BASE,digsIter);
		digsIter--;
	}
	//distribute data over processors
	MPI_Scatter( array, width,  MPI_UNSIGNED, x->significand, width,  MPI_UNSIGNED, 0, MPI_COMM_WORLD); 
	free(array);
}


char* makeFileName(char *input, int keyLength)
{
	char* first = "KEYS/";
	char middle[1];
	int i = keyLength;
	sprintf(middle,"%d", i);
	char* second = ".txt";
	char* both = malloc(strlen(first) + strlen(input) + strlen(middle) + strlen(second) + 1);
	strcpy(both, first);
	strcat(both, input);
	strcat(both, middle);
	strcat(both, second);
	return both;
}

int min( int a, int b){
    if( a < b ) return a;
    return b;
}

int max( int a, int b){
    if( a > b ) return a;
    return b;
}

BOOLEAN isODD(IEEE_754_FloatNum *x,
			  UINT32 numProcs, 
			  UINT32 procID)
{
	UINT32 width = BLOCK_SIZE(procID, numProcs, x->precision);
	BOOLEAN rest = FALSE;
	UINT32 digs = x->exp;
	
	UINT32 i, j, k;
	for(i=0, j=procID*width; i<width && j<digs; i++, j++)
	{
		;
	}
	
	
	//if x uneven set rest to 1
	if(i < width && j<=digs)
	{
		if(x->significand[i]%2 != 0)
		{
			rest = TRUE;
		}
	}
	
	MPI_Allreduce(&rest,&rest,1,MPI_UNSIGNED,MPI_SUM,MPI_COMM_WORLD);
	
	return rest;
}


void memCPY(IEEE_754_FloatNum *a, IEEE_754_FloatNum *b, 
			UINT32 width)
{
	//same as "void * memcpy ( void * destination, const void * source, size_t num );"
	memcpy(a->significand, b->significand, width*sizeof(UINT32));
	a->sign		= b->sign;
	a->base		= b->base; 
	a->exp		= b->exp; 
	a->precision= b->precision;
	
}

//a == b
BOOLEAN equal_IEEE(IEEE_754_FloatNum *a, IEEE_754_FloatNum *b,
				   UINT32 numProcs, 
				   UINT32 procID)
{
	UINT32 i=0;
	UINT32 width = BLOCK_SIZE(procID, numProcs, a->precision);

    if( a->exp != b->exp )
	{
		return FALSE;
	}
	
	if( a->sign != b->sign )
	{
		return FALSE;
	}
	
	//determine amout of same digits
	while(a->significand[i] == b->significand[i] && i < width-1)
	{
		i++;
	}
	
	i++;
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Allreduce(&i,&i,1,MPI_UNSIGNED,MPI_SUM,MPI_COMM_WORLD);
	
	if(i != a->precision)
	{
		return FALSE;
	}
	else
	{
		return TRUE;
	}
	
}

BOOLEAN smallerOne_IEEE(IEEE_754_FloatNum *a, UINT32 numProcs, UINT32 procID)
{
	UINT32 firstDig = 0;
	// if a->sign == TRUE
	if( a->sign == TRUE )
	{
		return FALSE;
	}
	
	if( a->exp < 0 && a->sign == FALSE )
	{
		//a < b
		return TRUE;
	}
	
	if( a->exp > 0 && a->sign == FALSE )
	{
		//a > b
		return FALSE;
	}
	
	if(a->exp == 0 && a->sign == FALSE)
	{
		if(procID==0)
		{
			firstDig = a->significand[0];
		}
		
		MPI_Bcast( &firstDig, 1,  MPI_UNSIGNED, 0, MPI_COMM_WORLD);
		if(firstDig == 0)
		{
			return TRUE;
		}
		else
		{
			return FALSE;
		}
	}
}


//a < b
BOOLEAN smaller_IEEE(IEEE_754_FloatNum *a, IEEE_754_FloatNum *b,
				   UINT32 numProcs, 
				   UINT32 procID)
{
	UINT32 i=0;
	UINT32 width = BLOCK_SIZE(procID, numProcs, a->precision);
	UINT32 *arrayA;
	UINT32 *arrayB;
		
    if( a->exp < b->exp && a->sign == b->sign )
	{
		//a < b
		return TRUE;
	}
	
	if( a->exp > b->exp && a->sign == b->sign )
	{
		//a > b
		return FALSE;
	}
	
	// if a->sign == FALSE && b->sign == TRUE
	if( a->sign < b->sign )
	{
		return FALSE;
	}
	
	// if a->sign == TRUE && b->sign == FALSE
	if( a->sign > b->sign )
	{
		return TRUE;
	}
	
	if(a->exp == b->exp && a->sign == b->sign)
	{
		arrayA = (UINT32 *)calloc(numProcs, sizeof(UINT32));
		arrayB = (UINT32 *)calloc(numProcs, sizeof(UINT32));
		
		//determine amout of same digits
		while(a->significand[i] == b->significand[i] && i < width-1)
		{
			i++;
		}

		//broadcast digit which is different or last digit which is same
		MPI_Allgather(&a->significand[i], 1, MPI_UNSIGNED, arrayA, 1, MPI_UNSIGNED, MPI_COMM_WORLD );
		MPI_Allgather(&b->significand[i], 1, MPI_UNSIGNED, arrayB, 1, MPI_UNSIGNED, MPI_COMM_WORLD );
		
		int j=0;
		while(arrayA[j] == arrayB[j] && j < numProcs-1)
		{
			j++;
		}		
		
		if( arrayA[j] < arrayB[j] )
		{ 
			//*min = *a;
			free(arrayA);
			free(arrayB);
			return TRUE;
		}
		else 
		{
			if( arrayA[j] > arrayB[j] )
			{ 
				//*min = *b;
				free(arrayA);
				free(arrayB);
				return FALSE;
			}
			
			if(j == (numProcs-1)) //numbers are same
			{
				//*min = *b;
				free(arrayA);
				free(arrayB);
				return FALSE;
			}
		}		
	}
	
	//*min = *a;
	free(arrayA);
	free(arrayB);
	return TRUE;
	
}


BOOLEAN minMax_IEEE(IEEE_754_FloatNum *min, IEEE_754_FloatNum *max,
				 IEEE_754_FloatNum *a, IEEE_754_FloatNum *b, 
				 BOOLEAN sign,
				 UINT32 numProcs, 
				 UINT32 procID){
	int i=0;
	UINT32 width = BLOCK_SIZE(procID, numProcs, a->precision);

	UINT32 *arrayA = (UINT32 *)calloc(numProcs, sizeof(UINT32));
	UINT32 *arrayB = (UINT32 *)calloc(numProcs, sizeof(UINT32));

	//obtain the min and change sign accordingly
    if( a->exp < b->exp )
	{
		//*min = *a;
		memCPY(min, a, width);
		memCPY(max, b, width);
		if(sign)
			min->sign = TRUE;
		
		free(arrayA);
		free(arrayB);
		return FALSE;
	}
	
	if( a->exp > b->exp)
	{
		//*min = *b;
		memCPY(min, b, width);
		memCPY(max, a, width);
		free(arrayA);
		free(arrayB);
		return FALSE;
	}
	
	if(a->exp == b->exp) 
	{

		//determine amout of same digits
		while(a->significand[i] == b->significand[i] && i < width-1)
		{
			i++;
		}
		//broadcast digit which is different or last digit which is same
		MPI_Allgather(&a->significand[i], 1, MPI_UNSIGNED, arrayA, 1, MPI_UNSIGNED, MPI_COMM_WORLD );
		MPI_Allgather(&b->significand[i], 1, MPI_UNSIGNED, arrayB, 1, MPI_UNSIGNED, MPI_COMM_WORLD );

		int j=0;
		while(arrayA[j] == arrayB[j] && j < numProcs-1)
		{
			j++;
		}
		

		if( arrayA[j] < arrayB[j] )
		{ 
			//*min = *a;
			memCPY(min, a, width);
			memCPY(max, b, width);
			if(sign)
				min->sign = TRUE;
			
			free(arrayA);
			free(arrayB);
			return FALSE;
		}
		else 
		{
			if( arrayA[j] > arrayB[j] )
			{ 
				//*min = *b;
				memCPY(min, b, width);
				memCPY(max, a, width);
				free(arrayA);
				free(arrayB);
				return FALSE;
			}
			
			if(j == (numProcs-1)) //numbers are same
			{
				//*min = *b;
				memCPY(min, b, width);
				memCPY(max, a, width);
				free(arrayA);
				free(arrayB);
				return TRUE;
			}
		}		
	}

	memCPY(min, a, width); 
	memCPY(max, b, width); 

	free(arrayA);
	free(arrayB);
	return FALSE;
}

void comparePI(IEEE_754_FloatNum *number, IEEE_754_FloatNum *number2, UINT32 numProcs, UINT32 procID)
{
	UINT32 *array = (UINT32 *)calloc(numProcs, sizeof(UINT32));
	UINT32 width = BLOCK_SIZE(procID, numProcs, number->precision);

	UINT32 i, j=0;
	for(i=0;i<width;i++)
	{
		if(number->significand[i] == number2->significand[i])
			j++;
		else
			break;
	}
	array[procID]=j;
	MPI_Allgather (&array[procID], 1, MPI_UNSIGNED, array, 1, MPI_UNSIGNED, MPI_COMM_WORLD );
	j=0;
	for(i=0;i<numProcs;i++)
	{
		j+=array[i];
	}
	
	free(array);

}

void read_IEEE_754_FloatNum(IEEE_754_FloatNum *number, UINT32 numProcs, UINT32 procID, char *fileName)
{
	UINT32 i,j,m;
	unsigned int k;
	UINT32 width = BLOCK_SIZE(procID, numProcs, number->precision);
	const unsigned int size = (unsigned int)log10((double)BASE);
	UINT32 *array = (UINT32 *)calloc(number->precision, sizeof(UINT32));
	memset(array, 0, number->precision);
	FILE *fh = NULL;
	char drei=0;
	char num=0;
	char CR=0;
	if(procID == 0)
	{
		fh = fopen( fileName, "rt" );
		while( ! feof(fh) )
		{
			fscanf (fh, "%c", &CR);
			if(CR == '\n')
				break;
			number->exp *= 10;
			number->exp += (INT32)(CR -'0');

		}
		
		
		fscanf (fh, "%c", &drei);
		k=0;
		
		while(k<(number->precision))
		{
			
			j=BASE;
			if( ! feof(fh) )
			{
				for(i = 0; i < size; i++)
				{
					fscanf (fh, "%c", &num);
					j=j/10;
					
					array[k] += j * (UINT32)(num-'0');
			
				}
			}
			else
			{
				array[k] = 0;
			}
			
			k++;
			
		}
		fclose(fh);
		m = (UINT32)(drei-'0');
	}

	MPI_Bcast( &m, 1,  MPI_UNSIGNED, 0, MPI_COMM_WORLD);
	MPI_Bcast( &number->exp, 1,  MPI_UNSIGNED, 0, MPI_COMM_WORLD);

	if(m>0)
		number->sign = TRUE;
	else
		number->sign = FALSE;

	MPI_Scatter( array, width,  MPI_UNSIGNED, number->significand, width,  MPI_UNSIGNED, 0, MPI_COMM_WORLD); 
	
	free(array);
}

void write_IEEE_754_FloatNum(IEEE_754_FloatNum *number, UINT32 numProcs, UINT32 procID, char *fileName)
{
	unsigned int k;
	UINT32 width = BLOCK_SIZE(procID, numProcs, number->precision);
	int numLength = number->exp;
	UINT32 *array = (UINT32 *)calloc(number->precision, sizeof(UINT32));
	memset(array, 0, number->precision);
	FILE *fh = NULL;
	
	MPI_Gather(number->significand, width, MPI_UNSIGNED, array, width, MPI_UNSIGNED, 0, MPI_COMM_WORLD );

	if(procID == 0)
	{
		fh = fopen( fileName, "wt" );
		fprintf (fh, "%d\n", number->exp); 
		fprintf (fh, "%d", number->sign); 
		k=0;
				
		while(k <= numLength)
		{
			
			if( ! feof(fh) )
			{
				fprintf (fh, "%d", array[k]); 
			}
			else
			{
				break;
			}
			
			k++;
			
		}

		fprintf (fh, "%d", 0); 
		fprintf (fh, "%d", 0); 
		fclose(fh);
	}
	
	free(array);
}

void read_PI_1000000(IEEE_754_FloatNum *number, UINT32 numProcs, UINT32 procID){
	UINT32 i,j,m;
	unsigned int k;
	UINT32 width = BLOCK_SIZE(procID, numProcs, number->precision);
	
	UINT32 *array = (UINT32 *)calloc(number->precision, sizeof(UINT32));
	memset(array, 0, number->precision);
	const char *szFileName = "pi.dat";
	FILE *fh = NULL;
    unsigned int nLine = 0;
	const unsigned int size = (unsigned int)log10((double)BASE);
	char drei=0;
	char num=0;
	char * pEnd;
	if(procID == 0)
	{
	fh = fopen( szFileName, "rt" );
	
	k=0;

	fscanf (fh, "%c", &drei);
	k=1;
		
	while( ! feof(fh) && k<number->precision)
    {
	
		j=BASE;
		for(i = 0; i < size; i++)
		{
			fscanf (fh, "%c", &num);
			j=j/10;
			
			array[k] += j * (UINT32)(num-'0');

		}
		
		k++;

	}
	fclose(fh);
	array[0]=(UINT32)(drei-'0');
	}
	MPI_Scatter( array, width,  MPI_UNSIGNED, number->significand, width,  MPI_UNSIGNED, 0, MPI_COMM_WORLD); 
	
	free(array);
	number->sign = FALSE;
	number->exp = 0;
	
}

void print_IEEE_754_FloatNum(IEEE_754_FloatNum *number, char *name, UINT32 numProcs, UINT32 procID){
	UINT32 i,j,k,num;
	UINT32 width = BLOCK_SIZE(procID, numProcs, number->precision);
	UINT32 *array = (UINT32 *)calloc(number->precision, sizeof(UINT32));

	char *sign = "";
	MPI_Allgather(number->significand, width, MPI_UNSIGNED, array, width, MPI_UNSIGNED, MPI_COMM_WORLD );
	UINT32 prec = number->precision;
	if(prec < 64)
		prec = number->precision;
	else
		prec = 64;
	
	if(procID == 0)
	{
		if(number->sign == TRUE)
			sign = "-";
	
		printf("%s %s", name, sign);
	
		for(i = 0; i < prec; i++)
		{
			if(i == 1)
				printf(".");
			num=array[i];
			for(j=BASE; j>=10;)
			{ 
				j=j/10;
				printf("%d", (num)/j);
				k=(num)/j;
				num=num-k*j;
				
			}

		}
		printf("*%d^%d\n", number->base, number->exp);
	}
	
	free(array);
}

/* These functions can be used to allocate and deallocate vectors and matrices.
 If not enough memory available, one processor halts them all.
 */

double *vecallocd(int n){
    /* This function allocates a vector of doubles of length n */
    double *pd;
	
    if (n==0){
        pd= NULL;
    } else {
        pd= (double *)malloc(n*SZDBL);
        if (pd==NULL)
            MPI_Abort(MPI_COMM_WORLD,-2);
    }
    return pd;
	
} /* end vecallocd */

int *vecalloci(int n){
    /* This function allocates a vector of integers of length n */
    int *pi;
	
    if (n==0){
        pi= NULL; 
    } else { 
        pi= (int *)malloc(n*SZINT);
        if (pi==NULL)
            MPI_Abort(MPI_COMM_WORLD,-3);
    }
    return pi;
	
} /* end vecalloci */

double **matallocd(int m, int n){
    /* This function allocates an m x n matrix of doubles */
    int i;
    double *pd, **ppd;
	
    if (m==0){
        ppd= NULL;  
    } else { 
        ppd= (double **)malloc(m*sizeof(double *));
        if (ppd==NULL)
            MPI_Abort(MPI_COMM_WORLD,-4);
        if (n==0){
            for (i=0; i<m; i++)
                ppd[i]= NULL;
        } else {  
            pd= (double *)malloc(m*n*SZDBL); 
            if (pd==NULL)
                MPI_Abort(MPI_COMM_WORLD,-4);
            ppd[0]=pd;
            for (i=1; i<m; i++)
                ppd[i]= ppd[i-1]+n;
        }
    }
    return ppd;
	
} /* end matallocd */

void vecfreed(double *pd){
    /* This function frees a vector of doubles */
	
    if (pd!=NULL)
        free(pd);
	
} /* end vecfreed */

void vecfreei(int *pi){
    /* This function frees a vector of integers */
	
    if (pi!=NULL)
        free(pi);
	
} /* end vecfreei */

void matfreed(double **ppd){
    /* This function frees a matrix of doubles */
	
    if (ppd!=NULL){
        if (ppd[0]!=NULL)
            free(ppd[0]);
        free(ppd);
    }
	
} /* end matfreed */

void toCyclic(UINT32 *x, 
			  double *compX,
			  UINT32 lengthOfx,
			  UINT32 numProcs,
			  UINT32 procID) { 
	UINT32 i,j, g, l, offset, dst;
		
	double *temp = (double *)calloc(lengthOfx, sizeof(double));
	
	
	int packetSize = lengthOfx/numProcs;
	for(i=0; i<lengthOfx; i++)
	{
		g=i+procID*lengthOfx;
		dst=g%numProcs;
		l=(g-dst)/numProcs;
	
		offset = packetSize*dst-packetSize*procID;
		temp[offset+l]=(double)x[i];
				
	}

	MPI_Alltoall(temp, packetSize, MPI_DOUBLE, compX, packetSize, MPI_DOUBLE, MPI_COMM_WORLD );
	free(temp);

}

void toBlock(double *a,
			 double *b,
			 UINT32 p, 
			 UINT32 s,
			 UINT32 n) 
{ 	
	double *temp = (double *)calloc(n, sizeof(double));

	UINT32 g, dst, l, i, j, offset;

	int packetSize = n/p;
	
	MPI_Alltoall(a, packetSize, MPI_DOUBLE, temp, packetSize, MPI_DOUBLE, MPI_COMM_WORLD );

	
	int k, h, m;
	for (i = 0; i < n; i++) 
	{ 
		g = i * p + s; 
		dst = g / n; 
		l = g % n;
		
		k=i+s*n;
		h=k%p;
		m=(k-dst)/p;
		
		b[l+dst-s]=temp[i];	
	} 
	
	free(temp);	

}

IEEE_754_FloatNum create_IEEE_NULL(UINT32 precision, UINT32 numProcs, UINT32 procID)
{
	IEEE_754_FloatNum x;
	x.base	= BASE;
	x.exp = 0;
	x.sign	= FALSE;
	x.significand = (UINT32 *)calloc(precision/numProcs, sizeof(UINT32));
	memset( x.significand, (UINT32)0, precision/numProcs*sizeof(UINT32) );
	x.precision = precision;
	return x;
}

void kill_IEEE_NULL(IEEE_754_FloatNum x)
{
	free(x.significand);
}

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "sort.h"
#define  TRUE 1

//main(int argc,char *argv[])
//{
//int DataSize;
//int *data;
//int	MyID, SumID;
//int i, j;
//int m, r;
//
//MPI_Status status;
//MPI_Init(&argc,&argv);
//
//MPI_Comm_rank(MPI_COMM_WORLD,&MyID);	
//
//
//MPI_Comm_size(MPI_COMM_WORLD,&SumID);	
//
//if(MyID==0)
//{
//	DataSize=GetDataSize();
//	data=(int *)malloc(DataSize*sizeof(int));
//	
//	if(data==0) 
//		ErrMsg("Malloc memory error!");
//	
//	srand(396);
//	for(i=0;i&lt;DataSize;i++)
//	{
//		data[i]=(int)rand();
//		printf("%10d",data[i]);
//	}
//	printf("\n");
//}
//
//m=log2(SumID);
//
//MPI_Bcast(&DataSize,1,MPI_INT,0,MPI_COMM_WORLD);
//
//para_QuickSort(data,0,DataSize-1,m,0,MyID);
//
//if(MyID==0)
//{
//	for(i=0;i&lt;DataSize;i++)
//	{
//		printf("%10d",data[i]);
//	}
//	printf("\n");
//}
//
//MPI_Finalize();		
//}


void para_QuickSort(int *data,int start,int end,int m,int id,int MyID)
{
	int i, j;
	int r;
	int MyLength;
	int *tmp;
	MPI_Status status;
	
	MyLength=-1;
	
	if(m==0)
	{
		if(MyID==id)
			QuickSort(data,start,end);
		return;
	}	
		
	if(MyID==id)
	{
		r=Partition(data,start,end);	
		MyLength=end-r;	
	
		MPI_Send(&MyLength,1,MPI_INT,id+exp_2(m-1),MyID,MPI_COMM_WORLD);
		if(MyLength!=0)
			MPI_Send(data+r+1,MyLength,MPI_INT,id+exp_2(m-1),MyID,MPI_COMM_WORLD);
	}
	
	if(MyID==id+exp_2(m-1))
	{
		MPI_Recv(&MyLength,1,MPI_INT,id,id,MPI_COMM_WORLD,&status);
		if(MyLength!=0)
		{
			tmp=(int *)malloc(MyLength*sizeof(int));
			if(tmp==0) ErrMsg("Malloc memory error!");
			MPI_Recv(tmp,MyLength,MPI_INT,id,id,MPI_COMM_WORLD,&status);
		}
	}
	
	j=r-1-start;	
	MPI_Bcast(&j,1,MPI_INT,id,MPI_COMM_WORLD);
	if(j>0)
		para_QuickSort(data,start,r-1,m-1,id,MyID);
	
	j=MyLength;	
	MPI_Bcast(&j,1,MPI_INT,id,MPI_COMM_WORLD);
	/*(1.5)	para_quicksort(data,r+1,j,m-1,id+2m-1)*/
	if(j>0)
		para_QuickSort(tmp,0,MyLength-1,m-1,id+exp_2(m-1),MyID);
	
	/*(1.6)	P(id+2m-1) send data[r+1,m-1] back to Pid */
	if((MyID==id+exp_2(m-1)) && (MyLength!=0))
		MPI_Send(tmp,MyLength,MPI_INT,id,id+exp_2(m-1),MPI_COMM_WORLD);
	
	if((MyID==id) && (MyLength!=0))
		MPI_Recv(data+r+1,MyLength,MPI_INT,id+exp_2(m-1),id+exp_2(m-1),MPI_COMM_WORLD,&status);
	
}

void QuickSort(int *data,int start,int end)
{
	int r;
	int i;
	
	if(start<end)
	{
		r=Partition(data,start,end);
		QuickSort(data,start,r-1);
		QuickSort(data,r+1,end);
	}
}


int Partition(int *data,int start,int end)
{
	int pivo;
	int i, j;
	int tmp;
	
	pivo=data[end];			
	i=start-1;				
	
	for(j=start;j<end;j++)
		if(data[j]<=pivo)
		{
			i++;			
			tmp=data[i];
			data[i]=data[j];
			data[j]=tmp;
		}
	
	tmp=data[i+1];
	data[i+1]=data[end];
	data[end]=tmp;			
	
	return i+1;
}


int exp_2(int num)
{
	int i;
	
	i=1;
	
	while(num>0)
	{
		num--;
		i=i*2;
	}
	
	return i;
}


int log_2(int num)
{
	int i, j;
	
	i=1;	
	j=2;
	
	while(j<num)
	{
		j=j*2;
		i++;
	}
	
	if(j>num)
		i--;
	
	return i;
}


int GetDataSize(void)
{
	int i;
	
	while(TRUE)
	{
		printf("Input the Data Size :");
		scanf("%d",&i);
		if((i>0) && (i<=65535))		
			
			break;
		ErrMsg("Wrong Data Size, must between [1..65535]");
	}
	return i;
}


void ErrMsg(char *msg)
{
	printf("Error: %s \n",msg);
}
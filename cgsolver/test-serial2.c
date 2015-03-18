#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

void cgsolver(double **, double *, double *, int, double, int);
void show_vec(double *,int);
void show_mat(double **,int,int);
void mtx_vec_mul(double *,double **,double *,int);

int main(int argc, char *argv[])
{
	// Test
	double ** A, ** B;
	double *b, *x, *tmpvec;
	double size = 25;
	double criterion = 1.0e-6;
	int i,j,k;
	
    A = (double **) malloc( size*sizeof( double *) );
	B = (double **) malloc( size*sizeof( double *) );
	b = (double * ) malloc( size*sizeof( double ) );
	x = (double * ) malloc( size*sizeof( double ) );
	tmpvec = (double *) malloc( size*sizeof( double ) );
	srand(time(NULL));
	for(i=0;i<size;i++)
	{
		*(A+i) = (double *)malloc( size*sizeof( double ) );
		*(B+i) = (double *)malloc( size*sizeof( double ) );
		*(b+i) = (double)(rand()-RAND_MAX/2.0)/(double)(RAND_MAX);
		for(j=0;j<size;j++)
			*(*(A+i)+j) = (double)(rand()-RAND_MAX/2.0)/(double)(RAND_MAX);
	}
	
	//calculate B = A * A^T
	for(i=0;i<size;i++)
	{
		for(j=0;j<size;j++)
		{
			*(*(B+i)+j) = 0.0;
			for(k=0;k<size;k++)
			{
				*(*(B+i)+j) += ( *(*(A+i)+k) ) * ( *(*(A+j)+k) );
			}
		}
	}

	printf (" A= ");
	show_mat(A,size,size);

	printf (" B= ");
	show_mat(B,size,size);

	printf ("\n b=");
	show_vec(b,size);

	cgsolver(B, x, b, size, criterion, 1);

	printf ("\n B*x =");
    mtx_vec_mul(tmpvec,B,x,size);
	show_vec(tmpvec,size);

	printf ("\n solution: x=");
	show_vec(x,size);

    return 0;
}

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void cgsolver(double **, double *, double *, int, double, int);
double vec_vec_mul(double *,double *,int);
double vec_mtx_vec_mul(double *,double **,double *,int);
void mtx_vec_mul(double *,double **,double *,int);
void vec_vec_combine(double *,double,double *,double,double *,int);
void vec_vec_copy(double *,double *,int);
void show_vec(double *,int);

int main(int argc, char *argv[])
{
	// Test
	double A0[3][3]={{2,-1,0},{-1,2,-1},{0,-1,2}};
	double ** A;
	double b[3] = {1,2.0,5.0}, x[3], tmpvec[3];
	double size = 3, tmpv;
	double criterion = 1.0e-6;
	int i,j;
	
    printf("*");
    A = (double **) malloc( size*sizeof( double *) );
	for(i=0;i<size;i++)
	{
		*(A+i) = (double *)malloc( size*sizeof( double ) );
		for(j=0;j<size;j++)
			*(*(A+i)+j) = A0[i][j];
	}   
	printf ("*");
	tmpv = vec_vec_mul(b,b,size);
	printf (" %8.5f \n ", tmpv);
    
	tmpv = vec_mtx_vec_mul(b,A,b,size);
	printf (" %8.5f \n", tmpv);

	mtx_vec_mul(tmpvec, A, b, size);

	show_vec(tmpvec,size);

	vec_vec_combine(tmpvec,3.0,b,-1.0,b,size);

	show_vec(tmpvec,size);

    vec_vec_copy(tmpvec,b,size);

	show_vec(tmpvec,size);

	cgsolver(A, x, b, size, criterion, 1);

    printf ("x= [");
	for(i=0;i<size;i++)
		printf (" %10.5f ", x[i]);
	printf (" ]\n");
    return 0;
}

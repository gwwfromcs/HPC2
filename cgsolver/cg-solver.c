#include <stdio.h>
#include <stdlib.h>
#include <math.h>



// display a vector:
void show_vec(double * x,int size)
{
	int i;
	printf ("\n[");
	for(i=0;i<size;i++)
	{
		printf (" %8.3f ", *(x+i));
		if(i<size-1) printf ("\n");
 	}
	printf ("]\n");
}

// display a mattrix:
void show_mat(double ** mtx, int nrow, int ncol)
{
	int i,j;
	printf ("\n[");
	for(i=0;i<nrow;i++)
	{
		printf ("[");
		for(j=0;j<ncol;j++)
		{
			printf (" %8.3f ", *(*(mtx+i)+j));
		}
		printf ("]");
		if(i<nrow-1) printf (";\n");
	}
	printf ("]\n");
}


//vec_vec_mul compute:
// dot product of vector x1 and vector x2 with length "size"
double vec_vec_mul(double * x1, double * x2, int size)
{
	int i;
	double dotproduct=0.0;
	for(i=0;i<size;i++)
		dotproduct += (*(x1+i)) * (*(x2+i));
	return dotproduct;
}

//mtx_vec_mul compute:
// y = mtxA * x
void mtx_vec_mul(double * y, double **mtxA, double * x, int size)
{
	int i, j;
	for(i=0;i<size;i++)
	{
		*(y+i) = 0.0;
		for(j=0;j<size;j++)
			*(y+i) += (*(*(mtxA+i)+j))  *  (*(x+j));
	}
}


//vec_mtx_vec compute:
double vec_mtx_vec_mul(double * x1, double ** mtxA, double * x2, int size)
{
	int i,j;
	double dotproduct=0.0;
	for(i=0;i<size;i++)
		for(j=0;j<size;j++)
			dotproduct += (*(x1+i)) * (*(*(mtxA+i)+j)) * (*(x2+j));
	return dotproduct;
}

//vec_vec_minus computes:
// y = a * x1 + b * x2
void vec_vec_combine(double * y, double a, double * x1, double b, double * x2, int size)
{
	int i;
	for(i=0;i<size;i++)
		*(y+i) = a * (*(x1+i)) + b * (*(x2+i));
}

//vec_vec_copy:
// copy x to y
void vec_vec_copy(double * y, double * x, int size)
{
	int i;
	for(i=0;i<size;i++)
		*(y+i) = *(x+i);
}

void cgsolver(double ** mtxA, double * x, double * b, int size, double criterion, int verbose)
{
	int i, j, loop;
	double averageA, averageb;
	double *p, *r, *tmpvec;
	double res = 1000.0, res_old, alpha, beta, product1, product2;

	p = (double *)malloc( size*sizeof( double ) );
	r = (double *)malloc( size*sizeof( double ) );
	tmpvec = (double *)malloc( size*sizeof( double ) );
	//The average value of mtxA
	averageA = 0.0;
	for(i=0;i<size;i++)
	{
		averageb += *(b+i);
		for(j=0;j<size;j++)
			averageA += *(*(mtxA+i)+j);
	}
	//initial guess of solution x
	for(i=0;i<size;i++)
		//*(x+i)=averageb/averageA;
		*(x+i) = 0.0;
	//initial setup of residue r
	// r_0 = b - A * x_0
	// p_0 = r_0
	mtx_vec_mul(r, mtxA, x, size);
	vec_vec_combine(r, 1.0, b, -1.0, r, size);
	res = vec_vec_mul( r, r, size);
	vec_vec_copy(p, r, size);
	loop = 0;
	//
	printf ("Start conjugate gradient method iteration:\n");
	while ( loop < 200 )
	{
		product1 = vec_vec_mul( r, r, size);
		product2 = vec_mtx_vec_mul( p, mtxA, p, size);
		alpha = product1/product2;
		
		vec_vec_combine(x, 1.0, x, alpha, p, size);
		mtx_vec_mul(tmpvec, mtxA, p, size);
		vec_vec_combine(r, 1.0, r, -alpha, tmpvec, size);
		
		res_old = res;
		res = vec_vec_mul( r, r, size);
		if(verbose) printf("loop: %4d, residue: %15.7f \n ",loop,sqrt(res/(double)size));

		if(sqrt(res)/(double)size < criterion) break;

		beta = res/res_old;

		vec_vec_combine(p, 1.0, r, beta, p, size);

		loop += 1;
	}
}


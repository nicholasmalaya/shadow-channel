/**********************************************************************
Author: Vanessa Lopez
	Department of Computer Science
	University of Illinois at Urbana-Champaign
	1304 W. Springfield Ave.
	Urbana, IL 61801-2987
	Email: vlopez@cse.uiuc.edu

Copyright 1999.  This code represents preliminary work toward the
author's thesis.  The code is not to be redistributed without the
author's permission.  Any corrections should be communicated to
the author.  Any modifications or reuse of the code must retain
acknowledgement of the original source.
**********************************************************************/

#include<string.h>
#include "minChnl.h"
#include "mcomplex.h"

/* Do Gaussian elimination to solve Ay=B, where A contains xdim banded matrices
   with dl+du+1 nonzero diagonals.  There is a zero diagonal between each
   nonzero diagonal of the matrices A.  Each matrix A has real entries and N 
   rows.  b is an N-by-xdim complex array with right hand side vectors to solve
   for.  dl is the number of nonzero super diagonals; du is the number of
   nonzero sub-diagonals.  x0 indicates from which matrix A we should start
   solving for.  A is overwritten with the LU factorization of the matrices;
   the result is returned in b.
*/
void bsolve(double ***A, mcomplex **b, int dl, int du, int N, int xdim, int x0)
{
    int i, j, k, x;

    /* LU factorization */
    for (i = 2; i < N; ++i)
    {
	/* coefficients in L */
	for (j = MAX(0, dl-(i/2)); j < dl; ++j)
	{
            for (k = 1; k <= j-MAX(0, dl-(i/2)); ++k)
	    {
		for (x = x0; x < xdim; ++x)
		{
		    A[i][j][x] = A[i][j][x] - 
		       A[i][j-k][x]*A[i-2*(dl-j+k)][dl+k][x]; 
		}
	    }
	    for (x = x0; x < xdim; ++x)
	    {
	        A[i][j][x] = A[i][j][x] / A[i-2*(dl-j)][dl][x]; 

	        /* Solve Lz=b (overwritting vector b). */  
	        Re(b[i][x]) = Re(b[i][x]) - A[i][j][x]*Re(b[i-2*(dl-j)][x]);  
	        Im(b[i][x]) = Im(b[i][x]) - A[i][j][x]*Im(b[i-2*(dl-j)][x]);  
	    }
	}

	/* coefficients in U */
	for (j = 0; j < MIN(du, (N-i+1)/2); ++j)
	{
            for (k = 1; k <= MIN(i/2, MIN(dl, du-j)); ++k)
	    {
	        for (x = x0; x < xdim; ++x)
		{
		    A[i][j+dl][x] = A[i][j+dl][x] - 
		       A[i][dl-k][x]*A[i-2*k][dl+j+k][x];
		}
	    }
	}
    }   /* end for i ... */

    /* Solve Ux=z */
    for (x = x0; x < xdim; ++x)
    {
        Re(b[N-1][x]) = Re(b[N-1][x]) / A[N-1][dl][x];  
        Im(b[N-1][x]) = Im(b[N-1][x]) / A[N-1][dl][x];  
    }
    for (x = x0; x < xdim; ++x)
    {
        Re(b[N-2][x]) = Re(b[N-2][x]) / A[N-2][dl][x];  
        Im(b[N-2][x]) = Im(b[N-2][x]) / A[N-2][dl][x];  
    }
    for (i = N-3; i >= 0; --i)
    {
	for (j = 1; j <= MIN(du, ((N-i+1)/2-1)); ++j)
	{
            for (x = x0; x < xdim; ++x)
            {
                Re(b[i][x]) = Re(b[i][x]) - A[i][dl+j][x]*Re(b[i+2*j][x]);  
                Im(b[i][x]) = Im(b[i][x]) - A[i][dl+j][x]*Im(b[i+2*j][x]);  
	    }
	}
        for (x = x0; x < xdim; ++x)
        {
            Re(b[i][x]) = Re(b[i][x]) / A[i][dl][x];  
            Im(b[i][x]) = Im(b[i][x]) / A[i][dl][x];  
	}

    }   /* end for i... */
}


/* Same as for  bsolve, but here A is a single matrix and the rhs vector is
   stored in the first column of b. 
*/
void bsolve0(double **A, mcomplex **b, int dl, int du, int N)
{
    int i, j, k;

    /* LU factorization */
    for (i = 2; i < N; ++i)
    {
	/* coefficients in L */
	for (j = MAX(0, dl-(i/2)); j < dl; ++j)
	{
            for (k = 1; k <= j-MAX(0, dl-(i/2)); ++k)
	    {
		A[i][j] = A[i][j] - A[i][j-k]*A[i-2*(dl-j+k)][dl+k]; 
	    }
	    A[i][j] = A[i][j] / A[i-2*(dl-j)][dl]; 

	    /* Solve Lz=b (overwritting vector b). */
	    Re(b[i][0]) = Re(b[i][0]) - A[i][j]*Re(b[i-2*(dl-j)][0]);  
	    Im(b[i][0]) = Im(b[i][0]) - A[i][j]*Im(b[i-2*(dl-j)][0]);  
	}

	/* coefficients in U */
	for (j = 0; j < MIN(du, (N-i+1)/2); ++j)
	{
            for (k = 1; k <= MIN(i/2, MIN(dl, du-j)); ++k)
	    {
	        A[i][j+dl] = A[i][j+dl] - A[i][dl-k]*A[i-2*k][dl+j+k];
	    }
	}
    }   /* end for i ... */

    /* Solve Ux=z (again, overwrite vector b) */
    Re(b[N-1][0]) = Re(b[N-1][0]) / A[N-1][dl];  
    Im(b[N-1][0]) = Im(b[N-1][0]) / A[N-1][dl];  
    
    Re(b[N-2][0]) = Re(b[N-2][0]) / A[N-2][dl];  
    Im(b[N-2][0]) = Im(b[N-2][0]) / A[N-2][dl];  
    
    for (i = N-3; i >= 0; --i)
    {
	for (j = 1; j <= MIN(du, ((N-i+1)/2-1)); ++j)
	{
            Re(b[i][0]) = Re(b[i][0]) - A[i][dl+j]*Re(b[i+2*j][0]);  
            Im(b[i][0]) = Im(b[i][0]) - A[i][dl+j]*Im(b[i+2*j][0]);  
	}
        Re(b[i][0]) = Re(b[i][0]) / A[i][dl];  
        Im(b[i][0]) = Im(b[i][0]) / A[i][dl];  

    }   /* end for i... */
}
/* Same as for  bsolve, but here A is a single matrix and the rhs vector is
   stored in the first column of b. 
*/
void bsolve_0(double **A, mcomplex *b, int dl, int du, int N)
{
    int i, j, k;

    /* LU factorization */
    for (i = 2; i < N; ++i)
    {
	/* coefficients in L */
	for (j = MAX(0, dl-(i/2)); j < dl; ++j)
	{
            for (k = 1; k <= j-MAX(0, dl-(i/2)); ++k)
	    {
		A[i][j] = A[i][j] - A[i][j-k]*A[i-2*(dl-j+k)][dl+k]; 
	    }
	    A[i][j] = A[i][j] / A[i-2*(dl-j)][dl]; 

	    /* Solve Lz=b (overwritting vector b). */
	    Re(b[i]) = Re(b[i]) - A[i][j]*Re(b[i-2*(dl-j)]);  
	    Im(b[i]) = Im(b[i]) - A[i][j]*Im(b[i-2*(dl-j)]);  
	}

	/* coefficients in U */
	for (j = 0; j < MIN(du, (N-i+1)/2); ++j)
	{
            for (k = 1; k <= MIN(i/2, MIN(dl, du-j)); ++k)
	    {
	        A[i][j+dl] = A[i][j+dl] - A[i][dl-k]*A[i-2*k][dl+j+k];
	    }
	}
    }   /* end for i ... */

    /* Solve Ux=z (again, overwrite vector b) */
    Re(b[N-1]) = Re(b[N-1]) / A[N-1][dl];  
    Im(b[N-1]) = Im(b[N-1]) / A[N-1][dl];  
    
    Re(b[N-2]) = Re(b[N-2]) / A[N-2][dl];  
    Im(b[N-2]) = Im(b[N-2]) / A[N-2][dl];  
    
    for (i = N-3; i >= 0; --i)
    {
	for (j = 1; j <= MIN(du, ((N-i+1)/2-1)); ++j)
	{
            Re(b[i]) = Re(b[i]) - A[i][dl+j]*Re(b[i+2*j]);  
            Im(b[i]) = Im(b[i]) - A[i][dl+j]*Im(b[i+2*j]);  
	}
        Re(b[i]) = Re(b[i]) / A[i][dl];  
        Im(b[i]) = Im(b[i]) / A[i][dl];  

    }   /* end for i... */
}


/* Multiply Ay, where A and y are as in the comment for bsolve. */
void smMult(double ***A, mcomplex **y, mcomplex **b, int dl, int du, int N,
    int xdim, int x0)
{
    int i, j, x;

    for (i = 0; i < N; ++i)
    {
        for (x = x0; x < xdim; ++x) 
	{
	    Re(b[i][x]) = 0.0;
	    Im(b[i][x]) = 0.0;
	}
    }
    for (i = 0; i < N; ++i)
    {
        for (j = MAX(0, dl-(i/2)); j < MIN(dl+du+1, (N-i+1)/2+dl); ++j)
	{
            for (x = x0; x < xdim; ++x) 
	    {
	        Re(b[i][x]) = Re(b[i][x]) + A[i][j][x]*Re(y[i-2*(dl-j)][x]);
	        Im(b[i][x]) = Im(b[i][x]) + A[i][j][x]*Im(y[i-2*(dl-j)][x]);
	    }
	}
    }
}

/* Multiply Ay, where A and y are as in the comment for bsolve0. */
void smMult0(double **A, mcomplex **y, mcomplex *b, int dl, int du, int N)
{
    int i, j;

    for (i = 0; i < N; ++i)
    {
	Re(b[i]) = 0.0;
	Im(b[i]) = 0.0;
        for (j = MAX(0, dl-(i/2)); j < MIN(dl+du+1, (N-i+1)/2+dl); ++j)
	{
	    Re(b[i]) = Re(b[i]) + A[i][j]*Re(y[i-2*(dl-j)][0]);
	    Im(b[i]) = Im(b[i]) + A[i][j]*Im(y[i-2*(dl-j)][0]);
	}
    }
}

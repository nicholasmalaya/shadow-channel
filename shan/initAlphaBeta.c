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

#include <stdio.h>
#include <string.h>
#include "minChnl.h"
#include "mvOps.h"
#include "arrays.h"

/* Compute initial coefficients in expansion, given the values for ux_hat,
   uy_hat, uz_hat */
void initAlphaBeta(void)
{
    /* function prototypes */
    void Alpha0(void);
    void Beta0(void);
    void Alpha(int z, int x0, int xdim);
    void Beta(int z, int x0, int xdim);

    /* variables */
    extern int Nz, Nx;
    int z;

    /* (Kx,Kz) = (0,0) */
    Alpha0();
    Beta0();

    /* Kz = 0; Kx != 0 */
    Alpha(0, 1, Nx/2);
    Beta(0, 1, Nx/2);

    /* Kz != 0; Kx = 0,... */
    for (z = 1; z < Nz; ++z)
    {
        Alpha(z, 0, Nx/2);
        Beta(z, 0, Nx/2);
    }


}    /* end initAlphaBeta */


void Alpha0(void)
{
    extern int qpts, dimR;
    extern double **Rw, **Rs, **MZ, *Rp0;
    extern mcomplex ****U, ****C;
    extern mcomplex **Uxb, **Uzb;
    int i, j;

    for (i = 0; i < dimR; ++i)
    {
        Re(C[0][ALPHA][i][0]) = 0.0;
        Im(C[0][ALPHA][i][0]) = 0.0;
        for (j = 0; j < qpts; ++j)
	{
            Re(C[0][ALPHA][i][0]) += Rw[i][j]*Re(U[0][XEL][j][0]);
            Im(C[0][ALPHA][i][0]) += Rw[i][j]*Im(U[0][XEL][j][0]);
	}
    }

    memcpy(MZ[0], Rs[0], dimR*T_RSDIAG*sizeof(double));
    bsolve0(MZ, C[0][ALPHA], RSDIAG-1, RSDIAG-1, dimR);

 /* Uxb=Rp0*a, Uzb=Rp0*b, used as boundary condition for the incremental system*/

    for (j=0; j< dimR; ++j)
      {
	Re(Uxb[0][0]) += Rp0[j]*Re(C[0][ALPHA][j][0]);
	Im(Uxb[0][0]) += Rp0[j]*Im(C[0][ALPHA][j][0]);
      }


}    /* end Alpha0() */


void Beta0(void)
{
    extern int qpts, dimR;
    extern double **Rw, **Rs, **MZ,*Rp0;
    extern mcomplex ****U, ****C;
    extern mcomplex **Uxb, **Uzb;
    int i, j;

    for (i = 0; i < dimR; ++i)
    {
        Re(C[0][BETA][i][0]) = 0.0;
        Im(C[0][BETA][i][0]) = 0.0;
        for (j = 0; j < qpts; ++j)
	{
            Re(C[0][BETA][i][0]) += Rw[i][j]*Re(U[0][ZEL][j][0]);
            Im(C[0][BETA][i][0]) += Rw[i][j]*Im(U[0][ZEL][j][0]);
	}
    }

    memcpy(MZ[0], Rs[0], dimR*T_RSDIAG*sizeof(double));
    bsolve0(MZ, C[0][BETA], RSDIAG-1, RSDIAG-1, dimR);
 
 /* Uxb=Rp0*a, Uzb=Rp0*b, used as boundary condition for the incremental system*/

    for (j=0; j< dimR; ++j)
      {
	Re(Uzb[0][0]) += Rp0[j]*Re(C[0][BETA][j][0]);
	Im(Uzb[0][0]) += Rp0[j]*Im(C[0][BETA][j][0]);
      }



}    /* end Beta0() */


void Alpha(int z, int x0, int xdim)
{
    void bsolveI(double **A, mcomplex **b, int dl, int du, int N, int xdim,
       int x0);

    extern int qpts, dimQ;
    extern double **Qw, **Qs;
    extern mcomplex ****U, ****C;

    int i, j, x;
    double **M;

    if ( (M = dMatrix(dimQ, T_QSDIAG)) == NULL)
    {
	printf("Alpha: No memory for M array.\n");
	return;
    }

    for (i = 0; i < dimQ; ++i)
    {
        for (x = x0; x < xdim; ++x)
	{
            Re(C[z][ALPHA][i][x]) = 0.0;
            Im(C[z][ALPHA][i][x]) = 0.0;
	}
    }
    for (i = 0; i < dimQ; ++i)
    {
        for (j = 0; j < qpts; ++j)
	{
            for (x = x0; x < xdim; ++x)
	    {
                Re(C[z][ALPHA][i][x]) += Qw[i][j]*Re(U[z][YEL][j][x]);
                Im(C[z][ALPHA][i][x]) += Qw[i][j]*Im(U[z][YEL][j][x]);
	    }
	}
    }

    memcpy(M[0], Qs[0], dimQ*T_QSDIAG*sizeof(double));
    bsolveI(M, C[z][ALPHA], QSDIAG-1, QSDIAG-1, dimQ, xdim, x0);
    freedMatrix(M);
}    /* end Alpha() */


void Beta(int z, int x0, int xdim)
{
    void bsolveI(double **A, mcomplex **b, int dl, int du, int N, int xdim,
       int x0);

    extern int qpts, dimR;
    extern double *Kx, *Kz, **K2;
    extern double **Rw, **Rs, **MZ, *Rp0;
    extern mcomplex ****U, ****C;
    extern mcomplex **Uxb, **Uzb;
    int i, j, x;

    for (i = 0; i < dimR; ++i)
    {
        for (x = x0; x < xdim; ++x)
	{
            Re(C[z][BETA][i][x]) = 0.0;
            Im(C[z][BETA][i][x]) = 0.0;
	}
    }
    for (i = 0; i < dimR; ++i)
    {
        for (j = 0; j < qpts; ++j)
	{
            for (x = x0; x < xdim; ++x)
	    {
                Re(C[z][BETA][i][x]) += Rw[i][j] * 
		   (Kx[x]*Im(U[z][ZEL][j][x]) - Kz[z]*Im(U[z][XEL][j][x]));
                Im(C[z][BETA][i][x]) += Rw[i][j] * 
		   (Kz[z]*Re(U[z][XEL][j][x]) - Kx[x]*Re(U[z][ZEL][j][x]));
	    }
	}
    }

    memcpy(MZ[0], Rs[0], dimR*T_RSDIAG*sizeof(double));
    bsolveI(MZ, C[z][BETA], RSDIAG-1, RSDIAG-1, dimR, xdim, x0);

 /* compute dux duz at y=-1 as the boundary condition for the incremental systems*/

    for (x=x0; x< xdim; ++x)
      {
	for (j=0; j < dimR; ++j)
	  {
	    Re(Uxb[z][x]) += Rp0[j]*Im(C[z][BETA][j][x])*Kz[z]/K2[z][x];
	    Im(Uxb[z][x]) +=-Rp0[j]*Re(C[z][BETA][j][x])*Kz[z]/K2[z][x];

	    Re(Uzb[z][x]) +=-Rp0[j]*Im(C[z][BETA][j][x])*Kx[x]/K2[z][x];
	    Im(Uzb[z][x]) += Rp0[j]*Re(C[z][BETA][j][x])*Kx[x]/K2[z][x];
	  }
      }
 
}    /* end Beta() */


void bsolveI(double **A, mcomplex **b, int dl, int du, int N, int xdim, int x0)
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
	        A[i][j] = A[i][j] - A[i][j-k]*A[i-2*(dl-j+k)][dl+k]; 
	    }
	    A[i][j] = A[i][j] / A[i-2*(dl-j)][dl]; 

	    /* Solve Lz=b (overwritting vector b). */  
	    for (x = x0; x < xdim; ++x)
	    {
	        Re(b[i][x]) = Re(b[i][x]) - A[i][j]*Re(b[i-2*(dl-j)][x]);  
	        Im(b[i][x]) = Im(b[i][x]) - A[i][j]*Im(b[i-2*(dl-j)][x]);  
	    }
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
    for (x = x0; x < xdim; ++x)
    {
        Re(b[N-1][x]) = Re(b[N-1][x]) / A[N-1][dl];  
        Im(b[N-1][x]) = Im(b[N-1][x]) / A[N-1][dl];  
    }
    for (x = x0; x < xdim; ++x)
    {
        Re(b[N-2][x]) = Re(b[N-2][x]) / A[N-2][dl];  
        Im(b[N-2][x]) = Im(b[N-2][x]) / A[N-2][dl];  
    }
    for (i = N-3; i >= 0; --i)
    {
	for (j = 1; j <= MIN(du, ((N-i+1)/2-1)); ++j)
	{
            for (x = x0; x < xdim; ++x)
            {
                Re(b[i][x]) = Re(b[i][x]) - A[i][dl+j]*Re(b[i+2*j][x]);  
                Im(b[i][x]) = Im(b[i][x]) - A[i][dl+j]*Im(b[i+2*j][x]);  
	    }
	}
        for (x = x0; x < xdim; ++x)
        {
            Re(b[i][x]) = Re(b[i][x]) / A[i][dl];  
            Im(b[i][x]) = Im(b[i][x]) / A[i][dl];  
	}

    }   /* end for i... */
}       /* end bsolveI */

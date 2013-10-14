/*******************************************************************
This function reads in the stored Fourier coefficients a, b, alpha, beta
and compute the velocity, ux, uy, uz, dux, duz.
******************************************************************/

#include <stdio.h>
#include <string.h>
#include "minChnl.h"
#include "mvOps.h"
#include "arrays.h"

void initAlphaBeta2(void)
{

  void initAlphaBeta_2(int z, int x0);

  extern int qpts, dimR, dimQ, Nx, Nz;
  extern mcomplex ****U, ****C;
  extern mcomplex **Uxb, **Uzb; 
  extern mcomplex **HUxb, **HUzb;
  extern double **Q, **Qp, **Qpp, **R, **Rp, *Rp0, *Rpp0;
extern double *Kx, *Kz, **K2;
 int i, j, x, z;
 double t[2];

   /* NOW COMPUTE U HATS */
    /* ux_hat = R*a */
    for (i = 0; i < qpts; ++i)
    {
	Re(U[0][XEL][i][0]) = 0.0;
	Im(U[0][XEL][i][0]) = 0.0;
        for (j = 0; j < dimR; ++j)
	{
	    Re(U[0][XEL][i][0]) += R[i][j]*Re(C[0][ALPHA][j][0]);
	    Im(U[0][XEL][i][0]) += R[i][j]*Im(C[0][ALPHA][j][0]);
	}
    }

    /* dux_hat = Rp*a */
    for (i = 0; i < qpts; ++i)
    {
	Re(U[0][DXEL][i][0]) = 0.0;
	Im(U[0][DXEL][i][0]) = 0.0;
        for (j = 0; j < dimR; ++j)
	{
	    Re(U[0][DXEL][i][0]) += Rp[i][j]*Re(C[0][ALPHA][j][0]);
	    Im(U[0][DXEL][i][0]) += Rp[i][j]*Im(C[0][ALPHA][j][0]);
	}
    }
    /* uy_hat = 0 */
    for (i = 0; i < qpts; ++i)
    {
	Re(U[0][YEL][i][0]) = 0.0;
	Im(U[0][YEL][i][0]) = 0.0;
    }

    /* uz_hat = R*b */
    for (i = 0; i < qpts; ++i)
    {
	Re(U[0][ZEL][i][0]) = 0.0;
	Im(U[0][ZEL][i][0]) = 0.0;
        for (j = 0; j < dimR; ++j)
	{
	    Re(U[0][ZEL][i][0]) += R[i][j]*Re(C[0][BETA][j][0]);
	    Im(U[0][ZEL][i][0]) += R[i][j]*Im(C[0][BETA][j][0]);
	}
    }

    /* duz_hat = Rp*b */
    for (i = 0; i < qpts; ++i)
    {
	Re(U[0][DZEL][i][0]) = 0.0;
	Im(U[0][DZEL][i][0]) = 0.0;
        for (j = 0; j < dimR; ++j)
	{
	    Re(U[0][DZEL][i][0]) += Rp[i][j]*Re(C[0][BETA][j][0]);
	    Im(U[0][DZEL][i][0]) += Rp[i][j]*Im(C[0][BETA][j][0]);
	}
    }

    /* Uxb=Rp0*a, Uzb=Rp0*b, used as boundary condition for the incremental system*/
    memset(Uxb[0], 0, Nz*(Nx/2)*sizeof(fftw_complex));
    memset(Uzb[0], 0, Nz*(Nx/2)*sizeof(fftw_complex));
    memset(HUxb[0], 0, Nz*(Nx/2)*sizeof(mcomplex));
    memset(HUzb[0], 0, Nz*(Nx/2)*sizeof(mcomplex));
    for (j=0; j< dimR; ++j)
      {
	Re(Uxb[0][0]) += Rp0[j]*Re(C[0][ALPHA][j][0]);
	Im(Uxb[0][0]) += Rp0[j]*Im(C[0][ALPHA][j][0]);

	Re(Uzb[0][0]) += Rp0[j]*Re(C[0][BETA][j][0]);
	Im(Uzb[0][0]) += Rp0[j]*Im(C[0][BETA][j][0]);
      }
    for (j=0; j< dimR; ++j)
      {
	Re(HUxb[0][0]) += Rpp0[j]*Re(C[0][ALPHA][j][0]);
	Im(HUxb[0][0]) += Rpp0[j]*Im(C[0][ALPHA][j][0]);

	Re(HUzb[0][0]) += Rpp0[j]*Re(C[0][BETA][j][0]);
	Im(HUzb[0][0]) += Rpp0[j]*Im(C[0][BETA][j][0]);
      }
    initAlphaBeta_2(0, 1);
    for (z=1; z< Nz; ++z)
      {
	initAlphaBeta_2(z, 0);
      }

}
void initAlphaBeta_2(int z, int x0)
{
  extern int qpts, dimR, dimQ, Nx, Nz;
  extern mcomplex ****U, ****C;
  extern mcomplex **Uxb, **Uzb;
  extern mcomplex **HUxb, **HUzb;
  extern double **Q, **Qp, **Qpp, **R, **Rp, *Rp0, *Rpp0;
extern double *Kx, *Kz, **K2;
 int i, j, x;
 double t[2];

  /* NOW COMPUTE U HATS */
    /* v = uy_hat. */
    for (i = 0; i < qpts; ++i)
    {
        memset(&U[z][YEL][i][x0], 0, (Nx/2-x0)*sizeof(mcomplex));
    }
    for (i = 0; i < qpts; ++i)
    {
        for (j = 0; j < dimQ; ++j)
	{
            for (x = x0; x < Nx/2; ++x)
	    {
	        Re(U[z][YEL][i][x]) += Q[i][j]*Re(C[z][ALPHA][j][x]);
	        Im(U[z][YEL][i][x]) += Q[i][j]*Im(C[z][ALPHA][j][x]);
	    }
	}
    }

    /* f = -dv/dy and store temporarily in XEL position of array U. */
    for (i = 0; i < qpts; ++i)
    {
        memset(&U[z][XEL][i][x0], 0, (Nx/2-x0)*sizeof(mcomplex));
    }
    for (i = 0; i < qpts; ++i)
    {
        for (j = 0; j < dimQ; ++j)
	{
            for (x = x0; x < Nx/2; ++x)
	    {
	        Re(U[z][XEL][i][x]) -= Qp[i][j]*Re(C[z][ALPHA][j][x]);
	        Im(U[z][XEL][i][x]) -= Qp[i][j]*Im(C[z][ALPHA][j][x]);
	    }
	}
    }

    /* sum(Q''alpha) and store in DXEL position. */
    for (i = 0; i < qpts; ++i)
    {
        memset(&U[z][DXEL][i][x0], 0, (Nx/2-x0)*sizeof(mcomplex));
    }
    for (i = 0; i < qpts; ++i)
    {
        for (j = 0; j < dimQ; ++j)
	{
            for (x = x0; x < Nx/2; ++x)
	    {
	        Re(U[z][DXEL][i][x]) += Qpp[i][j]*Re(C[z][ALPHA][j][x]);
	        Im(U[z][DXEL][i][x]) += Qpp[i][j]*Im(C[z][ALPHA][j][x]);
	    }
	}
    }

    /* Compute g = sum(beta*R) and store temporarily in ZEL position of 
       array U. */
    for (i = 0; i < qpts; ++i)
    {
        memset(&U[z][ZEL][i][x0], 0, (Nx/2-x0)*sizeof(mcomplex));
    }
    for (i = 0; i < qpts; ++i)
    {
        for (j = 0; j < dimR; ++j)
	{
            for (x = x0; x < Nx/2; ++x)
	    {
	        Re(U[z][ZEL][i][x]) += R[i][j]*Re(C[z][BETA][j][x]);
	        Im(U[z][ZEL][i][x]) += R[i][j]*Im(C[z][BETA][j][x]);
	    }
	}
    }

    /* Compute sum(beta*R') and store temporarily in DZEL position of
       array U. */ 
    for (i = 0; i < qpts; ++i)
    {
        memset(&U[z][DZEL][i][x0], 0, (Nx/2-x0)*sizeof(mcomplex));
    }
    for (i = 0; i < qpts; ++i)
    {
        for (j = 0; j < dimR; ++j)
	{
            for (x = x0; x < Nx/2; ++x)
	    {
	        Re(U[z][DZEL][i][x]) += Rp[i][j]*Re(C[z][BETA][j][x]);
	        Im(U[z][DZEL][i][x]) += Rp[i][j]*Im(C[z][BETA][j][x]);
	    }
	}
    }

    /* now compute ux hat, uz hat */
    for (i = 0; i < qpts; ++i)
    {
	for (x = x0; x < Nx/2; ++x)
	{
	   t[0] = Re(U[z][XEL][i][x]);	/* real part of f */
	   t[1] = Re(U[z][ZEL][i][x]);	/* real part of g */

	   Re(U[z][XEL][i][x]) = (Kx[x]*Im(U[z][XEL][i][x]) +
	      Kz[z]*Im(U[z][ZEL][i][x])) / K2[z][x];
	   Re(U[z][ZEL][i][x]) = (-Kx[x]*Im(U[z][ZEL][i][x]) +
	      Kz[z]*Im(U[z][XEL][i][x])) / K2[z][x];
	   Im(U[z][XEL][i][x]) = -(Kx[x]*t[0] + Kz[z]*t[1]) / K2[z][x];
	   Im(U[z][ZEL][i][x]) = (Kx[x]*t[1] - Kz[z]*t[0]) / K2[z][x];
	}
    }

    /* dux hat, duz hat */
    for (i = 0; i < qpts; ++i)
    {
	for (x = x0; x < Nx/2; ++x)
	{
	   t[0] = Re(U[z][DXEL][i][x]);	/* real part of Q''alpha */
	   t[1] = Re(U[z][DZEL][i][x]);	/* real part of R'beta */

	   Re(U[z][DXEL][i][x]) = (-Kx[x]*Im(U[z][DXEL][i][x]) +
	      Kz[z]*Im(U[z][DZEL][i][x])) / K2[z][x];
	   Re(U[z][DZEL][i][x]) = -(Kx[x]*Im(U[z][DZEL][i][x]) +
	      Kz[z]*Im(U[z][DXEL][i][x])) / K2[z][x];
	   Im(U[z][DXEL][i][x]) = (Kx[x]*t[0] - Kz[z]*t[1]) / K2[z][x];
	   Im(U[z][DZEL][i][x]) = (Kx[x]*t[1] + Kz[z]*t[0]) / K2[z][x];
	}
    }
  
    /* compute dux duz at y=-1 as the boundary condition for the incremental systems*/

    for (x=x0; x< Nx/2; ++x)
      {
	for (j=0; j < dimR; ++j)
	  {
	    Re(Uxb[z][x]) += Rp0[j]*Im(C[z][BETA][j][x])*Kz[z]/K2[z][x];
	    Im(Uxb[z][x]) +=-Rp0[j]*Re(C[z][BETA][j][x])*Kz[z]/K2[z][x];

	    Re(Uzb[z][x]) +=-Rp0[j]*Im(C[z][BETA][j][x])*Kx[x]/K2[z][x];
	    Im(Uzb[z][x]) += Rp0[j]*Re(C[z][BETA][j][x])*Kx[x]/K2[z][x];
	
	  }
      }
  for (x=x0; x< Nx/2; ++x)
      {
	for (j=0; j < dimR; ++j)
	  {
	    Re(HUxb[z][x]) += Rpp0[j]*Im(C[z][BETA][j][x])*Kz[z]/K2[z][x];
	    Im(HUxb[z][x]) +=-Rpp0[j]*Re(C[z][BETA][j][x])*Kz[z]/K2[z][x];

	    Re(HUzb[z][x]) +=-Rpp0[j]*Im(C[z][BETA][j][x])*Kx[x]/K2[z][x];
	    Im(HUzb[z][x]) += Rpp0[j]*Re(C[z][BETA][j][x])*Kx[x]/K2[z][x];

	  }
      }
}

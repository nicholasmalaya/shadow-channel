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
/***********************************************************
   Update on Oct. 2: add a new function 
   
       int CreateRp0(double *Rp0, int dimR, double y)
 which creates vector Rp0, with the first derivative of R evaluated at y=-1 and it
 will be used when computing the boundary conditions for the incremental state equations.
 
Note that there is no corresponding vector (Qp0) because the first order derivative
 of Q evaluated at y=-1 is always zero.

************************************************************/

#include <stdlib.h>
#include "minChnl.h"

/* Compute R and Q matrices and their derivatives.
   R(y) = (1-y^2)*P(y);
   Q(y) = (1-y^2)^2*P(y), where P represents the Legendre polynomials.

   In each case, the parameters dim1 and dim2 denote the number of rows and
   colums, respectively, of the matrices.  wrk is a double array of size dim1
   used for temporary computations, and y is a double array of size dim1 that
   contains the quadrature points (or values used in evaluating the functions).
   P is a dim1-by-dim2 matrix with P(m,n) = P_n(y(m)), P_n the nth Legendre
   polynomial.

   FYI... Recurrence relation for Legendre polynomials and derivatives
 	P_0(y)    = 1
 	P_1(y)    = y
 	P_n+1(y)  = (2n+1)/(n+1) yP_n  -  n/(n+1) P_n-1,  n=1,2,...

	P_n'(y)  = 1/(1-y^2) * (nP_n-1 - nyP_n)
	P_n''(y) = (2ny)/(1-y^2)^2 * P_n-1 
	    - (2ny^2 + n(n+1)(1-y^2))/(1-y^2)^2 * P_n
*/


/* Create matrix R, whose nth column is given by 
   (1 - y^2) * nth column of P (pointwise multiplication)
*/
void CreateR(double **R, double **P, int dim1, int dim2, double *y, double *wrk)
{
    int i, j;

    for (i=0; i<dim1; ++i)
    {
	wrk[i] = 1.0 - y[i]*y[i];
    }

    for (i=0; i<dim1; ++i)
    {
        for (j=0; j<dim2; ++j)
	{
	    R[i][j] = wrk[i]*P[i][j]; 
	}
    }
}


/* Create matrix Rp, with the 1st derivative values of R.
   Rp_m(y) = m*P_m-1 - (2+m)y*P_m, for m = 1,2,...
*/
void CreateRp(double **Rp, double **P, int dim1, int dim2, double *y)
{
    int i, j;

    for (i=0; i<dim1; ++i)
    {
	Rp[i][0] = -2.0*y[i];
    }

    for (i=0; i<dim1; ++i)
    {
        for (j=1; j<dim2; ++j)
	{
	    Rp[i][j] = j*P[i][j-1] - (j+2.0)*(y[i]*P[i][j]);
	}
    }
}


/* Create matrix Q, whose nth column is given by 
   (1 - y^2)^2 * nth column of P (pointwise multiplication), which is 
   (1 - y^2) * nth column of R.
*/
void CreateQ(double **Q, double **R, int dim1, int dim2, double *y,
    double *wrk)
{
    int i, j;

    for (i=0; i<dim1; ++i)
    {
	wrk[i] = 1.0 - y[i]*y[i];
    }

    for (i=0; i<dim1; ++i)
    {
        for (j=0; j<dim2; ++j)
	{
	    Q[i][j] = wrk[i]*R[i][j]; 
	}
    }
}


/* Create matrix Qp, with the 1st derivative values of Q.
   Qp_m(y) = m(1-y^2)*P_m-1 - (m+4)y(1-y^2)*P_m, for m = 1,2,...
*/
void CreateQp(double **Qp, double **P, int dim1, int dim2, 
   double *y, double *wrk)
{
    int i, j;

    for (i=0; i<dim1; ++i)
    {
	wrk[i] = 1.0 - y[i]*y[i];
	Qp[i][0] = -4.0*y[i]*wrk[i];
    }

    for (i=0; i<dim1; ++i)
    {
        for (j=1; j<dim2; ++j)
	{
	    Qp[i][j] = wrk[i] * (j*P[i][j-1] - (j+4.0)*(y[i]*P[i][j]));
	}
    }
}


/* Create matrix Qpp, with the 2nd derivative values of Q.
   Qpp_m(y) = -6m*y*P_m-1 - (2m*y^2 + m(m+1)(1-y^2))*P_m + 8m*y^2*P_m 
       + (12y^2 - 4)*P_m, m=1,2,...
*/
void CreateQpp(double **Qpp, double **P, int dim1, int dim2, 
   double *y, double *wrk)
{
    int i, j;

    for (i=0; i<dim1; ++i)
    {
	wrk[i] = y[i]*y[i];
	Qpp[i][0] = 12.0*wrk[i] - 4.0;
    }

    for (i=0; i<dim1; ++i)
    {
        for (j=1; j<dim2; ++j)
	{
            Qpp[i][j] = (12.0 + 7.0*j + j*j)*(wrk[i]*P[i][j]) -
		(4.0 + j*j + j)*P[i][j] - 6.0*j*(y[i]*P[i][j-1]);
	}
    }
}


/* create vector Rp0, with the first derivative of R evaluated at y=-1*/

int CreateRp0(double *Rp0, int dimR, double y)
{
  int n;
  double *P ;

  if ( (P = (double *)calloc(dimR, sizeof(double))) == NULL ) 
    {
	return(NO_MEM);
    }

  /* evaluate P_n(y) at given y*/

   P[0]=1.0;
   P[1]=y;
     
    for (n=1; n<(dimR-1); ++n)
     {
            P[n+1]=((2.*n+1.)*(y*P[n]) - n*P[n-1]) / (double)(n+1);
      }

 /* Rp_m(y) = m*P_m-1 - (2+m)y*P_m, for m = 1,2,..*/

     Rp0[0]=-2.0*y;

     for (n=1; n<dimR; ++n)
     {
	    Rp0[n] = n*P[n-1] - (n+2.0)*(y*P[n]);
     }

     free(P);
   return(NO_ERR);
}

int CreateRpp0(double *Rpp0, int dimR, double y)
{
  int n, sign;
  double *P ;

  if ( (P = (double *)calloc(dimR, sizeof(double))) == NULL ) 
    {
	return(NO_MEM);
    }

  /* evaluate P_n(y) at given y*/

   P[0]=1.0;
   P[1]=y;

    for (n=1; n<(dimR-1); ++n)
     {
            P[n+1]=((2.*n+1.)*(y*P[n]) - n*P[n-1]) / (double)(n+1);
      }

 /* Rpp_m(y) = m*P'_m-1 - (2+m)y*P'_m-(m+2)*P_m, for m = 1,2,..*/
    if( y!=-1 && y!=1)
      {
	Rpp0[0]=-2.0;
	for (n=1; n<dimR; ++n)
	  {
	    Rpp0[n] = n*(y*P[n-1]-P[n-2])/(y*y-1)*(n-1)-(n+2)*P[n]-(n+2)*y*(y*P[n]-P[n-1])/(y*y-1);
	  }
      }
    else if(y==-1)
      {
	Rpp0[0]=-2.0;
	sign=-1;
	for (n=1; n< dimR; ++n)
	  {
	    sign=sign*(-1);
	    Rpp0[n]=2*sign-2*y*n*(n+1)*sign;
	  }
      }
    else if(y==1)
      {
	Rpp0[0]=-2.0;
	for (n=1; n< dimR; ++n)
	  {
	    Rpp0[n]=-2-2*y*n*(n+1);
	  }
      }
     free(P);
   return(NO_ERR);
}

/* function for the adding term used in incremental solution U(0,0) */
void CreateUadd( double *Uadd, int n, double *y)
{
  int i;
  for (i=0; i< n; ++i)
    {
      Uadd[i]=1-y[i];
    }
}

/* function for the adding term used in incremental solution v(Kx, kz) */
void CreateVadd( double *Vadd, int n, double *y)
{
  int i;
  for (i=0; i< n; ++i)
    {
      Vadd[i]=(1-y[i])*(1-y[i])*(1+y[i]);
    }
}

/* function for the adding term used in incremental solution v(Kx, kz) */
void CreateVpadd( double *Vpadd, int n, double *y)
{
  int i;
  for (i=0; i< n; ++i)
    {
      Vpadd[i]=(1-y[i])*(1-y[i])-2*(1+y[i])*(1-y[i]);
    }
}

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
*****************************************************************
Update on Oct 2: when initialize Legendre matrix, create a new vector
 Rp0, see detailed info in RsQs.c

**********************************************************************/

/* Create the matrices to be used in Gaussian quadrature and the time advance 
   procedure. */
#include "LegendreSetup.h"

int LegendreSetup(void)
{
    /* External Variables */
    extern int qpts, dimR, dimQ;
    extern double **Q, **Qp, **Qpp, **R, **Rp, **Qw, **Qpw, **Rw, **Qs, **Qps,
      **Qpps, **Rs, **Rps, *Rp0, **Rpw, **Qppw, *Rpp0;
    extern double *Qy;
    extern double *Uadd, *Vadd, *Vpadd;
    extern double *W;
    /* local variables */
    int i;
    double *y, *w;		/* quadrature points and weights */
    double *wrk;		/* work array */
    double **P, **T;

    if ( (Q = dMatrix(qpts, dimQ)) == NULL)
    {
	printf("ERROR: No memory for Q Matrix.\n");
	return(NO_MEM);
    }
    if ( (Qp = dMatrix(qpts, dimQ)) == NULL)
    {
	printf("ERROR: No memory for Qp Matrix.\n");
	freedMatrix(Q);
	return(NO_MEM);
    }
    if ( (Qpp = dMatrix(qpts, dimQ)) == NULL)
    {
	printf("ERROR: No memory for Qpp Matrix.\n");
        freedMatrix(Q); freedMatrix(Qp);
	return(NO_MEM);
    }
    if ( (R = dMatrix(qpts, dimR)) == NULL)
    {
	printf("ERROR: No memory for R Matrix.\n");
        freedMatrix(Q); freedMatrix(Qp); freedMatrix(Qpp);
	return(NO_MEM);
    }
    if ( (Rp = dMatrix(qpts, dimR)) == NULL)
    {
	printf("ERROR: No memory for Rp Matrix.\n");
        freedMatrix(Q); freedMatrix(Qp); freedMatrix(Qpp); freedMatrix(R);
	return(NO_MEM);
    }
    if ( (Qw = dMatrix(dimQ, qpts)) == NULL)
    { 
	printf("ERROR: No memory for Qw Matrix.\n");
        freedMatrix(Q); freedMatrix(Qp); freedMatrix(Qpp); freedMatrix(R);
        freedMatrix(Rp);
	return(NO_MEM);
    }
    if ( (Qpw = dMatrix(dimQ, qpts)) == NULL)
    { 
	printf("ERROR: No memory for Qpw Matrix.\n");
        freedMatrix(Q); freedMatrix(Qp); freedMatrix(Qpp); freedMatrix(R);
        freedMatrix(Rp); freedMatrix(Qw);
	return(NO_MEM);
    }
    if ( (Rw = dMatrix(dimR, qpts)) == NULL)
    { 
	printf("ERROR: No memory for Rw Matrix.\n");
        freedMatrix(Q); freedMatrix(Qp); freedMatrix(Qpp); freedMatrix(R);
        freedMatrix(Rp); freedMatrix(Qw); freedMatrix(Qpw);
	return(NO_MEM);
    }
    if ( (Qs = dMatrix(dimQ, T_QSDIAG)) == NULL)
    { 
	printf("ERROR: No memory for Qs Matrix.\n");
        freedMatrix(Q); freedMatrix(Qp); freedMatrix(Qpp); freedMatrix(R);
        freedMatrix(Rp); freedMatrix(Qw); freedMatrix(Qpw); freedMatrix(Rw);
	return(NO_MEM);
    }
    if ( (Qps = dMatrix(dimQ, T_QSDIAG)) == NULL)
    { 
	printf("ERROR: No memory for Qps Matrix.\n");
        freedMatrix(Q); freedMatrix(Qp); freedMatrix(Qpp); freedMatrix(R);
        freedMatrix(Rp); freedMatrix(Qw); freedMatrix(Qpw); freedMatrix(Rw);
        freedMatrix(Qs);
	return(NO_MEM);
    }
    if ( (Qpps = dMatrix(dimQ, T_QSDIAG)) == NULL)
    { 
	printf("ERROR: No memory for Qpps Matrix.\n");
        freedMatrix(Q); freedMatrix(Qp); freedMatrix(Qpp); freedMatrix(R);
        freedMatrix(Rp); freedMatrix(Qw); freedMatrix(Qpw); freedMatrix(Rw);
        freedMatrix(Qs); freedMatrix(Qps);
	return(NO_MEM);
    }
    if ( (Rs = dMatrix(dimR, T_RSDIAG)) == NULL)
    { 
	printf("ERROR: No memory for Rs Matrix.\n");
        freedMatrix(Q); freedMatrix(Qp); freedMatrix(Qpp); freedMatrix(R);
        freedMatrix(Rp); freedMatrix(Qw); freedMatrix(Qpw); freedMatrix(Rw);
        freedMatrix(Qs); freedMatrix(Qps); freedMatrix(Qpps);
	return(NO_MEM);
    }
    if ( (Rps = dMatrix(dimR, T_RSDIAG)) == NULL)
    { 
	printf("ERROR: No memory for Rps Matrix.\n");
        freedMatrix(Q); freedMatrix(Qp); freedMatrix(Qpp); freedMatrix(R);
        freedMatrix(Rp); freedMatrix(Qw); freedMatrix(Qpw); freedMatrix(Rw);
        freedMatrix(Qs); freedMatrix(Qps); freedMatrix(Qpps);
        freedMatrix(Rs);
	return(NO_MEM);
    }
    if ( (P = dMatrix(qpts, dimR)) == NULL)
    {
	printf("ERROR: No memory for P Matrix.\n");
        freedMatrix(Q); freedMatrix(Qp); freedMatrix(Qpp); freedMatrix(R);
        freedMatrix(Rp); freedMatrix(Qw); freedMatrix(Qpw); freedMatrix(Rw);
        freedMatrix(Qs); freedMatrix(Qps); freedMatrix(Qpps);
        freedMatrix(Rs); freedMatrix(Rps);
	return(NO_MEM);
    }
    if ( (T = dMatrix(dimR, qpts)) == NULL)
    { 
	printf("ERROR: No memory for T Matrix.\n");
        freedMatrix(Q); freedMatrix(Qp); freedMatrix(Qpp); freedMatrix(R);
        freedMatrix(Rp); freedMatrix(Qw); freedMatrix(Qpw); freedMatrix(Rw);
        freedMatrix(Qs); freedMatrix(Qps); freedMatrix(Qpps);
        freedMatrix(Rs); freedMatrix(Rps); freedMatrix(P);
	return(NO_MEM);
    }


    /* add a new vector the store the first derivatives of R at boundary y=-1*/

    if ( (Rp0 = dVector(dimR)) == NULL)
    {
	printf("ERROR: No memory for Rp0.\n");
        freedMatrix(Q); freedMatrix(Qp); freedMatrix(Qpp); freedMatrix(R);
        freedMatrix(Rp); freedMatrix(Qw); freedMatrix(Qpw); freedMatrix(Rw);
        freedMatrix(Qs); freedMatrix(Qps); freedMatrix(Qpps);
        freedMatrix(Rs); freedMatrix(Rps); freedMatrix(P); freedMatrix(T);
	return(NO_MEM);
    }

    if ( (w = dVector(qpts*3)) == NULL)
    {
	printf("ERROR: No memory for quad points and weights.\n");
        freedMatrix(Q); freedMatrix(Qp); freedMatrix(Qpp); freedMatrix(R);
        freedMatrix(Rp); freedMatrix(Qw); freedMatrix(Qpw); freedMatrix(Rw);
        freedMatrix(Qs); freedMatrix(Qps); freedMatrix(Qpps);
        freedMatrix(Rs); freedMatrix(Rps); freedMatrix(P); freedMatrix(T);
        freedVector(Rp0); 
	return(NO_MEM);
    }

    if ( (Uadd = dVector(qpts))==NULL)
    {
	printf("ERROR: No memory for adding term in the incremental solution.\n");
        freedMatrix(Q); freedMatrix(Qp); freedMatrix(Qpp); freedMatrix(R);
        freedMatrix(Rp); freedMatrix(Qw); freedMatrix(Qpw); freedMatrix(Rw);
        freedMatrix(Qs); freedMatrix(Qps); freedMatrix(Qpps);
        freedMatrix(Rs); freedMatrix(Rps); freedMatrix(P); freedMatrix(T);
        freedVector(Rp0); freedVector(w);
	return(NO_MEM);
    }

    if (( Vadd=dVector(qpts))==NULL)
      {
	printf("ERROR: No memory for adding term in the incremental solution.\n");
        freedMatrix(Q); freedMatrix(Qp); freedMatrix(Qpp); freedMatrix(R);
        freedMatrix(Rp); freedMatrix(Qw); freedMatrix(Qpw); freedMatrix(Rw);
        freedMatrix(Qs); freedMatrix(Qps); freedMatrix(Qpps);
        freedMatrix(Rs); freedMatrix(Rps); freedMatrix(P); freedMatrix(T);
        freedVector(Rp0); freedVector(w); freedVector(Uadd);
	return(NO_MEM);
    }

    if (( Vpadd=dVector(qpts))==NULL)
      {
	printf("ERROR: No memory for adding term in the incremental solution.\n");
        freedMatrix(Q); freedMatrix(Qp); freedMatrix(Qpp); freedMatrix(R);
        freedMatrix(Rp); freedMatrix(Qw); freedMatrix(Qpw); freedMatrix(Rw);
        freedMatrix(Qs); freedMatrix(Qps); freedMatrix(Qpps);
        freedMatrix(Rs); freedMatrix(Rps); freedMatrix(P); freedMatrix(T);
        freedVector(Rp0); freedVector(w); freedVector(Vadd);freedVector(Uadd);
	return(NO_MEM);
    }
    
    if((Qy=dVector(qpts))==NULL)
      {
	printf("ERROR: No memory for quadrature points.\n");
        freedMatrix(Q); freedMatrix(Qp); freedMatrix(Qpp); freedMatrix(R);
        freedMatrix(Rp); freedMatrix(Qw); freedMatrix(Qpw); freedMatrix(Rw);
        freedMatrix(Qs); freedMatrix(Qps); freedMatrix(Qpps);
        freedMatrix(Rs); freedMatrix(Rps); freedMatrix(P); freedMatrix(T);
        freedVector(Rp0); freedVector(w); freedVector(Vadd);freedVector(Uadd);
	freedVector(Vpadd);
	return(NO_MEM);
    }
    if((W=dVector(qpts))==NULL)
      {
	printf("ERROR: No memory for quadrature points.\n");
        freedMatrix(Q); freedMatrix(Qp); freedMatrix(Qpp); freedMatrix(R);
        freedMatrix(Rp); freedMatrix(Qw); freedMatrix(Qpw); freedMatrix(Rw);
        freedMatrix(Qs); freedMatrix(Qps); freedMatrix(Qpps);
        freedMatrix(Rs); freedMatrix(Rps); freedMatrix(P); freedMatrix(T);
        freedVector(Rp0); freedVector(w); freedVector(Vadd);freedVector(Uadd);
	freedVector(Vpadd);
	return(NO_MEM);
    }
    y = w + qpts;
    wrk = y + qpts; 

    if ( (Rpw = dMatrix(dimR, qpts)) == NULL)
      { 
	printf("ERROR: No memory for Rpw.\n");
        freedMatrix(Q); freedMatrix(Qp); freedMatrix(Qpp); freedMatrix(R);
        freedMatrix(Rp); freedMatrix(Qw); freedMatrix(Qpw); freedMatrix(Rw);
        freedMatrix(Qs); freedMatrix(Qps); freedMatrix(Qpps);
        freedMatrix(Rs); freedMatrix(Rps); freedMatrix(P); freedMatrix(T);
        freedVector(Rp0); freedVector(w);freedVector(Vadd);freedVector(Vpadd);
	freedVector(Uadd); freedVector(Qy); freedVector(W);
	return(NO_MEM);
      }
    if ( (Qppw = dMatrix(dimQ, qpts)) == NULL)
      { 
	printf("ERROR: No memory for Qppw.\n");
        freedMatrix(Q); freedMatrix(Qp); freedMatrix(Qpp); freedMatrix(R);
        freedMatrix(Rp); freedMatrix(Qw); freedMatrix(Qpw); freedMatrix(Rw);
        freedMatrix(Qs); freedMatrix(Qps); freedMatrix(Qpps);
        freedMatrix(Rs); freedMatrix(Rps); freedMatrix(P); freedMatrix(T);
        freedVector(Rp0); freedVector(w);freedVector(Vadd);freedVector(Vpadd);
	freedVector(Uadd); freedVector(Qy); freedVector(W); freedMatrix(Rpw);
	return(NO_MEM);
      }


  if ( (Rpp0 = dVector(dimR)) == NULL)
    {
	printf("ERROR: No memory for Rpp0.\n");
        freedMatrix(Q); freedMatrix(Qp); freedMatrix(Qpp); freedMatrix(R);
        freedMatrix(Rp); freedMatrix(Qw); freedMatrix(Qpw); freedMatrix(Rw);
        freedMatrix(Qs); freedMatrix(Qps); freedMatrix(Qpps);
        freedMatrix(Rs); freedMatrix(Rps); freedMatrix(P); freedMatrix(T);
	freedVector(Rp0); freedVector(w);freedVector(Vadd);freedVector(Vpadd);
	freedVector(Uadd); freedVector(Qy);
	return(NO_MEM);
    }
    /* Compute quadrature points */
    if ( (i = Qpoints(qpts, y, qpts, Legendre)) != NO_ERR )
    {
	printf("ERROR: Qpoints error code %d.\n", i);
        freedMatrix(Q); freedMatrix(Qp); freedMatrix(Qpp); freedMatrix(R);
        freedMatrix(Rp); freedMatrix(Qw); freedMatrix(Qpw); freedMatrix(Rw);
        freedMatrix(Qs); freedMatrix(Qps); freedMatrix(Qpps);
        freedMatrix(Rs); freedMatrix(Rps); freedMatrix(P); freedMatrix(T);
        freedVector(Rp0); freedVector(w);freedVector(Vadd);freedVector(Vpadd);
	freedVector(Uadd); freedVector(Qy);freedVector(Rpp0);
	return(i);
    }

    memcpy(Qy, y, qpts*sizeof(double));
  
    //  printf("the Quadrature points are: \n");
    /* for (i=0; i<qpts; ++i)
      {
	printf(" %f\n", Qy[i]);
	}*/

    /* Compute quadrature weights and create the matrices */
    Qweights(qpts, y, w, Legendre);
    memcpy(W, w, qpts*sizeof(double));
    LegendreMatrix(P, qpts, dimR, y);

    /* create several new vector used in incremental computation */
    CreateUadd(Uadd, qpts, Qy);
    CreateVadd(Vadd, qpts, Qy);
    CreateVpadd(Vpadd, qpts, Qy);

    CreateR(R, P, qpts, dimR, y, wrk);
    CreateQ(Q, R, qpts, dimQ, y, wrk); 
    CreateQp(Qp, P, qpts, dimQ, y, wrk); 
    CreateQpp(Qpp, P, qpts, dimQ, y, wrk); 
    CreateRp(Rp, P, qpts, dimR, y); 
    CreateRp0(Rp0, dimR, -1.0);             /* new added vector */
    CreateRpp0(Rpp0, dimR, -1.0);

    weightsToM(Qw, Q, dimQ, qpts, w); 
    weightsToM(Qpw, Qp, dimQ, qpts, w); 
    weightsToM(Qppw, Qpp, dimQ, qpts, w); 
    weightsToM(Rw, R, dimR, qpts, w); 
     weightsToM(Rpw, Rp, dimR, qpts, w); 
    LdiagsMatrix(Qs, dimQ, T_QSDIAG, 0, QSDIAG, Qw, qpts, Q);
    LdiagsMatrix(Qps, dimQ, T_QSDIAG, 1, QPSDIAG, Qpw, qpts, Qp);
    weightsToM(T, Qpp, dimQ, qpts, w); 
    LdiagsMatrix(Qpps, dimQ, T_QSDIAG, 2, QPPSDIAG, T, qpts, Qpp);

    LdiagsMatrix(Rs, dimR, T_RSDIAG, 0, RSDIAG, Rw, qpts, R);
    weightsToM(T, Rp, dimR, qpts, w); 
    LdiagsMatrix(Rps, dimR, T_RSDIAG, 1, RPSDIAG, T, qpts, Rp);

    freedVector(w);
    freedMatrix(P); 
    freedMatrix(T); 

    return(NO_ERR);
}


/* function to evaluate the Legendre polynomials */
/* Parameters:  n : degree of the polynomial
		x : point at which the poly is to be evaluated
   Return value: Legendre polynomial of degree n evaluated at x 
*/
double Legendre(int n, double x) 
{
    int j;
    double y, prev, curr;

    prev = 0.;
    curr = 1.;
    for (j=0; j<n; ++j)
    {
        y = ((double)(2*j+1)*x*curr - (double)j*prev) / (double)(j+1);
        prev = curr;
        curr = y;
    }
    return(y);
}


/* Compute Gauss quadrature points.
   (I translated this into C from a routine that was written in vectoral.
   Note from vectoral code: this is phillip's initial guess, search through the
   Chebychev roots to find crossings, start with 4*n chebychev points, use more 
   if needed store bounds of each root.)
Parameters:
   n : degree of the polynomial
   x : 1D array to store the points in
   len : size of the array x
   poly : pointer to the function that evaluates the polynomial.  The arguments
      to poly are: the degree of the polynomial and the value at which to  
      evaluate the polynomial.
Return value: code indicating successful completion or no convergence.
*/      
int Qpoints(int n, double *x, int len, double (*poly)(int, double))
{
    int k, m, r, nd2;
    double *xe, p1, p2, p3, x1, x2, x3, test;

    k = n; m = 0;
    nd2 = n/2;

    if ( (xe = (double *)calloc(n, sizeof(double))) == NULL ) 
    {
	return(NO_MEM);
    }

    while (m < nd2)
    {
        k = 2*k;
        x1 = -1.0;
        p1 = (*poly)(n, x1);
        m = 0;
        for (r=1; r<=k; ++r)
        {
            x2 = -cos(0.5*M_PI*(r-0.5)/k);    	/* ith CHEBYSHEV POINT */
            p2 = (*poly)(n, x2);
            if (p1*p2 < 0) 
                {  x[m] = x1;  xe[m] = x2;  m = m+1;  }
            p1 = p2;   x1 = x2;
        }
    }
    
    /* loop through half the points iterating to convergence */
    for (m=0; m<nd2; ++m)
    {
        x1 = x[m];
        x2 = xe[m];
        p1 = (*poly)(n, x1);
        p2 = (*poly)(n, x2);
        r = 1;
    
        /* this is a modified secant iteration to avoid getting stuck */
        test = 1.0 + (x2-x1)/TOL;
        while  (test != 1.0  &&  r < ITMAX)
	{
          if  (r%10 == 0)
	  {
            x3 = 0.5 * (x1+x2);
	  }
          else
	  {
	    if  (fabs(p1) < fabs(p2))
	      {  x3 = x1 + MAX(0.05,p1/(p1-p2))*(x2-x1);  }
	    else
	      {  x3 = x2 + MIN(-0.05,p2/(p1-p2))*(x2-x1);  }
          }
          p3 = (*poly)(n, x3);
          if  (p3*p1 <= 0)
	  {
            x2 = x3;	 p2 = p3;
	  }
          else
	  {
            x1 = x3;	 p1 = p3;
          }
          r = r+1;
          test = 1.0 + (x2-x1)/TOL;	/* not in vectoral code, but... */
        }     /* while test... */
        if  (r == ITMAX)
	  {  return(NO_CONVERGENCE);  }

        x[m] = x3;
        x[n-1-m] = -x3;
    }   /* for m = 1 to nd2  */
    
    if (n%2 == 1)
      {  x[nd2] = 0.0;  }

    free(xe);

    return(NO_ERR);
}


/* Create a rows-by-cols matrix with entries P(m,n) = P_n(x_m), where P_n is
 the nth Legendre polynomial and x_m is the mth value at which the polynomial
 must be evaluated.
 Parameters
	P    : ptr to double array
	rows : number of rows for P
	cols : number of columns for P
	x    : rows-by-1 array of values for evaluating the polynomials.

 Recurrence relation for Legendre Polynomials
 	P_0(x)    = 1
 	P_1(x)    = x
 	P_n+1(x)  = (2n+1)/(n+1) xP_n  -  n/(n+1) P_n-1,  n=1,2,...
*/
void LegendreMatrix(double **P, int rows, int cols, double *x)
{
    int m, n;

    for (m=0; m<rows; ++m)
    {
	P[m][0] = 1.0;
	P[m][1] = x[m];
        for (n=1; n<(cols-1); ++n)
	{
            P[m][n+1]=((2.*n+1.)*(x[m]*P[m][n]) - n*P[m][n-1]) / (double)(n+1);
	}
    }
}


/* Create a symmetric, banded matrix D with entries D(i,j)=<A(m),B(n)>, for 
   appropriate row m of A, column n of B.  D has a zero diagonal between each
   nonzero diagonal and the number of nonzero diagonals in the upper part of
   D is given by the parameter diags.
   Parameters
	D     : n-by-m double array for storing the resulting matrix
	n     : number of rows in D
	m     : number of columns in D.  Should equal 2*diags-1...
	start : used to identify which column of D will contain the main 
		diagonal (and, thus, the other nonzero diagonals).
		0 => main diagonal is in column diags-1 (and D will not contain
		     any zero colums)
		1 => main diagonal is in column diags (and D will contain
		     two zero colums)
		2 => main diagonal is in column diags+1 (and D will contain
		     four zero colums)
                This is done to simplify the rest of the code...
	diags : number of nonzero diagonals in upper part of D
	Adim2 : number of colums in matrix A
	A     : n-by-Adim2 double array
	B     : Adim2-by-n double array
*/
void LdiagsMatrix(double **D, int n, int m, int start, int diags, double **A,
    int Adim2, double **B)
{
    int i, j, k, s;

    for (i=0; i<n; ++i)
    {
        for (j=0; j<m; ++j)
	{
	    D[i][j] = 0.0;
	}
    }

    for (i=start, s=0; s<diags; ++i, ++s)
    {
        for (j=0; j<(n-2*s); ++j)
	{
            for (k=0; k<Adim2; ++k)
	    {
		D[j][diags-1+i] = D[j][diags-1+i] + A[j][k]*B[k][j+2*s]; 
	    }
	    if (i != start)
	    {
	        D[2*s+j][diags-1+i-2*s] = D[j][diags-1+i]; 
	    }
	}
    }
}


/* Compute Gauss quadrature weights.  
   Parameters
 	n    : number of quadrature points
 	x    : n x 1 array containing the quadrature points
 	w    : n x 1 array to store the weights
 	poly : string containing the name of the function that evaluates
 	       the polynomial.  Its parameters are the degree of the polynomial
	       and the value at which to evaluate the polynomial.
*/
void Qweights(int n, double *x, double *w, double (*poly)(int, double))
{
    int lim, k;
    double v1;

    lim = (n+1)/2;

    for (k=0; k<lim; ++k)
    {
	v1 = (double)n * ((*poly)(n-1,x[k]) - x[k]*(*poly)(n,x[k]));
        w[k] = 2.*(1.-x[k]*x[k]) / (v1*v1);
	w[n-1-k] = w[k];
    }
}


/* Create matrix Mw with entries Mw(i,j) = w(j)*M(j,i).
   Parameters
	Mw   : double array for storing the resulting matrix
	dim1 : number of rows in Mw
	dim2 : number of columns in Mw
	M    : dim2-by-dim1 double array used in contructing Mw
	w    : dim2 double vector with the weights
*/
void weightsToM(double **Mw, double **M, int dim1, int dim2, double *w)
{
    int i, j;

    for (i=0; i<dim1; ++i)
    {
        for (j=0; j<dim2; ++j)
	{
	    Mw[i][j] = w[j]*M[j][i]; 
	}
    }
}

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
#include <stdlib.h> 
#include <math.h> 
#include "minChnl.h"
#include "arrays.h"

int waveNums(int x, int z, double Lx, double Lz)
{
    /* External Variables */
    extern double *Kx, *Kz, **K2;

    int i, j;

    if ( (Kx = dVector(x)) == NULL)
    {
	printf("ERROR: No memory for Kx array.\n");
	return(NO_MEM);
    }
    if ( (Kz = dVector(z)) == NULL)
    {
	printf("ERROR: No memory for Kz array.\n");
	freedVector(Kx);
	return(NO_MEM);
    }
    if ( (K2 = dMatrix(z, x)) == NULL)      
    { 
	printf("ERROR: No memory for K2 Matrix.\n");
	freedVector(Kx); freedVector(Kz);
	return(NO_MEM);
    }

    for (i = 0; i < x; ++i)
        {   Kx[i] = i*(2.*M_PI/Lx);   }
    for (i = 0; i < z/2+1; ++i)
        {   Kz[i] = i*(2.*M_PI/Lz);   }
    for (i = z/2+1; i < z; ++i)
        {   Kz[i] = (i-z)*(2.*M_PI/Lz);   }

    for (i = 0; i < z; ++i)
    {
        for (j = 0; j < x; ++j)
	    {   K2[i][j] = Kz[i]*Kz[i] + Kx[j]*Kx[j];   }
    }

    return(NO_ERR);
}


int cflVars(double Lx, double Lz)
{
    double Legendre(int, double);
    int Qpoints(int n, double *x, int len, double (*poly)(int, double));

    extern int    qpts, Nz, Nx;
    extern double cfl1, cfl3;
    extern double *cfl2;

    int j;
    double *y;

    if ( (y = dVector(qpts)) == NULL )
    {
	printf("INITIAL FIELD: No memory for quad points.\n");
	return(NO_MEM);
    }
    if ( (j = Qpoints(qpts, y, qpts, Legendre)) != NO_ERR )
    {
	printf("INITIAL FIELD: Qpoints error code %d.\n", j);
	return(j);
    }
  
    cfl1 = M_PI / (Lx/Nx);
    cfl3 = M_PI / (Lz/Nz);
    for (j = 1; j < qpts-1; ++j)
    {
	cfl2[j] = 2.*M_PI / (y[j+1]-y[j-1]);
    }
    cfl2[0] = cfl2[qpts-1] = M_PI / (y[0]+1.0);

    freedVector(y);

    return(0);
}     /* end cflVars */

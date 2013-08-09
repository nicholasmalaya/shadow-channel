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

#ifndef ARRAYS_H
#define ARRAYS_H

#include "mcomplex.h"

/* allocate space for a vector of doubles; complex numbers */
double *dVector(int len);
mcomplex *cVector(int len);

/* allocate space for a matrix of doubles; complex numbers */
double **dMatrix(int dim1, int dim2);
mcomplex **cMatrix(int dim1, int dim2);

/* allocate space for a 3D array of doubles; complex numbers */
double ***d3Darray(int dim1, int dim2, int dim3);
mcomplex ***c3Darray(int dim1, int dim2, int dim3);

/* allocate space for a 4D array of doubles; complex numbers */
double ****d4Darray(int dim1, int dim2, int dim3, int dim4);
mcomplex ****c4Darray(int dim1, int dim2, int dim3, int dim4);

/* allocate space for a 5D array of doubles; complex numbers */
mcomplex *****c5Darray(int dim1, int dim2, int dim3, int dim4, int dim5);

/* free space allocated */
void freedVector(double *v);
void freecVector(mcomplex * v);
void freedMatrix(double **m);
void freecMatrix(mcomplex ** m);
void freed3Darray(double ***m);
void freec3Darray(mcomplex *** m);
void freed4Darray(double ****m);
void freec4Darray(mcomplex **** m);
void freec5Darray(mcomplex ***** m);
#endif

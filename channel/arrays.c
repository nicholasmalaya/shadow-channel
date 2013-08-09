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

#include <stdlib.h>
#include "mcomplex.h"

/* allocate space for a vector of doubles; complex numbers */
double *dVector(int len)
{
    double *v = (double *) calloc(len, sizeof(double));
    return (v);
}

mcomplex *cVector(int len)
{
    mcomplex *v = (mcomplex *) calloc(len, sizeof(mcomplex));
    return (v);
}

/* allocate space for a matrix of doubles; mcomplex numbers */
double **dMatrix(int dim1, int dim2)
{
    int i;
    double **m;
    if ((m = (double **) malloc(dim1 * sizeof(double *))) == NULL) {
        return (NULL);
    }
    if ((m[0] = (double *) calloc(dim1 * dim2, sizeof(double))) == NULL) {
        free(m);
        return (NULL);
    }

    for (i = 1; i < dim1; ++i) {
        m[i] = m[i - 1] + dim2;
    }
    return (m);
}

mcomplex **cMatrix(int dim1, int dim2)
{
    int i;
    mcomplex **m;
    if ((m = (mcomplex **) malloc(dim1 * sizeof(mcomplex *))) == NULL) {
        return (NULL);
    }
    if ((m[0] =
         (mcomplex *) calloc(dim1 * dim2, sizeof(mcomplex))) == NULL) {
        free(m);
        return (NULL);
    }

    for (i = 1; i < dim1; ++i) {
        m[i] = m[i - 1] + dim2;
    }
    return (m);
}

/* allocate space for a 3D array of doubles; complex numbers */
double ***d3Darray(int dim1, int dim2, int dim3)
{
    int i;
    double ***m;
    if ((m = (double ***) malloc(dim1 * sizeof(double **))) == NULL) {
        return (NULL);
    }
    if ((m[0] = dMatrix(dim1 * dim2, dim3)) == NULL) {
        free(m);
        return (NULL);
    }
    for (i = 1; i < dim1; ++i) {
        m[i] = m[i - 1] + dim2;
    }
    return (m);
}

mcomplex ***c3Darray(int dim1, int dim2, int dim3)
{
    int i;
    mcomplex ***m;
    if ((m = (mcomplex ***) malloc(dim1 * sizeof(mcomplex **))) == NULL) {
        return (NULL);
    }
    if ((m[0] = cMatrix(dim1 * dim2, dim3)) == NULL) {
        free(m);
        return (NULL);
    }
    for (i = 1; i < dim1; ++i) {
        m[i] = m[i - 1] + dim2;
    }
    return (m);
}

/* allocate space for a 4D array of doubles; complex numbers */
double ****d4Darray(int dim1, int dim2, int dim3, int dim4)
{
    int i;
    double ****m;
    if ((m = (double ****) malloc(dim1 * sizeof(double ***))) == NULL) {
        return (NULL);
    }
    if ((m[0] = d3Darray(dim1 * dim2, dim3, dim4)) == NULL) {
        free(m);
        return (NULL);
    }
    for (i = 1; i < dim1; ++i) {
        m[i] = m[i - 1] + dim2;
    }
    return (m);
}

mcomplex ****c4Darray(int dim1, int dim2, int dim3, int dim4)
{
    int i;
    mcomplex ****m;
    if ((m = (mcomplex ****) malloc(dim1 * sizeof(mcomplex ***))) == NULL) {
        return (NULL);
    }
    if ((m[0] = c3Darray(dim1 * dim2, dim3, dim4)) == NULL) {
        free(m);
        return (NULL);
    }
    for (i = 1; i < dim1; ++i) {
        m[i] = m[i - 1] + dim2;
    }
    return (m);
}

mcomplex *****c5Darray(int dim1, int dim2, int dim3, int dim4, int dim5)
{
    int i;
    mcomplex *****m;
    if ((m =
         (mcomplex *****) malloc(dim1 * sizeof(mcomplex ****))) == NULL) {
        return (NULL);
    }
    if ((m[0] = c4Darray(dim1 * dim2, dim3, dim4, dim5)) == NULL) {
        free(m);
        return (NULL);
    }
    for (i = 1; i < dim1; ++i) {
        m[i] = m[i - 1] + dim2;
    }
    return (m);
}

/* free space allocated */
void freedVector(double *v)
{
    free(v);
}

void freecVector(mcomplex * v)
{
    free(v);
}

void freedMatrix(double **m)
{
    free(m[0]);
    free(m);
}

void freecMatrix(mcomplex ** m)
{
    free(m[0]);
    free(m);
}

void freed3Darray(double ***m)
{
    freedMatrix(m[0]);
    free(m);
}

void freec3Darray(mcomplex *** m)
{
    freecMatrix(m[0]);
    free(m);
}

void freed4Darray(double ****m)
{
    freed3Darray(m[0]);
    free(m);
}

void freec4Darray(mcomplex **** m)
{
    freec3Darray(m[0]);
    free(m);
}

void freec5Darray(mcomplex ***** m)
{
    freec4Darray(m[0]);
    free(m);
}

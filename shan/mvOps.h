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

#ifndef MVOPS_H
#define MVOPS_H

#include "mcomplex.h"

void bsolve(double ***A, mcomplex **b, int dl, int du, int N, int xdim, 
    int x0);
void bsolve0(double **A, mcomplex **b, int dl, int du, int N);
void bsolve_0(double **A, mcomplex *b, int dl, int du, int N);
void smMult(double ***A, mcomplex **y, mcomplex **b, int dl, int du, int N,
    int xdim, int x0);
void smMult0(double **A, mcomplex **y, mcomplex *b, int dl, int du, int N);

#endif

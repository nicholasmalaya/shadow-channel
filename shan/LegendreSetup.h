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

#ifndef LEGENDRESETUP_H
#define LEGENDRESETUP_H

#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 
#include "minChnl.h"
#include "arrays.h"
#include "RsQs.h"

/* function declarations */
void weightsToM(double **, double **, int, int, double *);
void LdiagsMatrix(double **, int, int, int, int, double **, int, double **);
double Legendre(int, double);
int Qpoints(int n, double *x, int len, double (*poly)(int, double));
void Qweights(int n, double *x, double *w, double (*poly)(int, double));
void LegendreMatrix(double **P, int rows, int cols, double *x);

#endif  /* LEGENDRESETUP_H */

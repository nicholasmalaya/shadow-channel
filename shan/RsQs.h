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

#ifndef RSQS_H
#define RSQS_H

void CreateR(double **R, double **P, int dim1, int dim2, double *y, 
    double *wrk);
void CreateRp(double **Rp, double **P, int dim1, int dim2, double *y);
void CreateQ(double **Q, double **R, int dim1, int dim2, double *y, 
    double *wrk);
void CreateQp(double **Qp, double **P, int dim1, int dim2, double *y,
    double *wrk);
void CreateQpp(double **Qpp, double **P, int dim1, int dim2, double *y,
    double *wrk);
int CreateRp0(double *Rp0, int dimR, double y);
int CreateRpp0(double *Rpp0, int dimR, double y);
void CreateUadd( double *Uadd, int n, double *y);
void CreateVadd( double *Vadd, int n, double *y);
void CreateVpadd( double *Vpadd, int n, double *y);
#endif   /* RSQS_H */

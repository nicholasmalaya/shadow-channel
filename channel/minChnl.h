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

#ifndef MINCHNL_H
#define MINCHNL_H

#include "mcomplex.h"

/* Error codes */
#define NO_ERR          0
#define ERROR          -1
#define NO_MEM         -2
#define NO_CONVERGENCE -3

/* Parameters for Qpoints function */
#define ITMAX  200
#define TOL    1E-15

/* Other defines */
#define CFL_MAX  2.8

/* diagonals for Legendre matrices (in terms of upper triangle of matrix) */
#define QSDIAG    5
#define QPSDIAG   4
#define QPPSDIAG  3
#define RSDIAG    3
#define RPSDIAG   2
/* total */
#define T_QSDIAG    (2*QSDIAG - 1)
#define T_RSDIAG    (2*RSDIAG - 1)

/* indices for 4D arrays */
/* array U */
#define XEL    0
#define YEL    1
#define ZEL    2
#define DXEL   3
#define DZEL   4
#define MAXU   5
#define HXEL   1
#define HYEL   2
#define HZEL   0
#define WXEL   4
#define WYEL   5
#define WZEL   3
#define MAXT   6
#define MAXTT  9
#define AXEL   6
#define AYEL   7
#define AZEL   8
#define SXEL   1
#define SYEL   2
#define SZEL   0
#define TXEL   5
#define TYEL   3
#define TZEL   4
#define MAXSTEP 100
#define FT  0.5
/* array C */
#define ALPHA  0
#define BETA   1

#define MAX(x,y)  ((x)>(y) ? (x):(y))
#define MIN(x,y)  ((x)<(y) ? (x):(y))
#define MAGNITUDE(x)  (sqrt(Re((x))*Re((x)) + Im((x))*Im((x))))
typedef void (*func_force_t) (int, int, int, mcomplex *, mcomplex *);

#endif                          /* MINCHNL_H */

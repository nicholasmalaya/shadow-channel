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

#ifndef CHANNEL_MAIN_H
#define CHANNEL_MAIN_H

#include <math.h>
#include <stdio.h>
#include <string.h>
#include "fftw.h"
#include "rfftw.h"
#include "hdf5.h"
#include "minChnl.h"
#include "arrays.h"
#include "myforce.h"
#include "myadjointforce.h"

/* function prototypes */
int LegendreSetup(void);
int waveNums(int x, int z, double Lx, double Lz);
int cflVars(double Lx, double Lz);
int getMem(void);
void project0(int k, int n, func_force_t force);
void project(int n, int k, int z, int x0, func_force_t force);
void adjproject(int n, int k, int z, int x0, int count,
                func_force_t force);
void adjproject0(int k, int n, int count, func_force_t force);
void increproject(int k, int z, int x0, int n,
                  func_force_tt force);
void increproject0(int k, int n, int flag, func_force_t force);
void increadjproject0(int k, int n, int count, func_force_t force);
void increadjproject(int n, int k, int z, int x0, int count,
                     func_force_t force);
int comp_stat(double **us);


FILE *fp, *fp2;
FILE *fp3, *fp4, *fp5, *fp6, *fp7, *fp8, *fp9, *fp10, *fp11;

int pass1(int dctr, int n);
int pass2(int dctr, int n);
int increBoundary();
int increBoundary_s();
int InitialField1(double, double);
int InitialField2(double, double);
int InitialField3(double, double);
int InitialField3b(double, double);
int InitialField4(double, double);
int InitialField4b(double, double);
int InitialField5(double, double);
double compute_err(mcomplex **** A, mcomplex **** B);
double compute_err2(mcomplex **** A, mcomplex **** B, int n, int k);
void initAlphaBeta(void);
void initAlphaBeta2(void);
void incre_initAlphaBeta(void);
void incre_initAlphaBeta2(void);
void restart(int restart_flag);
void restart2(int restart_flag);
int write_data2(int n);
int write_data(int n);
int comp_gradient(int n, int k);
int comp_hess(int n, int k);
int comp_stats(double **us);

/* External Variables */

int qpts,                       /* 3Ny/2, Ny specified by user */
 dimR,                          /* Ny-2.  # of columns of matrix R(i,j) = (1-y(i)^2)P_j(y(i)) */
 dimQ,                          /* Ny-4.  # of cols of matrix Q(i,j) = (1-y(i)^2)^2 P_j(y(i)) */
 Nx,                            /* specified by user */
 Nz;                            /* specified by user */

int nsteps;                     /* number of time steps that are saved in memory */

double dt,                      /* time step */
 re,                            /* 1/RE, Reynolds number RE specified by user */
 mpg,                           /* mean pressure gradient (used when (Kx,Kz)=(0,0)) */
 flux,                          /* total mass flux in the domain. */
 cfl1,                          /* constant for cfl computation in x */
 cfl3;                          /* constant for cfl computation in z */
double *Kx,                     /* Nx/2 vector of wave numbers in x-dir */
*Kz,                            /* Nz vector of wave numbers in z-dir */
**K2,                           /* Nz-by-(Nx/2) matrix for sum of squares of wave numbers */
*cfl2;                          /* 3Ny/2 vector for cfl computation in y */
double **Q,                     /* qpts-by-dimQ matrix. Q(i,j) = (1-y(i)^2)^2 P_j(y(i)) */
**Qp,                           /* qpts-by-dimQ matrix. 1st derivatives of Q */
**Qpp,                          /* qpts-by-dimQ matrix. 2nd derivatives of Q */
**R,                            /* qpts-by-dimR matrix. R(i,j) = (1-y(i)^2)P_j(y(i)) */
**Rp,                           /* qpts-by-dimR matrix. 1st derivatives of R */
/* new defined variable, there is no corresponding vector for Q, because the first 
   order derivative of Q evaluated at y=-1 is always zero */
*Rp0,                           /*dimR vector.  1st derivative of R evaluated at boundary y=-1 */
*Rpp0,                          /*dimR vector.  2nd derivative of R evaluated at boundary y=-1 */
**Qw,                           /* dimQ-by-qpts matrix. Includes weights for Gaussian quad */
**Qpw,                          /* dimQ-by-qpts matrix. Includes weights for Gaussian quad */
**Qppw,                         /* dimQ-by-qpts matrix. Includes weights for Gaussian quad */
**Rw,                           /* dimR-by-qpts matrix. Includes weights for Gaussian quad */
**Rpw,                          /* dimR-by-qpts matrix. Includes weights for Gaussian quad */
**Qs,                           /* dimQ-by-9 matrix.  Qs(i,j) has integral(-1,1) of Q_i*Q_j. 
                                   This matrix is sparse with 9 nonzero diagonals: there is
                                   a diagonal of all zeros between each nonzero diagonal. */
**Qps,                          /* dimQ-by-9 matrix.  Qps(i,j) has integral(-1,1) of Qp_i*Qp_j. 
                                   This matrix is sparse with 7 nonzero diagonals, but I store
                                   two diagonals with zeros for simplicity.  If this becomes a
                                   problem, the project and Legendre setup routines must be
                                   modified */
**Qpps,                         /* dimQ-by-9 matrix.  Qpps(i,j) has integral(-1,1) of 
                                   Qpp_i*Qpp_j. This matrix is sparse with 5 nonzero diagonals,
                                   but I store diagonals with zeros for simplicity.  Same 
                                   comment as above. */
**Rs,                           /* dimR-by-5 matrix.  Rs(i,j) has integral(-1,1) of R_i*R_j. 
                                   This matrix is sparse with 5 nonzero diagonals: there is
                                   a diagonal of all zeros between each nonzero diagonal. */
**Rps;                          /* dimR-by-5 matrix.  Rps(i,j) has integral(-1,1) of 
                                   Rp_i*Rp_j. This matrix is sparse with 3 nonzero diagonals,
                                   but I store diagonals with zeros for simplicity.  Same 
                                   comment as above. */

double *Uadd, *Vadd, *Vpadd;    /* qpts vectors. new vector added used in incremental state system
                                   for the lift function to satisfy the boundary
                                   conditions. */
double *Qy, *W;

double ***M;                    /* (Ny-2)-by-9-by-(Nx/2). Used in the discretization process. */
double **MZ;                    /* (Ny-2)-by-5. Used in the discretization process. */
mcomplex ****U,                 /* Nz-by-5-by-(3Ny/2)-by-(Nx/2).  Used to store ux_hat, uy_hat,
                                   uz_hat, dux_hat, duz_hat. */
****C,                          /* Nz-by-2-by-(Ny-2)-by-(Nx/2).  Used for alphas and betas. */
****IU,                         /* Nz-by-5-by-(3Ny/3)-by-(Nx/2). used to store incremental solutions
                                   Iux_hat, Iuy_hat, Iuz_hat, dIux_hat, dIuz_hat. */
****IC,                         /*Nz-by-2-by-(Ny-2)-by-(Nx/2) used for alphas and betas in the incremenal system. */
****AU,                         /* Nz-by-5-by-(3Ny/2)-by-(Nx/2).  Used to store adjoint solution \lambda_x_hat, \lambda_y_hat,
                                   \lambda_z_hat, d\lambda_x_hat, d\lambda_z_hat. */
****AC,                         /* Nz-by-2-by-(Ny-2)-by-(Nx/2).  Used for alphas and betas. */
****IAU,                        /* Nz-by-5-by-(3Ny/3)-by-(Nx/2). used to store incremental solutions
                                   Iux_hat, Iuy_hat, Iuz_hat, dIux_hat, dIuz_hat. */
****IAC, ****LU,                /* Nz-by-5-by-(3Ny/2)-by-(Nx/2).  Used to store ux_hat, uy_hat,
                                   uz_hat, dux_hat, duz_hat for previous time stage, used in adjoint solver */
****LIU;                        /* Nz-by-5-by-(3Ny/3)-by-(Nx/2). used to store incremental solutions
                                   Iux_hat, Iuy_hat, Iuz_hat, dIux_hat, dIuz_hat at previous time stage. */
mcomplex *****MC,               /*(3T+1)* Nz-by-2-by-(Ny-2)-by-(Nx/2).  Used for alphas and betas for adjoint computations */
*****MIC;                       /*(3T+1)* Nz-by-2-by-(Ny-2)-by-(Nx/2).  Used for incremental  alphas and betas for adjoint computations */

mcomplex **Fa,                  /* (Ny-2)-by-(Nx/2).  Used for Fs computed with alphas. */
**Fb,                           /* (Ny-2)-by-(Nx/2).  Used for Fs computed with betas. */
**TM;                           /* (Ny-2)-by-(Nx/2).  Work array (SEE IF I CAN GET AWAY 
                                   WITHOUT IT...? */
mcomplex **test;
mcomplex **IFa,                 /* (Ny-2)-by-(Nx/2).  Used for Fs computed with alphas in the incremental system. */
**IFb,                          /* (Ny-2)-by-(Nx/2).  Used for Fs computed with betas in the incremental system. */
**ITM;                          /* (Ny-2)-by-(Nx/2).  Work array  */

mcomplex *fa,                   /* (Ny-2) vector.  Used for Fs computed when (Kx,Kz)=(0,0). */
*fb,                            /* (Ny-2) vector.  Used for Fs computed when (Kx,Kz)=(0,0). */
*tm;                            /* (Ny-2) vector.  Work array */

mcomplex *Ifa,                  /* (Ny-2) vector.  Used for Fs computed when (Kx,Kz)=(0,0) in the incremental system. */
*Ifb,                           /* (Ny-2) vector.  Used for Fs computed when (Kx,Kz)=(0,0) in the incremental system. */
*Itm;                           /* (Ny-2) vector.  Work array in the incremental system */

mcomplex **Uxb,                 /* Nz-by-(Nx/2) vector, Store dux at y=-1 as the boundary condition of incremental system */
**Uzb,                          /*  Nz-by-(Nx/2) vector, Store duz at y=-1 as the boundary condition of incremental system */
**Uxbt,                         /* Nz-by-(Nx/2) vector, Store dux at y=-1 as the boundary condition of incremental system */
**Uzbt;                         /*  Nz-by-(Nx/2) vector, Store duz at y=-1 as the boundary condition of incremental system */
mcomplex **IUxb, **IUzb;
mcomplex **IAUxb, **IAUzb;
mcomplex **GUxb, **GUzb;
mcomplex **GIUxb, **GIUzb;
mcomplex **HUxb, **HUzb;
mcomplex **HAUxb, **HAUzb;
mcomplex **AUxb, **AUzb;
mcomplex **grad, **hess, ishear;
fftw_plan pf1,                  /* complex to complex backward transform */
 pf2;                           /* complex to complex forward transform */
rfftwnd_plan pr1,               /* complex to real transform */
 pr2;                           /* real to complex transform */

fftw_plan Ipf1,                 /* complex to complex backward transform */
 Ipf2;                          /* complex to complex forward transform */


mcomplex ****MU, ****MIU;       /* variables used to store manufacture solutions */
fftw_complex ***CT;             /* 9-by-(3Nz/2)-by-(3*Nx/4+1).  For Fourier transforms. */
fftw_complex ***ICT;            /* 9-by-(3Nz/2)-by-(3*Nx/4+1).  For Fourier transforms in the incremental system. */

double tau, itau, maxtstep, atau, iatau;
double **uu, **us;

#endif

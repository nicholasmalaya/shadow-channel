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
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "minChnl.h"
#include "mvOps.h"

/* project when (Kx,Kz) = (0,0) */
void project0(int k, int n, func_force_t force)
{
    /* External Variables */
    extern int qpts, dimR, Nz, Nx;
    extern double dt, re, flux;
    extern mcomplex *fa, *fb, *tm;
    extern double **R, **Rp, **Rw, **Rs, **Rps, *Rp0;
    extern double **MZ;
    extern mcomplex **Uxb, **Uzb;
    extern mcomplex ****U, ****C;
    extern mcomplex *****MC;
    /* Local variables */
    int i, j;
    mcomplex tmp[dimR], tmp2[dimR], add[dimR];
    double flux_t;

    /* static double a[3] = {29./96., -3./40.,   1./6.};
       static double b[3] = {37./160., 5./24.,   1./6.};
       static double c[3] = { 8./15.,  5./12.,   3./4.};
       static double d[3] = { 0.,    -17./60.,  -5./12.};
       static double e[3] = { 0., 8./15., 2./3.}; */

    static double a[3] = { 1. / 3., -1. / 2, 1. / 3. };
    static double b[3] = { 1. / 6, 2. / 3, 0. };
    static double c[3] = { 1. / 2, 1. / 3., 1., };
    static double d[3] = { 0., -1. / 6, -2. / 3 };

    if (force != NULL) {
        force(n, k, 0, tmp, tmp2);
    }

    /* Create matrices for solving linear system. 
       Right hand side of system:  If this is the first step in the Runge-Kutta
       scheme, compute 
       [M0 + (1/RE)a[k]dt*D0]*C[0][ALPHA] + dt*c[k]Fu. 
       [M0 + (1/RE)a[k]dt*D0]*C[0][BETA] + dt*c[k]Fw. 
       Otherwise, compute C[0][ALPHA] + dt*c[k]Fu, 
       C[0][BETA] + dt*c[k]Fw, 
       Lhs matrix:  M0 - (1/RE)b[k]dt*D0                                */

    if (k == 0) {               /* first step */
        for (i = 0; i < dimR; ++i) {    /* MZ = M0 + (1/RE)a[k]dt*D0 */
            for (j = 0; j < T_RSDIAG; ++j) {
                MZ[i][j] = Rs[i][j] - re * a[0] * dt * Rps[i][j];
            }
        }

        /* compute [M0 + (1/RE)a[k]dt*D0]*C[0][ALPHA] for x=0 and store the
           result in tm.  Then transfer the result back to C[0][ALPHA].  
           Do the same for the BETAS */
        smMult0(MZ, C[0][ALPHA], tm, RSDIAG - 1, RSDIAG - 1, dimR);
        for (i = 0; i < dimR; ++i) {
            Re(C[0][ALPHA][i][0]) = Re(tm[i]);
            Im(C[0][ALPHA][i][0]) = Im(tm[i]);
        }

        smMult0(MZ, C[0][BETA], tm, RSDIAG - 1, RSDIAG - 1, dimR);
        for (i = 0; i < dimR; ++i) {
            Re(C[0][BETA][i][0]) = Re(tm[i]);
            Im(C[0][BETA][i][0]) = Im(tm[i]);
        }
    }

    /* MZ = M0 - (1/RE)b[k]dt*D0 */
    for (i = 0; i < dimR; ++i) {
        for (j = 0; j < T_RSDIAG; ++j) {
            MZ[i][j] = Rs[i][j] + re * b[k] * dt * Rps[i][j];
        }
    }

    memset(add, 0, dimR * sizeof(mcomplex));
    for (i = 0; i < dimR; ++i) {
        for (j = 0; j < qpts; ++j) {
            Re(add[i]) += Rw[i][j];
        }
    }
    /* Compute add  */
    bsolve_0(MZ, add, RSDIAG - 1, RSDIAG - 1, dimR);


    /* MZ = M0 - (1/RE)b[k]dt*D0 */
    for (i = 0; i < dimR; ++i) {
        for (j = 0; j < T_RSDIAG; ++j) {
            MZ[i][j] = Rs[i][j] + re * b[k] * dt * Rps[i][j];
        }
    }

    /* Finish computing the right hand size */
    /* array fa */
    memset(fa, 0, dimR * sizeof(mcomplex));
    for (i = 0; i < dimR; ++i) {
        for (j = 0; j < qpts; ++j) {
            Re(fa[i]) += Rw[i][j] * (Re(U[0][HXEL][j][0]));
        }
        if (force != NULL) {
            Re(fa[i]) += Re(tmp[i]);
            Im(fa[i]) += Im(tmp[i]);
        }
        Re(C[0][ALPHA][i][0]) += dt * c[k] * Re(fa[i]);
    }

    /* Compute a's */
    bsolve0(MZ, C[0][ALPHA], RSDIAG - 1, RSDIAG - 1, dimR);


    /* MZ = M0 - (1/RE)b[k]dt*D0 */
    for (i = 0; i < dimR; ++i) {
        for (j = 0; j < T_RSDIAG; ++j) {
            MZ[i][j] = Rs[i][j] + re * b[k] * dt * Rps[i][j];
        }
    }

    /* array fb */
    memset(fb, 0, dimR * sizeof(mcomplex));
    for (i = 0; i < dimR; ++i) {
        for (j = 0; j < qpts; ++j) {
            Re(fb[i]) += Rw[i][j] * Re(U[0][HZEL][j][0]);
        }

        if (force != NULL) {
            Re(fb[i]) += Re(tmp2[i]);
            Im(fb[i]) += Im(tmp2[i]);
        }

        Re(C[0][BETA][i][0]) += dt * c[k] * Re(fb[i]);
    }

    /* Compute b's */
    bsolve0(MZ, C[0][BETA], RSDIAG - 1, RSDIAG - 1, dimR);

    /* NOW COMPUTE U HATS */
    /* ux_hat = R*a */
    flux_t = 0;

    for (i = 0; i < dimR; ++i) {
        for (j = 0; j < qpts; ++j) {
            flux_t += Rw[i][j] * Re(C[0][ALPHA][i][0]);
        }
    }
    printf("flux_t=%f\n", flux_t);

    flux_t = flux - flux_t;

    /* for (j = 0; j < dimR; ++j)
       {
       Re(C[0][ALPHA][j][0])= Re(C[0][ALPHA][j][0])+flux_t*Re(add[j]);
       Im(C[0][ALPHA][j][0])= Im(C[0][ALPHA][j][0])+flux_t*Im(add[j]);
       } */

    for (i = 0; i < qpts; ++i) {
        Re(U[0][XEL][i][0]) = 0.0;
        Im(U[0][XEL][i][0]) = 0.0;
        for (j = 0; j < dimR; ++j) {
            Re(U[0][XEL][i][0]) += R[i][j] * Re(C[0][ALPHA][j][0]);
            Im(U[0][XEL][i][0]) += R[i][j] * Im(C[0][ALPHA][j][0]);
        }
    }

    /* dux_hat = Rp*a */
    for (i = 0; i < qpts; ++i) {
        Re(U[0][DXEL][i][0]) = 0.0;
        Im(U[0][DXEL][i][0]) = 0.0;
        for (j = 0; j < dimR; ++j) {
            Re(U[0][DXEL][i][0]) += Rp[i][j] * Re(C[0][ALPHA][j][0]);
            Im(U[0][DXEL][i][0]) += Rp[i][j] * Im(C[0][ALPHA][j][0]);
        }
    }

    /* uy_hat = 0 */
    for (i = 0; i < qpts; ++i) {
        Re(U[0][YEL][i][0]) = 0.0;
        Im(U[0][YEL][i][0]) = 0.0;
    }

    /* uz_hat = R*b */
    for (i = 0; i < qpts; ++i) {
        Re(U[0][ZEL][i][0]) = 0.0;
        Im(U[0][ZEL][i][0]) = 0.0;
        for (j = 0; j < dimR; ++j) {
            Re(U[0][ZEL][i][0]) += R[i][j] * Re(C[0][BETA][j][0]);
            Im(U[0][ZEL][i][0]) += R[i][j] * Im(C[0][BETA][j][0]);
        }
    }

    /* duz_hat = Rp*b */
    for (i = 0; i < qpts; ++i) {
        Re(U[0][DZEL][i][0]) = 0.0;
        Im(U[0][DZEL][i][0]) = 0.0;
        for (j = 0; j < dimR; ++j) {
            Re(U[0][DZEL][i][0]) += Rp[i][j] * Re(C[0][BETA][j][0]);
            Im(U[0][DZEL][i][0]) += Rp[i][j] * Im(C[0][BETA][j][0]);
        }
    }


    /* Uxb=Rp0*a, Uzb=Rp0*b, used as boundary condition for the incremental system */
    memset(Uxb[0], 0, Nz * (Nx / 2) * sizeof(mcomplex));
    memset(Uzb[0], 0, Nz * (Nx / 2) * sizeof(mcomplex));

    for (j = 0; j < dimR; ++j) {
        Re(Uxb[0][0]) += Rp0[j] * Re(C[0][ALPHA][j][0]);
        Im(Uxb[0][0]) += Rp0[j] * Im(C[0][ALPHA][j][0]);

        Re(Uzb[0][0]) += Rp0[j] * Re(C[0][BETA][j][0]);
        Im(Uzb[0][0]) += Rp0[j] * Im(C[0][BETA][j][0]);
    }

    /*
    for (i = 0; i < dimR; i++) {
        Re(MC[count][0][ALPHA][i][0]) = Re(C[0][ALPHA][i][0]);
        Im(MC[count][0][ALPHA][i][0]) = Im(C[0][ALPHA][i][0]);
        Re(MC[count][0][BETA][i][0]) = Re(C[0][BETA][i][0]);
        Im(MC[count][0][BETA][i][0]) = Im(C[0][BETA][i][0]);
    }
    */

    /* UPDATE RHS FOR NEXT TIME */
    if (k != 2) {               /* not last step */
        /* MZ = M0 + (1/RE)a[k+1]dt*D0 */
        for (i = 0; i < dimR; ++i) {
            for (j = 0; j < T_RSDIAG; ++j) {
                MZ[i][j] = Rs[i][j] - re * a[k + 1] * dt * Rps[i][j];
            }
        }

        /* compute [M0 + (1/RE)a[k+1]dt*D0]*C[0][ALPHA].  Then update
           C[0][ALPHA].  Do the same for the b's. */
        smMult0(MZ, C[0][ALPHA], tm, RSDIAG - 1, RSDIAG - 1, dimR);
        for (i = 0; i < dimR; ++i) {
            Re(C[0][ALPHA][i][0]) = Re(tm[i]) + dt * d[k + 1] * Re(fa[i]);
            Im(C[0][ALPHA][i][0]) = Im(tm[i]) + dt * d[k + 1] * Im(fa[i]);
        }

        /* now b's */
        smMult0(MZ, C[0][BETA], tm, RSDIAG - 1, RSDIAG - 1, dimR);
        for (i = 0; i < dimR; ++i) {
            Re(C[0][BETA][i][0]) = Re(tm[i]) + dt * d[k + 1] * Re(fb[i]);
            Im(C[0][BETA][i][0]) = Im(tm[i]) + dt * d[k + 1] * Im(fb[i]);
        }
    }
}                               /* end project zero (Kx,Kz) */

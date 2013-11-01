/* tangent forcing for Re variation when (Kx,Kz) != (0,0) */
#include <stdio.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include "minChnl.h"
#include "mvOps.h"


void tangent_forcing(int n, int k, int z,mcomplex * fb, mcomplex * fa)
{
    /* External Variables */
    extern int qpts, dimR, dimQ, Nx;
    extern double dt, re;
    extern double *Kx, *Kz, **K2;
    extern mcomplex **Fa, **Fb, **TM;
    extern double **Q, **Qp, **Qpp, **R, **Rp, **Qw, **Qpw, **Rw, **Qs,
        **Qps, **Qpps, **Rs, **Rps, *Rp0;
    extern double ***M;
    extern mcomplex ****U, ****C;
    extern mcomplex *****MC;
   /* Local variables */
    int i, j, x;
    double s;

    /* ALPHA Forcing */

    /* compute diffusion matrix M = -(1/RE^2) * Dv for ALPHAS  */
    if (k == 0) {               /* first step */
        memset(M[0][0], 0, dimR * 9 * (Nx / 2) * sizeof(double));
        for (i = 0; i < dimQ; ++i) {    /* M = -(1/RE^2)*Dv */
            for (j = 0; j < T_QSDIAG; ++j) {
                for (x = x0; x < Nx / 2; ++x) {
                    s = K2[z][x] * K2[z][x];
                    M[i][j][x] = - re * re * (s * Qs[i][j] + 2. * K2[z][x] * Qps[i][j] + Qpps[i][j]);
                }
            }
        }
        /* compute M*C[z][ALPHA] and store the result in fa. */
        smMult(M, C[z][ALPHA], fa, QSDIAG - 1, QSDIAG - 1, dimQ, Nx / 2, x0);
        
    }

    /* diffusion matrix M = -(1/RE^2)*Dv for ALPHAS */
    memset(M[0][0], 0, dimR * 9 * (Nx / 2) * sizeof(double));
    for (i = 0; i < dimQ; ++i) {
        for (j = 0; j < T_QSDIAG; ++j) {
            for (x = x0; x < Nx / 2; ++x) {
                s = K2[z][x] * K2[z][x];
                M[i][j][x] = - re * re * (s * Qs[i][j] + 2. * K2[z][x] * Qps[i][j] + Qpps[i][j]);
            }
        }
    }

    /* Compute alpha forcing  */
    smMult(M, C[z][ALPHA], fa, QSDIAG - 1, QSDIAG - 1, dimQ, Nx / 2, x0);
    

    /* BETA Forcing */

    if (k == 0) {               /* first step */
        /* M = -(1/RE)^2 * Dg */
        memset(M[0][0], 0, dimR * 9 * (Nx / 2) * sizeof(double));
        for (i = 0; i < dimR; ++i) {
            for (j = 0; j < T_RSDIAG; ++j) {
                for (x = x0; x < Nx / 2; ++x) {
                    M[i][j][x] = - re * re * (Rps[i][j] + K2[z][x] * Rs[i][j]);
                }
            }
        }

        /* compute M*C[z][BETA] and store the result in fb.*/
        smMult(M, C[z][BETA], fb, RSDIAG - 1, RSDIAG - 1, dimR, Nx / 2, x0);
    }

    /* M = -(1/RE^2)*Dg */
    memset(M[0][0], 0, dimR * 9 * (Nx / 2) * sizeof(double));
    for (i = 0; i < dimR; ++i) {        /* M = - (1/RE)^2*Dg */
        for (j = 0; j < T_RSDIAG; ++j) {
            for (x = x0; x < Nx / 2; ++x) {
                M[i][j][x] = - re * re * (Rps[i][j] + K2[z][x] * Rs[i][j]);
            }
        }
    }

    /* compute M*C[z][BETA] and store the result in fb.*/
    smMult(M, C[z][BETA], fb, RSDIAG - 1, RSDIAG - 1, dimR, Nx / 2, x0);


}


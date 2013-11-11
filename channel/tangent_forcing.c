/* tangent forcing for Re variation when (Kx,Kz) != (0,0) */
#include <stdio.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include "minChnl.h"
#include "mvOps.h"


void tangent_forcing(int count, int k, int z, mcomplex ** f_a, mcomplex ** f_b)
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
    int i, j, x, x0;
    double s;

    /* first x */
    x0 = (z == 0) ? 1 : 0;
    /* ALPHA Forcing */

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
    smMult(M, MC[count][z][ALPHA], f_a, QSDIAG - 1, QSDIAG - 1, dimQ, Nx / 2, x0);

    double f_a_norm = 0.0;
    for (i = 0; i < dimR; ++i) {
        for (x = x0; x < Nx / 2; ++x) {
            f_a_norm += pow(MAGNITUDE(f_a[i][x]), 2);
        }
    }
    printf("z = %d, f_a_norm = %e\n", z, sqrt(f_a_norm));
    

    /* BETA Forcing */
    /* M = -(1/RE^2)*Dg */
    memset(M[0][0], 0, dimR * 9 * (Nx / 2) * sizeof(double));
    for (i = 0; i < dimR; ++i) {        /* M = - (1/RE)^2*Dg */
        for (j = 0; j < T_RSDIAG; ++j) {
            for (x = x0; x < Nx / 2; ++x) {
                M[i][j][x] = re * re * (Rps[i][j] + K2[z][x] * Rs[i][j]);
            }
        }
    }

    /* compute M*C[z][BETA] and store the result in f_b.*/
    smMult(M, MC[count][z][BETA], f_b, RSDIAG - 1, RSDIAG - 1, dimR, Nx / 2, x0);

    double f_b_norm = 0.0;
    for (i = 0; i < dimR; ++i) {
        for (x = x0; x < Nx / 2; ++x) {
            f_b_norm += pow(MAGNITUDE(f_b[i][x]), 2);
        }
    }
    printf("z = %d, f_b_norm = %e\n", z, sqrt(f_b_norm));

}


/**************************************************************************
This function is similar to initAlphabeta2.c. It uses the stored Fourier
coefficients a, b, alpha,beta from the adjoint state system to compute
the adjoint veolcity Aux, Auy, Auz, dAux, dAuz.
**************************************************************************/
#include <stdio.h>
#include <string.h>
#include "minChnl.h"
#include "mvOps.h"
#include "arrays.h"

void adj_initAlphaBeta2(void)
{

    void adj_initAlphaBeta_2(int z, int x0);
    extern int qpts, dimR, dimQ, Nx, Nz;
    extern mcomplex ****AU, ****AC;
    extern double **R, **Rp, *Rp0;
    int i, j, z;
    extern double *Uadd;
    /* Aux_hat = R*a */
    /* ux_hat = R*a */
    for (i = 0; i < qpts; ++i) {
        Re(AU[0][XEL][i][0]) = 0.0;
        Im(AU[0][XEL][i][0]) = 0.0;
        for (j = 0; j < dimR; ++j) {
            Re(AU[0][XEL][i][0]) += R[i][j] * Re(AC[0][ALPHA][j][0]);
            Im(AU[0][XEL][i][0]) += R[i][j] * Im(AC[0][ALPHA][j][0]);
        }

    }

    /* dux_hat = Rp*a */
    for (i = 0; i < qpts; ++i) {
        Re(AU[0][DXEL][i][0]) = 0.0;
        Im(AU[0][DXEL][i][0]) = 0.0;
        for (j = 0; j < dimR; ++j) {
            Re(AU[0][DXEL][i][0]) += Rp[i][j] * Re(AC[0][ALPHA][j][0]);
            Im(AU[0][DXEL][i][0]) += Rp[i][j] * Im(AC[0][ALPHA][j][0]);
        }
    }
    /* uy_hat = 0 */
    for (i = 0; i < qpts; ++i) {
        Re(AU[0][YEL][i][0]) = 0.0;
        Im(AU[0][YEL][i][0]) = 0.0;
    }

    /* uz_hat = R*b */
    for (i = 0; i < qpts; ++i) {
        Re(AU[0][ZEL][i][0]) = 0.0;
        Im(AU[0][ZEL][i][0]) = 0.0;
        for (j = 0; j < dimR; ++j) {
            Re(AU[0][ZEL][i][0]) += R[i][j] * Re(AC[0][BETA][j][0]);
            Im(AU[0][ZEL][i][0]) += R[i][j] * Im(AC[0][BETA][j][0]);
        }
    }

    /* duz_hat = Rp*b */
    for (i = 0; i < qpts; ++i) {
        Re(AU[0][DZEL][i][0]) = 0.0;
        Im(AU[0][DZEL][i][0]) = 0.0;
        for (j = 0; j < dimR; ++j) {
            Re(AU[0][DZEL][i][0]) += Rp[i][j] * Re(AC[0][BETA][j][0]);
            Im(AU[0][DZEL][i][0]) += Rp[i][j] * Im(AC[0][BETA][j][0]);
        }
    }
  
    adj_initAlphaBeta_2(0, 1);
    for (z = 1; z < Nz; ++z) {
        adj_initAlphaBeta_2(z, 0);
    }
}


void adj_initAlphaBeta_2(int z, int x0)
{
    extern int qpts, dimR, dimQ, Nx, Nz;
    extern mcomplex ****AU, ****AC;
    extern double **Q, **Qp, **Qpp, **R, **Rp, *Rp0;
    extern double *Kx, *Kz, **K2;
    int i, j, x;
    double t[2];
    extern double *Uadd, *Vadd, *Vpadd;
    int sign;
    /* NOW COMPUTE AU HATS */
    /* v = uy_hat. */
    for (i = 0; i < qpts; ++i) {
        memset(&AU[z][YEL][i][x0], 0, (Nx / 2 - x0) * sizeof(mcomplex));
    }
    for (i = 0; i < qpts; ++i) {
        for (j = 0; j < dimQ; ++j) {
            for (x = x0; x < Nx / 2; ++x) {
                Re(AU[z][YEL][i][x]) += Q[i][j] * Re(AC[z][ALPHA][j][x]);
                Im(AU[z][YEL][i][x]) += Q[i][j] * Im(AC[z][ALPHA][j][x]);
            }
        }
    }

    /* f = -dv/dy and store temporarily in XEL position of array U. */
    for (i = 0; i < qpts; ++i) {
        memset(&AU[z][XEL][i][x0], 0, (Nx / 2 - x0) * sizeof(mcomplex));
    }
    for (i = 0; i < qpts; ++i) {
        for (j = 0; j < dimQ; ++j) {
            for (x = x0; x < Nx / 2; ++x) {
                Re(AU[z][XEL][i][x]) -= Qp[i][j] * Re(AC[z][ALPHA][j][x]);
                Im(AU[z][XEL][i][x]) -= Qp[i][j] * Im(AC[z][ALPHA][j][x]);
            }
        }
    }

    /* sum(Q''alpha) and store in DXEL position. */
    for (i = 0; i < qpts; ++i) {
        memset(&AU[z][DXEL][i][x0], 0, (Nx / 2 - x0) * sizeof(mcomplex));
    }
    for (i = 0; i < qpts; ++i) {
        for (j = 0; j < dimQ; ++j) {
            for (x = x0; x < Nx / 2; ++x) {
                Re(AU[z][DXEL][i][x]) +=
                    Qpp[i][j] * Re(AC[z][ALPHA][j][x]);
                Im(AU[z][DXEL][i][x]) +=
                    Qpp[i][j] * Im(AC[z][ALPHA][j][x]);
            }
        }
    }

    /* Compute g = sum(beta*R) and store temporarily in ZEL position of 
       array U. */
    for (i = 0; i < qpts; ++i) {
        memset(&AU[z][ZEL][i][x0], 0, (Nx / 2 - x0) * sizeof(mcomplex));
    }
    for (i = 0; i < qpts; ++i) {
        for (j = 0; j < dimR; ++j) {
            for (x = x0; x < Nx / 2; ++x) {
                Re(AU[z][ZEL][i][x]) += R[i][j] * Re(AC[z][BETA][j][x]);
                Im(AU[z][ZEL][i][x]) += R[i][j] * Im(AC[z][BETA][j][x]);
            }
        }
    }

    /* Compute sum(beta*R') and store temporarily in DZEL position of
       array U. */
    for (i = 0; i < qpts; ++i) {
        memset(&AU[z][DZEL][i][x0], 0, (Nx / 2 - x0) * sizeof(mcomplex));
    }
    for (i = 0; i < qpts; ++i) {
        for (j = 0; j < dimR; ++j) {
            for (x = x0; x < Nx / 2; ++x) {
                Re(AU[z][DZEL][i][x]) += Rp[i][j] * Re(AC[z][BETA][j][x]);
                Im(AU[z][DZEL][i][x]) += Rp[i][j] * Im(AC[z][BETA][j][x]);
            }
        }
    }

    /* now compute ux hat, uz hat */
    for (i = 0; i < qpts; ++i) {
        for (x = x0; x < Nx / 2; ++x) {
            t[0] = Re(AU[z][XEL][i][x]);        /* real part of f */
            t[1] = Re(AU[z][ZEL][i][x]);        /* real part of g */

            Re(AU[z][XEL][i][x]) = (Kx[x] * Im(AU[z][XEL][i][x]) +
                                    Kz[z] * Im(AU[z][ZEL][i][x])) /
                K2[z][x];
            Re(AU[z][ZEL][i][x]) =
                (-Kx[x] * Im(AU[z][ZEL][i][x]) +
                 Kz[z] * Im(AU[z][XEL][i][x])) / K2[z][x];
            Im(AU[z][XEL][i][x]) =
                -(Kx[x] * t[0] + Kz[z] * t[1]) / K2[z][x];
            Im(AU[z][ZEL][i][x]) =
                (Kx[x] * t[1] - Kz[z] * t[0]) / K2[z][x];
        }
    }

    /* dux hat, duz hat */
    for (i = 0; i < qpts; ++i) {
        for (x = x0; x < Nx / 2; ++x) {
            t[0] = Re(AU[z][DXEL][i][x]);       /* real part of Q''alpha */
            t[1] = Re(AU[z][DZEL][i][x]);       /* real part of R'beta */

            Re(AU[z][DXEL][i][x]) = (-Kx[x] * Im(AU[z][DXEL][i][x]) +
                                     Kz[z] * Im(AU[z][DZEL][i][x])) /
                K2[z][x];
            Re(AU[z][DZEL][i][x]) =
                -(Kx[x] * Im(AU[z][DZEL][i][x]) +
                  Kz[z] * Im(AU[z][DXEL][i][x])) / K2[z][x];
            Im(AU[z][DXEL][i][x]) =
                (Kx[x] * t[0] - Kz[z] * t[1]) / K2[z][x];
            Im(AU[z][DZEL][i][x]) =
                (Kx[x] * t[1] + Kz[z] * t[0]) / K2[z][x];
        }
    }

}

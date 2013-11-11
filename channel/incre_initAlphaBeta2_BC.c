/**************************************************************************
This function is similar to initAlphabeta2.c. It uses the stored Fourier
coefficients a, b, alpha,beta from the incremental state system to compute
the veolcity Iux, Iuy, Iuz, dIux, dIuz.
**************************************************************************/
#include <stdio.h>
#include <string.h>
#include "minChnl.h"
#include "mvOps.h"
#include "arrays.h"

void incre_initAlphaBeta2(void)
{

    void incre_initAlphaBeta_2(int z, int x0);
    extern int qpts, dimR, dimQ, Nx, Nz;
    extern mcomplex **Uxb, **Uzb;
    extern mcomplex **IUxb, **IUzb;
    extern mcomplex ****IU, ****IC;
    extern double **R, **Rp, *Rp0;
    int i, j, z;
    extern double *Uadd;
    /* Iux_hat = R*a+c3/2*(1-y) */
    for (i = 0; i < qpts; ++i) {
        Re(IU[0][XEL][i][0]) = 0.0;
        Im(IU[0][XEL][i][0]) = 0.0;
        for (j = 0; j < dimR; ++j) {
            Re(IU[0][XEL][i][0]) += R[i][j] * Re(IC[0][ALPHA][j][0]);
            Im(IU[0][XEL][i][0]) += R[i][j] * Im(IC[0][ALPHA][j][0]);
        }
        Re(IU[0][XEL][i][0]) += Re(Uxb[0][0]) * Uadd[i] / 2.;
        Im(IU[0][XEL][i][0]) += Im(Uxb[0][0]) * Uadd[i] / 2.;
    }

    /* dIux_hat = Rp*a-c3/2 */
    for (i = 0; i < qpts; ++i) {
        Re(IU[0][DXEL][i][0]) = 0.0;
        Im(IU[0][DXEL][i][0]) = 0.0;
        for (j = 0; j < dimR; ++j) {
            Re(IU[0][DXEL][i][0]) += Rp[i][j] * Re(IC[0][ALPHA][j][0]);
            Im(IU[0][DXEL][i][0]) += Rp[i][j] * Im(IC[0][ALPHA][j][0]);
        }
        Re(IU[0][DXEL][i][0]) += -Re(Uxb[0][0]) / 2.;
        Im(IU[0][DXEL][i][0]) += -Im(Uxb[0][0]) / 2.;
    }

    /* Iuy_hat = 0 */
    for (i = 0; i < qpts; ++i) {
        Re(IU[0][YEL][i][0]) = 0.0;
        Im(IU[0][YEL][i][0]) = 0.0;
    }

    /* Iuz_hat = R*b +c4/2*(1-y) */
    for (i = 0; i < qpts; ++i) {
        Re(IU[0][ZEL][i][0]) = 0.0;
        Im(IU[0][ZEL][i][0]) = 0.0;
        for (j = 0; j < dimR; ++j) {
            Re(IU[0][ZEL][i][0]) += R[i][j] * Re(IC[0][BETA][j][0]);
            Im(IU[0][ZEL][i][0]) += R[i][j] * Im(IC[0][BETA][j][0]);
        }
        Re(IU[0][ZEL][i][0]) += Re(Uzb[0][0]) * Uadd[i] / 2.;
        Im(IU[0][ZEL][i][0]) += Im(Uzb[0][0]) * Uadd[i] / 2.;
    }

    /* duz_hat = Rp*b-c4/2 */
    for (i = 0; i < qpts; ++i) {
        Re(IU[0][DZEL][i][0]) = 0.0;
        Im(IU[0][DZEL][i][0]) = 0.0;
        for (j = 0; j < dimR; ++j) {
            Re(IU[0][DZEL][i][0]) += Rp[i][j] * Re(IC[0][BETA][j][0]);
            Im(IU[0][DZEL][i][0]) += Rp[i][j] * Im(IC[0][BETA][j][0]);
        }
        Re(IU[0][DZEL][i][0]) += -Re(Uzb[0][0]) / 2.;
        Im(IU[0][DZEL][i][0]) += -Im(Uzb[0][0]) / 2.;
    }

    memset(IUxb[0], 0, Nz * (Nx / 2) * sizeof(fftw_complex));
    memset(IUzb[0], 0, Nz * (Nx / 2) * sizeof(fftw_complex));
    for (j = 0; j < dimR; ++j) {
        Re(IUxb[0][0]) += Rp0[j] * Re(IC[0][ALPHA][j][0]);
        Im(IUxb[0][0]) += Rp0[j] * Im(IC[0][ALPHA][j][0]);

        Re(IUzb[0][0]) += Rp0[j] * Re(IC[0][BETA][j][0]);
        Im(IUzb[0][0]) += Rp0[j] * Im(IC[0][BETA][j][0]);
    }

    Re(IUxb[0][0]) -= Re(Uxb[0][0]) / 2.;
    Im(IUxb[0][0]) -= Im(Uxb[0][0]) / 2.;
    Re(IUzb[0][0]) -= Re(Uzb[0][0]) / 2.;
    Im(IUzb[0][0]) -= Im(Uzb[0][0]) / 2.;

    incre_initAlphaBeta_2(0, 1);
    for (z = 1; z < Nz; ++z) {
        incre_initAlphaBeta_2(z, 0);
    }
}


void incre_initAlphaBeta_2(int z, int x0)
{
    extern int qpts, dimR, dimQ, Nx, Nz;
    extern mcomplex **Uxb, **Uzb;
    extern mcomplex **IUxb, **IUzb;
    extern mcomplex ****IU, ****IC;
    extern double **Q, **Qp, **Qpp, **R, **Rp, *Rp0;
    extern double *Kx, *Kz, **K2;
    int i, j, x;
    double t[2];
    extern double *Uadd, *Vadd, *Vpadd;
    int sign;
    /* NOW COMPUTE IU HATS */
    /* v = uy_hat+c2/4*(1-y)^2*(1+y). */
    for (i = 0; i < qpts; ++i) {
        memset(&IU[z][YEL][i][x0], 0, (Nx / 2 - x0) * sizeof(mcomplex));
    }
    for (i = 0; i < qpts; ++i) {
        for (x = x0; x < Nx / 2; ++x) {
            for (j = 0; j < dimQ; ++j) {
                Re(IU[z][YEL][i][x]) += Q[i][j] * Re(IC[z][ALPHA][j][x]);
                Im(IU[z][YEL][i][x]) += Q[i][j] * Im(IC[z][ALPHA][j][x]);
            }
            Re(IU[z][YEL][i][x]) += Re(Uzb[z][x]) * Vadd[i] / 4.;
            Im(IU[z][YEL][i][x]) += Im(Uzb[z][x]) * Vadd[i] / 4.;
        }
    }


    /* f = -dv/dy and store temporarily in XEL position of array U. */
    /*f=-dv/dy-c2/4*(3y^2-2y-1) */
    for (i = 0; i < qpts; ++i) {
        memset(&IU[z][XEL][i][x0], 0, (Nx / 2 - x0) * sizeof(mcomplex));
    }
    for (i = 0; i < qpts; ++i) {
        for (x = x0; x < Nx / 2; ++x) {
            for (j = 0; j < dimQ; ++j) {
                Re(IU[z][XEL][i][x]) -= Qp[i][j] * Re(IC[z][ALPHA][j][x]);
                Im(IU[z][XEL][i][x]) -= Qp[i][j] * Im(IC[z][ALPHA][j][x]);
            }
            Re(IU[z][XEL][i][x]) -= Re(Uzb[z][x]) * Vpadd[i] / 4.;
            Im(IU[z][XEL][i][x]) -= Im(Uzb[z][x]) * Vpadd[i] / 4.;
        }
    }

    /* sum(Q''alpha) and store in DXEL position. */
    /*df/dy=d^2 v/dy^2+c2/4*(6y-2)=d^2 v/dy^2+c2/2*(-3+3y+2) */
    for (i = 0; i < qpts; ++i) {
        memset(&IU[z][DXEL][i][x0], 0, (Nx / 2 - x0) * sizeof(mcomplex));
    }
    for (i = 0; i < qpts; ++i) {
        for (x = x0; x < Nx / 2; ++x) {
            for (j = 0; j < dimQ; ++j) {
                Re(IU[z][DXEL][i][x]) +=
                    Qpp[i][j] * Re(IC[z][ALPHA][j][x]);
                Im(IU[z][DXEL][i][x]) +=
                    Qpp[i][j] * Im(IC[z][ALPHA][j][x]);
            }
            Re(IU[z][DXEL][i][x]) +=
                Re(Uzb[z][x]) * (-3 * Uadd[i] + 2) / 2.;
            Im(IU[z][DXEL][i][x]) +=
                Im(Uzb[z][x]) * (-3 * Uadd[i] + 2) / 2.;
        }
    }

    /* Compute g = sum(beta*R) and store temporarily in ZEL position of 
       array U. */
    for (i = 0; i < qpts; ++i) {
        memset(&IU[z][ZEL][i][x0], 0, (Nx / 2 - x0) * sizeof(mcomplex));
    }
    for (i = 0; i < qpts; ++i) {
        for (x = x0; x < Nx / 2; ++x) {
            for (j = 0; j < dimR; ++j) {
                Re(IU[z][ZEL][i][x]) += R[i][j] * Re(IC[z][BETA][j][x]);
                Im(IU[z][ZEL][i][x]) += R[i][j] * Im(IC[z][BETA][j][x]);
            }

            Re(IU[z][ZEL][i][x]) += Re(Uxb[z][x]) * Uadd[i] / 2.;
            Im(IU[z][ZEL][i][x]) += Im(Uxb[z][x]) * Uadd[i] / 2.;
        }
    }


    /* Compute sum(beta*R') and store temporarily in DZEL position of
       array IU. */
    for (i = 0; i < qpts; ++i) {
        memset(&IU[z][DZEL][i][x0], 0, (Nx / 2 - x0) * sizeof(mcomplex));
    }
    for (i = 0; i < qpts; ++i) {
        for (x = x0; x < Nx / 2; ++x) {
            for (j = 0; j < dimR; ++j) {
                Re(IU[z][DZEL][i][x]) += Rp[i][j] * Re(IC[z][BETA][j][x]);
                Im(IU[z][DZEL][i][x]) += Rp[i][j] * Im(IC[z][BETA][j][x]);
            }
            Re(IU[z][DZEL][i][x]) += -Re(Uxb[z][x]) / 2.;
            Im(IU[z][DZEL][i][x]) += -Im(Uxb[z][x]) / 2.;
        }
    }

    /* now compute Iux hat, Iuz hat */
    for (i = 0; i < qpts; ++i) {
        for (x = x0; x < Nx / 2; ++x) {
            t[0] = Re(IU[z][XEL][i][x]);        /* real part of f */
            t[1] = Re(IU[z][ZEL][i][x]);        /* real part of g */

            Re(IU[z][XEL][i][x]) = (Kx[x] * Im(IU[z][XEL][i][x]) +
                                    Kz[z] * Im(IU[z][ZEL][i][x])) /
                K2[z][x];
            Re(IU[z][ZEL][i][x]) =
                (-Kx[x] * Im(IU[z][ZEL][i][x]) +
                 Kz[z] * Im(IU[z][XEL][i][x])) / K2[z][x];
            Im(IU[z][XEL][i][x]) =
                -(Kx[x] * t[0] + Kz[z] * t[1]) / K2[z][x];
            Im(IU[z][ZEL][i][x]) =
                (Kx[x] * t[1] - Kz[z] * t[0]) / K2[z][x];
        }
    }

    /* dux hat, duz hat */
    for (i = 0; i < qpts; ++i) {
        for (x = x0; x < Nx / 2; ++x) {
            t[0] = Re(IU[z][DXEL][i][x]);       /* real part of Q''alpha */
            t[1] = Re(IU[z][DZEL][i][x]);       /* real part of R'beta */

            Re(IU[z][DXEL][i][x]) = (-Kx[x] * Im(IU[z][DXEL][i][x]) +
                                     Kz[z] * Im(IU[z][DZEL][i][x])) /
                K2[z][x];
            Re(IU[z][DZEL][i][x]) =
                -(Kx[x] * Im(IU[z][DZEL][i][x]) +
                  Kz[z] * Im(IU[z][DXEL][i][x])) / K2[z][x];
            Im(IU[z][DXEL][i][x]) =
                (Kx[x] * t[0] - Kz[z] * t[1]) / K2[z][x];
            Im(IU[z][DZEL][i][x]) =
                (Kx[x] * t[1] + Kz[z] * t[0]) / K2[z][x];
        }
    }

    for (x = x0; x < Nx / 2; ++x) {
        for (j = 0; j < dimR; ++j) {
            Re(IUxb[z][x]) +=
                (Rp0[j] * Im(IC[z][BETA][j][x])) * Kz[z] / K2[z][x];
            Im(IUxb[z][x]) +=
                (-Rp0[j] * Re(IC[z][BETA][j][x])) * Kz[z] / K2[z][x];

            Re(IUzb[z][x]) +=
                -(Rp0[j] * Im(IC[z][BETA][j][x])) * Kx[x] / K2[z][x];
            Im(IUzb[z][x]) +=
                (Rp0[j] * Re(IC[z][BETA][j][x])) * Kx[x] / K2[z][x];

        }
        Re(IUxb[z][x]) += -Im(Uxb[z][x]) / 2. * Kz[z] / K2[z][x];
        Im(IUxb[z][x]) += Re(Uxb[z][x]) / 2. * Kz[z] / K2[z][x];
        Re(IUzb[z][x]) += Im(Uxb[z][x]) / 2. * Kx[x] / K2[z][x];
        Im(IUzb[z][x]) += -Re(Uxb[z][x]) / 2. * Kx[x] / K2[z][x];
    }

    for (x = x0; x < Nx / 2; ++x) {
        sign = -1;
        for (j = 0; j < dimQ; ++j) {
            sign = sign * (-1);
            Re(IUxb[z][x]) -=
                (sign * 8. * Im(IC[z][ALPHA][j][x])) * Kx[x] / K2[z][x];
            Im(IUxb[z][x]) -=
                (-sign * 8. * Re(IC[z][ALPHA][j][x])) * 2 * Kx[x] /
                K2[z][x];

            Re(IUzb[z][x]) -=
                (sign * 8. * Im(IC[z][ALPHA][j][x])) * Kz[z] / K2[z][x];
            Im(IUzb[z][x]) -=
                (-sign * 8. * Re(IC[z][ALPHA][j][x])) * 2 * Kz[z] /
                K2[z][x];
        }
        Re(IUxb[z][x]) -= -Im(Uzb[z][x]) * 2. * Kx[x] / K2[z][x];
        Im(IUxb[z][x]) -= Re(Uzb[z][x]) * 2. * Kx[x] / K2[z][x];
        Re(IUzb[z][x]) -= -Im(Uzb[z][x]) * 2. * Kz[z] / K2[z][x];
        Im(IUzb[z][x]) -= Re(Uzb[z][x]) * 2 * Kz[z] / K2[z][x];
    }

}

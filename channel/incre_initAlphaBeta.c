
#include <stdio.h>
#include <string.h>
#include "minChnl.h"
#include "mvOps.h"
#include "arrays.h"

/* Compute initial coefficients in expansion, given the values for ux_hat,
   uy_hat, uz_hat */
void incre_initAlphaBeta(void)
{
    /* function prototypes */
    void incre_Alpha0(void);
    void incre_Beta0(void);
    void incre_Alpha(int z, int x0, int xdim);
    void incre_Beta(int z, int x0, int xdim);

    /* variables */
    extern int Nz, Nx;
    int z;

    /* (Kx,Kz) = (0,0) */
    incre_Alpha0();
    incre_Beta0();

    /* Kz = 0; Kx != 0 */
    incre_Alpha(0, 1, Nx / 2);
    incre_Beta(0, 1, Nx / 2);

    /* Kz != 0; Kx = 0,... */
    for (z = 1; z < Nz; ++z) {
        incre_Alpha(z, 0, Nx / 2);
        incre_Beta(z, 0, Nx / 2);
    }


}                               /* end incre_initAlphaBeta */


void incre_Alpha0(void)
{
    extern int qpts, dimR;
    extern double **Rw, **Rs, **MZ, *Qy;
    extern mcomplex ****IU, ****IC;
    extern mcomplex **Uxb;
    int i, j;

    for (i = 0; i < dimR; ++i) {
        Re(IC[0][ALPHA][i][0]) = 0.0;
        Im(IC[0][ALPHA][i][0]) = 0.0;
        for (j = 0; j < qpts; ++j) {
            Re(IC[0][ALPHA][i][0]) +=
                Rw[i][j] * (Re(IU[0][XEL][j][0])); // -
                          //  Re(Uxb[0][0]) * (1 - Qy[j]) / 2.);
            Im(IC[0][ALPHA][i][0]) +=
                Rw[i][j] * (Im(IU[0][XEL][j][0])); // -
                          //  Im(Uxb[0][0]) * (1 - Qy[j]) / 2.);
        }
    }

    memcpy(MZ[0], Rs[0], dimR * T_RSDIAG * sizeof(double));
    bsolve0(MZ, IC[0][ALPHA], RSDIAG - 1, RSDIAG - 1, dimR);

}                               /* end incre_Alpha0() */


void incre_Beta0(void)
{
    extern int qpts, dimR;
    extern double **Rw, **Rs, **MZ, *Qy;
    extern mcomplex ****IU, ****IC;
    extern mcomplex **Uxb, **Uzb;
    int i, j;

    for (i = 0; i < dimR; ++i) {
        Re(IC[0][BETA][i][0]) = 0.0;
        Im(IC[0][BETA][i][0]) = 0.0;
        for (j = 0; j < qpts; ++j) {
            Re(IC[0][BETA][i][0]) +=
                Rw[i][j] * (Re(IU[0][ZEL][j][0])); // -
                           // Re(Uzb[0][0]) * (1 - Qy[j]) / 2.);
            Im(IC[0][BETA][i][0]) +=
                Rw[i][j] * (Im(IU[0][ZEL][j][0])); // -
                           // Im(Uzb[0][0]) * (1 - Qy[j]) / 2.);
        }
    }

    memcpy(MZ[0], Rs[0], dimR * T_RSDIAG * sizeof(double));
    bsolve0(MZ, IC[0][BETA], RSDIAG - 1, RSDIAG - 1, dimR);
}                               /* end incre_Beta0() */


void incre_Alpha(int z, int x0, int xdim)
{
    void bsolveI(double **A, mcomplex ** b, int dl, int du, int N,
                 int xdim, int x0);

    extern int qpts, dimQ;
    extern double **Qw, **Qs, *Qy;
    extern mcomplex ****IU, ****IC;
    extern mcomplex **Uxb, **Uzb;
    int i, j, x;
    double **M;

    if ((M = dMatrix(dimQ, T_QSDIAG)) == NULL) {
        printf("incre-Alpha: No memory for M array.\n");
        return;
    }

    for (i = 0; i < dimQ; ++i) {
        for (x = x0; x < xdim; ++x) {
            Re(IC[z][ALPHA][i][x]) = 0.0;
            Im(IC[z][ALPHA][i][x]) = 0.0;
        }
    }
    for (i = 0; i < dimQ; ++i) {
        for (j = 0; j < qpts; ++j) {
            for (x = x0; x < xdim; ++x) {
                Re(IC[z][ALPHA][i][x]) +=
                    Qw[i][j] * (Re(IU[z][YEL][j][x]));/* -
                                Re(Uzb[z][x]) / 4. * (1 - Qy[j]) * (1 -
                                                                    Qy[j] *
                                                                    Qy
                                                                    [j])); */
                Im(IC[z][ALPHA][i][x]) +=
                    Qw[i][j] * (Im(IU[z][YEL][j][x]));/* -
                                Im(Uzb[z][x]) / 4. * (1 - Qy[j]) * (1 -
                                                                    Qy[j] *
                                                                    Qy
                                                                    [j])); */
            }
        }
    }

    memcpy(M[0], Qs[0], dimQ * T_QSDIAG * sizeof(double));
    bsolveI(M, IC[z][ALPHA], QSDIAG - 1, QSDIAG - 1, dimQ, xdim, x0);
    freedMatrix(M);
}                               /* end incre_Alpha() */


void incre_Beta(int z, int x0, int xdim)
{
    void bsolveI(double **A, mcomplex ** b, int dl, int du, int N,
                 int xdim, int x0);

    extern int qpts, dimR;
    extern double *Kx, *Kz, *Qy;
    extern double **Rw, **Rs, **MZ;
    extern mcomplex ****IU, ****IC;
    extern mcomplex **Uxb, **Uzb;
    int i, j, x;

    for (i = 0; i < dimR; ++i) {
        for (x = x0; x < xdim; ++x) {
            Re(IC[z][BETA][i][x]) = 0.0;
            Im(IC[z][BETA][i][x]) = 0.0;
        }
    }
    for (i = 0; i < dimR; ++i) {
        for (j = 0; j < qpts; ++j) {
            for (x = x0; x < xdim; ++x) {
                Re(IC[z][BETA][i][x]) += Rw[i][j] *
                    ((Kx[x] * Im(IU[z][ZEL][j][x]) -
                      Kz[z] * Im(IU[z][XEL][j][x]))); /* -
                     Re(Uxb[z][x]) / 2. * (1 - Qy[j]));*/
                Im(IC[z][BETA][i][x]) +=
                    Rw[i][j] *
                    ((Kz[z] * Re(IU[z][XEL][j][x]) -
                      Kx[x] * Re(IU[z][ZEL][j][x]))); /* -
                     Im(Uxb[z][x]) / 2. * (1 - Qy[j]));*/
            }
        }
    }

    memcpy(MZ[0], Rs[0], dimR * T_RSDIAG * sizeof(double));
    bsolveI(MZ, IC[z][BETA], RSDIAG - 1, RSDIAG - 1, dimR, xdim, x0);

}                               /* end Beta() */

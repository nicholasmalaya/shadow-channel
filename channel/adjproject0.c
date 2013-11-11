#include <math.h>
#include <stdio.h>
#include <string.h>
#include "minChnl.h"
#include "mvOps.h"

/* project when (Kx,Kz) = (0,0) */
void adjproject0(int count, int k, func_force_t force)
{
    /* External Variables */
    extern int qpts, dimR, Nx, Nz;
    extern double dt, re;
    extern mcomplex *fa, *fb, *tm;
    extern double **R, **Rp, **Rw, **Rpw, **Rs, **Rps, *Rp0, *Rpp0;
    extern double **MZ;
    extern mcomplex **Uxb, **Uzb;
    extern mcomplex **HAUxb, **HAUzb;
    extern double *Qy;
    extern mcomplex ****U;
    extern mcomplex ****AU, ****AC;
    extern mcomplex *****MC;
    extern mcomplex ****LU;

    /* Local variables */
    int i, j;
    mcomplex tmp[dimR], tmp2[dimR], add[dimR];
    double flux_t;

    static double a[3] = { 1. / 3., -1. / 2, 1. / 3. };
    static double b[3] = { 2. / 3, 1. / 6, 0. };
    static double c[3] = { 1., 1. / 3., 1. / 2., };
    static double d[3] = { 0., -2. / 3, -1. / 6 };

    if (force != NULL) {
        force(0, k, 0, tmp, tmp2);
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

        /* for alpha */
        smMult0(MZ, AC[0][ALPHA], tm, RSDIAG - 1, RSDIAG - 1, dimR);
        for (i = 0; i < dimR; ++i) {
            Re(AC[0][ALPHA][i][0]) = Re(tm[i]);
            Im(AC[0][ALPHA][i][0]) = Im(tm[i]);
        }

        smMult0(Rps, MC[count][0][ALPHA], tm, RSDIAG - 1, RSDIAG - 1,
                dimR);
        for (i = 0; i < dimR; ++i) {
            Re(AC[0][ALPHA][i][0]) -= Re(tm[i]) * dt * re * (a[0] + b[0]);
            Im(AC[0][ALPHA][i][0]) -= Im(tm[i]) * dt * re * (a[0] + b[0]);
        }

        /*beta */
        smMult0(MZ, AC[0][BETA], tm, RSDIAG - 1, RSDIAG - 1, dimR);
        for (i = 0; i < dimR; ++i) {
            Re(AC[0][BETA][i][0]) = Re(tm[i]);
            Im(AC[0][BETA][i][0]) = Im(tm[i]);
        }
        smMult0(Rps, MC[count][0][BETA], tm, RSDIAG - 1, RSDIAG - 1, dimR);
        for (i = 0; i < dimR; ++i) {
            Re(AC[0][BETA][i][0]) -= Re(tm[i]) * dt * re * (a[0] + b[0]);
            Im(AC[0][BETA][i][0]) -= Im(tm[i]) * dt * re * (a[0] + b[0]);
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
            Re(fa[i]) +=
                Rw[i][j] * (Re(U[0][HXEL][j][0])) +
                Rpw[i][j] * Re(U[0][DZEL][j][0]);
        }
        if (force != NULL) {
            Re(fa[i]) += Re(tmp[i]);
            Im(fa[i]) += Im(tmp[i]);
        }

        Re(AC[0][ALPHA][i][0]) += dt * c[k] * Re(fa[i]);
    }

    /* Compute a's */
    bsolve0(MZ, AC[0][ALPHA], RSDIAG - 1, RSDIAG - 1, dimR);


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
            Re(fb[i]) +=
                Rw[i][j] * Re(U[0][HZEL][j][0]) -
                Rpw[i][j] * Re(U[0][DXEL][j][0]);
        }

        if (force != NULL) {
            Re(fb[i]) += Re(tmp2[i]);
            Im(fb[i]) += Im(tmp2[i]);
        }
        Re(AC[0][BETA][i][0]) += dt * c[k] * Re(fb[i]);
    }

    /* Compute b's */
    bsolve0(MZ, AC[0][BETA], RSDIAG - 1, RSDIAG - 1, dimR);


    flux_t = 0;

    for (i = 0; i < dimR; ++i) {
        for (j = 0; j < qpts; ++j) {
            flux_t += Rw[i][j] * Re(AC[0][ALPHA][i][0]);
        }
    }

    flux_t = -flux_t;

    for (j = 0; j < dimR; ++j)
       {
       Re(AC[0][ALPHA][j][0])= Re(AC[0][ALPHA][j][0])+flux_t*Re(add[j]) * 0.5;
       Im(AC[0][ALPHA][j][0])= Im(AC[0][ALPHA][j][0])+flux_t*Im(add[j]) * 0.5;
    }


    /* NOW COMPUTE U HATS */
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

    /* Uxb=Rp0*a, Uzb=Rp0*b, used as boundary condition for the incremental system */
    memset(Uxb[0], 0, Nz * (Nx / 2) * sizeof(mcomplex));
    memset(Uzb[0], 0, Nz * (Nx / 2) * sizeof(mcomplex));
    memset(HAUxb[0], 0, Nz * (Nx / 2) * sizeof(mcomplex));
    memset(HAUzb[0], 0, Nz * (Nx / 2) * sizeof(mcomplex));
    for (j = 0; j < dimR; ++j) {
        Re(Uxb[0][0]) += Rp0[j] * Re(AC[0][ALPHA][j][0]);
        Im(Uxb[0][0]) += Rp0[j] * Im(AC[0][ALPHA][j][0]);

        Re(Uzb[0][0]) += Rp0[j] * Re(AC[0][BETA][j][0]);
        Im(Uzb[0][0]) += Rp0[j] * Im(AC[0][BETA][j][0]);
    }
    for (j = 0; j < dimR; ++j) {
        Re(HAUxb[0][0]) += Rpp0[j] * Re(AC[0][ALPHA][j][0]);
        Im(HAUxb[0][0]) += Rpp0[j] * Im(AC[0][ALPHA][j][0]);

        Re(HAUzb[0][0]) += Rpp0[j] * Re(AC[0][BETA][j][0]);
        Im(HAUzb[0][0]) += Rpp0[j] * Im(AC[0][BETA][j][0]);
    }
    /* UPDATE RHS FOR NEXT TIME */
    if (k != 2) {               /* not last step */

        if (force != NULL) {
            force(0, k + 1, 0, tmp, tmp2);
        }

        /* MZ = M0 + (1/RE)a[k+1]dt*D0 */
        for (i = 0; i < dimR; ++i) {
            for (j = 0; j < T_RSDIAG; ++j) {
                MZ[i][j] = Rs[i][j] - re * a[k + 1] * dt * Rps[i][j];
            }
        }

        /* compute [M0 + (1/RE)a[k+1]dt*D0]*C[0][ALPHA].  Then update
           C[0][ALPHA].  Do the same for the b's. */
        smMult0(MZ, AC[0][ALPHA], tm, RSDIAG - 1, RSDIAG - 1, dimR);

        memset(fa, 0, dimR * sizeof(mcomplex));
        for (i = 0; i < dimR; ++i) {
            for (j = 0; j < qpts; ++j) {
                Re(fa[i]) +=
                    Rw[i][j] * (Re(LU[0][HXEL][j][0])) +
                    Rpw[i][j] * Re(LU[0][DZEL][j][0]);
            }
            if (force != NULL) {
                Re(fa[i]) += Re(tmp[i]);
                Im(fa[i]) += Im(tmp[i]);
            }

        }
        for (i = 0; i < dimR; ++i) {
            Re(AC[0][ALPHA][i][0]) = Re(tm[i]) + dt * d[k + 1] * Re(fa[i]);
            Im(AC[0][ALPHA][i][0]) = Im(tm[i]) + dt * d[k + 1] * Im(fa[i]);
        }
        smMult0(Rps, MC[count - 1][0][ALPHA], tm, RSDIAG - 1, RSDIAG - 1,
                dimR);
        for (i = 0; i < dimR; ++i) {
            Re(AC[0][ALPHA][i][0]) -=
                Re(tm[i]) * dt * re * (a[k + 1] + b[k + 1]);
            Im(AC[0][ALPHA][i][0]) -=
                Im(tm[i]) * dt * re * (a[k + 1] + b[k + 1]);
        }

        /* now b's */
        smMult0(MZ, AC[0][BETA], tm, RSDIAG - 1, RSDIAG - 1, dimR);
        memset(fb, 0, dimR * sizeof(mcomplex));
        for (i = 0; i < dimR; ++i) {
            for (j = 0; j < qpts; ++j) {
                Re(fb[i]) +=
                    Rw[i][j] * Re(LU[0][HZEL][j][0]) -
                    Rpw[i][j] * Re(LU[0][DXEL][j][0]);
            }
            if (force != NULL) {
                Re(fb[i]) += Re(tmp2[i]);
                Im(fb[i]) += Im(tmp2[i]);
            }

        }


        for (i = 0; i < dimR; ++i) {
            Re(AC[0][BETA][i][0]) = Re(tm[i]) + dt * d[k + 1] * Re(fb[i]);
            Im(AC[0][BETA][i][0]) = Im(tm[i]) + dt * d[k + 1] * Im(fb[i]);
        }
        smMult0(Rps, MC[count - 1][0][BETA], tm, RSDIAG - 1, RSDIAG - 1,
                dimR);
        for (i = 0; i < dimR; ++i) {
            Re(AC[0][BETA][i][0]) -=
                Re(tm[i]) * dt * re * (a[k + 1] + b[k + 1]);
            Im(AC[0][BETA][i][0]) -=
                Im(tm[i]) * dt * re * (a[k + 1] + b[k + 1]);
        }
    }
}                               /* end project zero (Kx,Kz) */

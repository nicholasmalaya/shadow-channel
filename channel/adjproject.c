
#include <stdio.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include "minChnl.h"
#include "mvOps.h"

/* project when (Kx,Kz) != (0,0) */
void adjproject(int count, int k, int z, int x0, func_force_t force)
{
    /* External Variables */
    extern int qpts, dimR, dimQ, Nx;
    extern double dt, re;
    extern double *Kx, *Kz, **K2;
    extern mcomplex **Fa, **Fb, **TM;
    extern double **Q, **Qp, **Qpp, **R, **Rp, **Qw, **Qpw, **Rw, **Qs,
        **Qps, **Qpps, **Rs, **Rps, *Rp0, **Qppw, **Rpw, *Qy, *Rpp0;
    extern double ***M;
    extern mcomplex **Uxb, **Uzb;
    extern mcomplex **HAUxb, **HAUzb;
    extern mcomplex ****U, ****C;
    extern mcomplex ****AU, ****AC;
    extern mcomplex *****MC, ****LU;;
    /* Local variables */
    int i, j, x;
    double s, t[2];
    /*  static double e[3] = { 1./3., 7./15, 1.};
       static double a[3] = {1./6, -3./40., 29./96.};
       static double b[3] = { 5./24.,37./160.,   1./6.};
       static double c[3] = {   3./4., 5./12., 8./15.};
       static double d[3] = { 0.,  -5./12.,  -17./60.}; */

    static double a[3] = { 1. / 3., -1. / 2, 1. / 3. };
    static double b[3] = { 2. / 3, 1. / 6, 0. };
    static double c[3] = { 1., 1. / 3., 1. / 2., };
    static double d[3] = { 0., -2. / 3, -1. / 6 };
    static double e[3] = { 1. / 3., 1. / 2, 1. };
    mcomplex tmp[Nx / 2][dimR], tmp2[Nx / 2][dimQ];

    if (force != NULL) {
        force(0, k, z, tmp[0], tmp2[0]);
    }
    /* FIRST COMPUTE ALPHAS */
    /* Create matrices for solving linear system. 
       Right hand side of system:  If this is the first step in the Runge-Kutta
       scheme, compute 
       [Mv + (1/RE)a[k]dt*Dv]*C[z][ALPHA] + dt*c[k]Fv. 
       Otherwise, compute C[z][ALPHA] + dt*c[k]Fv. 
       Lhs matrix:  Mv - (1/RE)b[k]dt*Dv                                */
    if (k == 0) {               /* first step */
        memset(M[0][0], 0, dimR * 9 * (Nx / 2) * sizeof(double));
        for (i = 0; i < dimQ; ++i) {    /* M = Mv + (1/RE)a[k]dt*Dv */
            for (j = 0; j < T_QSDIAG; ++j) {
                for (x = x0; x < Nx / 2; ++x) {
                    s = K2[z][x] * K2[z][x];
                    M[i][j][x] = -(K2[z][x] * Qs[i][j] + Qps[i][j]) +
                        re * a[0] * dt * (s * Qs[i][j] +
                                          2. * K2[z][x] * Qps[i][j] +
                                          Qpps[i][j]);
                }
            }
        }
        /* compute [Mv + (1/RE)a[k]dt*Dv]*C[z][ALPHA] and store the result
           in TM.  Then transfer the result back to C[z][ALPHA]. */
        smMult(M, AC[z][ALPHA], TM, QSDIAG - 1, QSDIAG - 1, dimQ, Nx / 2,
               x0);
        for (i = 0; i < dimQ; ++i) {
            memcpy(&AC[z][ALPHA][i][x0], &TM[i][x0],
                   (Nx / 2 - x0) * sizeof(mcomplex));
        }

        memset(M[0][0], 0, dimR * 9 * (Nx / 2) * sizeof(double));
        for (i = 0; i < dimQ; ++i) {
            for (j = 0; j < T_QSDIAG; ++j) {
                for (x = x0; x < Nx / 2; ++x) {
                    s = K2[z][x] * K2[z][x];
                    M[i][j][x] =
                        (s * Qs[i][j] + 2. * K2[z][x] * Qps[i][j] +
                         Qpps[i][j]);
                }
            }
        }
        smMult(M, MC[count][z][ALPHA], TM, QSDIAG - 1, QSDIAG - 1, dimQ,
               Nx / 2, x0);
        for (i = 0; i < dimQ; ++i) {
            for (x = x0; x < Nx / 2; x++) {
                Re(AC[z][ALPHA][i][x]) +=
                    Re(TM[i][x]) * re * dt * (a[0] + b[0]);
                Im(AC[z][ALPHA][i][x]) +=
                    Im(TM[i][x]) * re * dt * (a[0] + b[0]);
            }
        }

    }

    /* Left hand side M = Mv - (1/RE)b[k]dt*Dv */
    memset(M[0][0], 0, dimR * 9 * (Nx / 2) * sizeof(double));
    for (i = 0; i < dimQ; ++i) {
        for (j = 0; j < T_QSDIAG; ++j) {
            for (x = x0; x < Nx / 2; ++x) {
                s = K2[z][x] * K2[z][x];
                M[i][j][x] = -(K2[z][x] * Qs[i][j] + Qps[i][j]) -
                    re * b[k] * dt * (s * Qs[i][j] +
                                      2. * K2[z][x] * Qps[i][j] +
                                      Qpps[i][j]);
            }
        }
    }

    /* Finish computing the right hand size */
    /* array Fa */
    memset(Fa[0], 0, dimR * (Nx / 2) * sizeof(mcomplex));
    for (i = 0; i < dimQ; ++i) {
        for (j = 0; j < qpts; ++j) {
            for (x = x0; x < Nx / 2; ++x) {
                Re(Fa[i][x]) += (Qpw[i][j] *
                                 (-Kx[x] * Im(U[z][HXEL][j][x]) -
                                  Kz[z] * Im(U[z][HZEL][j][x])) -
                                 K2[z][x] * (Qw[i][j] *
                                             Re(U[z][HYEL][j][x]))) +
                    Qppw[i][j] * (Kz[z] * Im(U[z][DXEL][j][x]) -
                                  Kx[x] * Im(U[z][DZEL][j][x]));

                Im(Fa[i][x]) += (Qpw[i][j] *
                                 (Kx[x] * Re(U[z][HXEL][j][x]) +
                                  Kz[z] * Re(U[z][HZEL][j][x])) -
                                 K2[z][x] * (Qw[i][j] *
                                             Im(U[z][HYEL][j][x]))) +
                    Qppw[i][j] * (-Kz[z] * Re(U[z][DXEL][j][x]) +
                                  Kx[x] * Re(U[z][DZEL][j][x]));
            }
        }
    }

    if (force != NULL) {
        for (x = x0; x < Nx / 2; ++x) {
            for (i = 0; i < dimQ; ++i) {
                Im(Fa[i][x]) += Im(tmp2[x][i]);
                Re(Fa[i][x]) += Re(tmp2[x][i]);
            }
        }
    }

    for (i = 0; i < dimQ; ++i) {
        for (x = x0; x < Nx / 2; ++x) {
            Re(AC[z][ALPHA][i][x]) += dt * c[k] * Re(Fa[i][x]);
            Im(AC[z][ALPHA][i][x]) += dt * c[k] * Im(Fa[i][x]);
        }
    }

    /* Compute alphas */
    bsolve(M, AC[z][ALPHA], QSDIAG - 1, QSDIAG - 1, dimQ, Nx / 2, x0);

    /* NOW COMPUTE BETAS */
    /* Create matrices for solving linear system. 
       Right hand side of system:  If this is the first step in the Runge-Kutta
       scheme, compute 
       [Mg + (1/RE)a[k]dt*Dg]*C[z][BETA] + dt*c[k]Fg. 
       Otherwise, compute C[z][BETA] + dt*c[k]Fg. 
       Lhs matrix:  Mg - (1/RE)b[k]dt*Dg                                */

    if (k == 0) {               /* first step */
        /* M = Mg + (1/RE)a[k]dt*Dg */
        memset(M[0][0], 0, dimR * 9 * (Nx / 2) * sizeof(double));
        for (i = 0; i < dimR; ++i) {
            for (j = 0; j < T_RSDIAG; ++j) {
                for (x = x0; x < Nx / 2; ++x) {
                    M[i][j][x] = Rs[i][j] -
                        re * a[0] * dt * (Rps[i][j] + K2[z][x] * Rs[i][j]);
                }
            }
        }

        /* compute [Mv + (1/RE)a[k]dt*Dv]*C[z][BETA] and store the result
           in TM.  Then transfer the result back to C[z][ALPHA]. */
        smMult(M, AC[z][BETA], TM, RSDIAG - 1, RSDIAG - 1, dimR, Nx / 2,
               x0);
        for (i = 0; i < dimR; ++i) {
            memcpy(&AC[z][BETA][i][x0], &TM[i][x0],
                   (Nx / 2 - x0) * sizeof(mcomplex));
        }

        memset(M[0][0], 0, dimR * 9 * (Nx / 2) * sizeof(double));
        for (i = 0; i < dimR; ++i) {
            for (j = 0; j < T_RSDIAG; ++j) {
                for (x = x0; x < Nx / 2; ++x) {
                    M[i][j][x] = -(Rps[i][j] + K2[z][x] * Rs[i][j]);
                }
            }
        }

        smMult(M, MC[count][z][BETA], TM, RSDIAG - 1, RSDIAG - 1, dimR,
               Nx / 2, x0);
        for (i = 0; i < dimR; ++i) {
            for (x = x0; x < Nx / 2; x++) {
                Re(AC[z][BETA][i][x]) +=
                    Re(TM[i][x]) * re * dt * (a[0] + b[0]);
                Im(AC[z][BETA][i][x]) +=
                    Im(TM[i][x]) * re * dt * (a[0] + b[0]);
            }
        }
    }

    /* left hand side M = Mg - (1/RE)b[k]dt*Dg */
    memset(M[0][0], 0, dimR * 9 * (Nx / 2) * sizeof(double));
    for (i = 0; i < dimR; ++i) {        /* M = Mg - (1/RE)b[k]dt*Dg */
        for (j = 0; j < T_RSDIAG; ++j) {
            for (x = x0; x < Nx / 2; ++x) {
                M[i][j][x] = Rs[i][j] +
                    re * b[k] * dt * (Rps[i][j] + K2[z][x] * Rs[i][j]);
            }
        }
    }

    /* Finish computing the right hand size */
    /* array Fb */
    memset(Fb[0], 0, dimR * (Nx / 2) * sizeof(mcomplex));
    for (i = 0; i < dimR; ++i) {
        for (j = 0; j < qpts; ++j) {
            for (x = x0; x < Nx / 2; ++x) {
                Re(Fb[i][x]) += Rw[i][j] *
                    (-Kz[z] * Im(U[z][HXEL][j][x]) +
                     Kx[x] * Im(U[z][HZEL][j][x])) +
                    Rpw[i][j] * (-Kz[z] * Im(U[z][DZEL][j][x]) -
                                 Kx[x] * Im(U[z][DXEL][j][x]));
                Im(Fb[i][x]) +=
                    Rw[i][j] * (Kz[z] * Re(U[z][HXEL][j][x]) -
                                Kx[x] * Re(U[z][HZEL][j][x])) +
                    Rpw[i][j] * (Kz[z] * Re(U[z][DZEL][j][x]) +
                                 Kx[x] * Re(U[z][DXEL][j][x]));
            }
        }
    }
    if (force != NULL) {
        for (x = x0; x < Nx / 2; ++x) {
            for (i = 0; i < dimR; ++i) {
                Im(Fb[i][x]) += Im(tmp[x][i]);
                Re(Fb[i][x]) += Re(tmp[x][i]);
            }
        }
    }


    for (i = 0; i < dimR; ++i) {
        for (x = x0; x < Nx / 2; ++x) {
            Re(AC[z][BETA][i][x]) += dt * c[k] * Re(Fb[i][x]);
            Im(AC[z][BETA][i][x]) += dt * c[k] * Im(Fb[i][x]);
        }
    }

    /* Compute betas */
    bsolve(M, AC[z][BETA], RSDIAG - 1, RSDIAG - 1, dimR, Nx / 2, x0);


    /* NOW COMPUTE U HATS */
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

    /* compute dux duz at y=-1 as the boundary condition for the incremental systems */

    for (x = x0; x < Nx / 2; ++x) {
        for (j = 0; j < dimR; ++j) {
            Re(Uxb[z][x]) +=
                Rp0[j] * Im(AC[z][BETA][j][x]) * Kz[z] / K2[z][x];
            Im(Uxb[z][x]) +=
                -Rp0[j] * Re(AC[z][BETA][j][x]) * Kz[z] / K2[z][x];

            Re(Uzb[z][x]) +=
                -Rp0[j] * Im(AC[z][BETA][j][x]) * Kx[x] / K2[z][x];
            Im(Uzb[z][x]) +=
                Rp0[j] * Re(AC[z][BETA][j][x]) * Kx[x] / K2[z][x];
        }
    }
    for (x = x0; x < Nx / 2; ++x) {
        for (j = 0; j < dimR; ++j) {
            Re(HAUxb[z][x]) +=
                Rpp0[j] * Im(AC[z][BETA][j][x]) * Kz[z] / K2[z][x];
            Im(HAUxb[z][x]) +=
                -Rpp0[j] * Re(AC[z][BETA][j][x]) * Kz[z] / K2[z][x];

            Re(HAUzb[z][x]) +=
                -Rpp0[j] * Im(AC[z][BETA][j][x]) * Kx[x] / K2[z][x];
            Im(HAUzb[z][x]) +=
                Rpp0[j] * Re(AC[z][BETA][j][x]) * Kx[x] / K2[z][x];
        }
    }
    /* UPDATE RHS FOR NEXT TIME */
    if (k != 2) {               /* not last step */
        if (force != NULL) {
            force(0, k + 1, z, tmp[0], tmp2[0]);
        }

        /* first alphas */
        memset(M[0][0], 0, dimR * 9 * (Nx / 2) * sizeof(double));
        for (i = 0; i < dimQ; ++i) {    /* M = Mv + (1/RE)a[k+1]dt*Dv */
            for (j = 0; j < T_QSDIAG; ++j) {
                for (x = x0; x < Nx / 2; ++x) {
                    s = K2[z][x] * K2[z][x];
                    M[i][j][x] = -(K2[z][x] * Qs[i][j] + Qps[i][j]) +
                        re * a[k + 1] * dt * (s * Qs[i][j] +
                                              2. * K2[z][x] * Qps[i][j] +
                                              Qpps[i][j]);
                }
            }
        }

        /* compute [Mv + (1/RE)a[k]dt*Dv]*C[z][ALPHA].  Then update
           C[z][ALPHA] */
        smMult(M, AC[z][ALPHA], TM, QSDIAG - 1, QSDIAG - 1, dimQ, Nx / 2,
               x0);
        memset(Fa[0], 0, dimR * (Nx / 2) * sizeof(mcomplex));
        for (i = 0; i < dimQ; ++i) {
            for (j = 0; j < qpts; ++j) {
                for (x = x0; x < Nx / 2; ++x) {
                    Re(Fa[i][x]) += (Qpw[i][j] *
                                     (-Kx[x] * Im(LU[z][HXEL][j][x]) -
                                      Kz[z] * Im(LU[z][HZEL][j][x])) -
                                     K2[z][x] * (Qw[i][j] *
                                                 Re(LU[z][HYEL][j][x]))) +
                        Qppw[i][j] * (Kz[z] * Im(LU[z][DXEL][j][x]) -
                                      Kx[x] * Im(LU[z][DZEL][j][x]));

                    Im(Fa[i][x]) += (Qpw[i][j] *
                                     (Kx[x] * Re(LU[z][HXEL][j][x]) +
                                      Kz[z] * Re(LU[z][HZEL][j][x])) -
                                     K2[z][x] * (Qw[i][j] *
                                                 Im(LU[z][HYEL][j][x]))) +
                        Qppw[i][j] * (-Kz[z] * Re(LU[z][DXEL][j][x]) +
                                      Kx[x] * Re(LU[z][DZEL][j][x]));
                }
            }
        }
        if (force != NULL) {
            for (x = x0; x < Nx / 2; ++x) {
                for (i = 0; i < dimQ; ++i) {
                    Im(Fa[i][x]) += Im(tmp2[x][i]);
                    Re(Fa[i][x]) += Re(tmp2[x][i]);
                }
            }
        }

        for (i = 0; i < dimQ; ++i) {
            for (x = x0; x < Nx / 2; ++x) {
                Re(AC[z][ALPHA][i][x]) =
                    Re(TM[i][x]) + dt * d[k + 1] * Re(Fa[i][x]);
                Im(AC[z][ALPHA][i][x]) =
                    Im(TM[i][x]) + dt * d[k + 1] * Im(Fa[i][x]);
            }
        }
        memset(M[0][0], 0, dimR * 9 * (Nx / 2) * sizeof(double));
        for (i = 0; i < dimQ; ++i) {
            for (j = 0; j < T_QSDIAG; ++j) {
                for (x = x0; x < Nx / 2; ++x) {
                    s = K2[z][x] * K2[z][x];
                    M[i][j][x] =
                        (s * Qs[i][j] + 2. * K2[z][x] * Qps[i][j] +
                         Qpps[i][j]);
                }
            }
        }
        smMult(M, MC[count - 1][z][ALPHA], TM, QSDIAG - 1, QSDIAG - 1,
               dimQ, Nx / 2, x0);

        for (i = 0; i < dimQ; ++i) {
            for (x = x0; x < Nx / 2; x++) {
                Re(AC[z][ALPHA][i][x]) +=
                    Re(TM[i][x]) * re * dt * (a[k + 1] + b[k + 1]);
                Im(AC[z][ALPHA][i][x]) +=
                    Im(TM[i][x]) * re * dt * (a[k + 1] + b[k + 1]);
            }
        }

        /* now betas */
        memset(M[0][0], 0, dimR * 9 * (Nx / 2) * sizeof(double));
        for (i = 0; i < dimR; ++i) {    /* M = Mg + (1/RE)a[k+1]dt*Dg */
            for (j = 0; j < T_RSDIAG; ++j) {
                for (x = x0; x < Nx / 2; ++x) {
                    M[i][j][x] = Rs[i][j] -
                        re * a[k + 1] * dt * (Rps[i][j] +
                                              K2[z][x] * Rs[i][j]);
                }
            }
        }

        /* compute [Mg + (1/RE)a[k+1]dt*Dg]*C[z][BETA] and then
           update C[z][BETA] */
        smMult(M, AC[z][BETA], TM, RSDIAG - 1, RSDIAG - 1, dimR, Nx / 2,
               x0);
        memset(Fb[0], 0, dimR * (Nx / 2) * sizeof(mcomplex));
        for (i = 0; i < dimR; ++i) {
            for (j = 0; j < qpts; ++j) {
                for (x = x0; x < Nx / 2; ++x) {
                    Re(Fb[i][x]) += Rw[i][j] *
                        (-Kz[z] * Im(LU[z][HXEL][j][x]) +
                         Kx[x] * Im(LU[z][HZEL][j][x])) +
                        Rpw[i][j] * (-Kz[z] * Im(LU[z][DZEL][j][x]) -
                                     Kx[x] * Im(LU[z][DXEL][j][x]));
                    Im(Fb[i][x]) +=
                        Rw[i][j] * (Kz[z] * Re(LU[z][HXEL][j][x]) -
                                    Kx[x] * Re(LU[z][HZEL][j][x])) +
                        Rpw[i][j] * (Kz[z] * Re(LU[z][DZEL][j][x]) +
                                     Kx[x] * Re(LU[z][DXEL][j][x]));
                }
            }
        }

        if (force != NULL) {
            for (x = x0; x < Nx / 2; ++x) {
                for (i = 0; i < dimR; ++i) {
                    Im(Fb[i][x]) += Im(tmp[x][i]);
                    Re(Fb[i][x]) += Re(tmp[x][i]);
                }
            }
        }

        for (i = 0; i < dimR; ++i) {
            for (x = x0; x < Nx / 2; ++x) {
                Re(AC[z][BETA][i][x]) =
                    Re(TM[i][x]) + dt * d[k + 1] * Re(Fb[i][x]);
                Im(AC[z][BETA][i][x]) =
                    Im(TM[i][x]) + dt * d[k + 1] * Im(Fb[i][x]);
            }
        }

        memset(M[0][0], 0, dimR * 9 * (Nx / 2) * sizeof(double));
        for (i = 0; i < dimR; ++i) {
            for (j = 0; j < T_RSDIAG; ++j) {
                for (x = x0; x < Nx / 2; ++x) {
                    M[i][j][x] = -(Rps[i][j] + K2[z][x] * Rs[i][j]);
                }
            }
        }
        smMult(M, MC[count - 1][z][BETA], TM, RSDIAG - 1, RSDIAG - 1, dimR,
               Nx / 2, x0);
        for (i = 0; i < dimR; ++i) {
            for (x = x0; x < Nx / 2; x++) {
                Re(AC[z][BETA][i][x]) +=
                    Re(TM[i][x]) * re * dt * (a[k + 1] + b[k + 1]);
                Im(AC[z][BETA][i][x]) +=
                    Im(TM[i][x]) * re * dt * (a[k + 1] + b[k + 1]);
            }
        }

    }
}                               /* end project for nonzero (Kx,Kz) */

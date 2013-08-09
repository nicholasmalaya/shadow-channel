#include <math.h>
#include <stdio.h>
#include <string.h>
#include "minChnl.h"

/**************forces used for manufacture solution 1 ***************************/
void adjforce0(int n, int k, int flag, mcomplex * f, mcomplex * g)
{
    extern double *Qy, **Rw;
    extern double dt, re;
    extern int dimR, qpts;
    static double e[3] = { 1. / 3., 1. / 2, 1. };

    int i, j;
    for (i = 0; i < dimR; ++i) {
        Re(f[i]) = 0;
        Im(f[i]) = 0;
        Re(g[i]) = 0;
        Im(g[i]) = 0;
        if (flag == 0) {
            for (j = 0; j < qpts; ++j) {
                Re(f[i]) +=
                    Rw[i][j] * (-2 * M_PI * (1 - Qy[j] * Qy[j]) *
                                cos(2 * M_PI *
                                    (n * dt - dt * e[k] - 0.01)))
                    +
                    Rw[i][j] * 2 * re * sin(2 * M_PI *
                                            (n * dt - dt * e[k] - 0.01));
            }
        }
    }
}


void increadjforce0(int n, int k, int flag, mcomplex * f, mcomplex * g)
{
    extern double *Qy, **Rw;
    extern double dt, re;
    extern int dimR, qpts;
    static double e[3] = { 1. / 3., 1. / 2, 1. };
    int i, j;
    for (i = 0; i < dimR; ++i) {
        Re(f[i]) = 0;
        Im(f[i]) = 0;
        Re(g[i]) = 0;
        Im(g[i]) = 0;
        for (j = 0; j < qpts; ++j) {
            Re(f[i]) +=
                Rw[i][j] * (-2 * M_PI * (2 - Qy[j] - Qy[j] * Qy[j]) *
                            cos(2 * M_PI * (n * dt - dt * e[k] - 0.01)))
                +
                Rw[i][j] * 2 * re * sin(2 * M_PI *
                                        (n * dt - dt * e[k] - 0.01));
        }
    }

}

  //  +Rw[i][j]*(-2*M_PI*(2-Qy[j]-Qy[j]*Qy[j])*cos(2*M_PI*(n*dt-dt*e[k]-0.01)))
            //+2*re  +2*re*sin(2*M_PI*(n*dt-dt*e[k]-0.01))) ;
            // +Rw[i][j]*2*re*sin(2*M_PI*(n*dt-dt*e[k]-0.01));
void increadjforce(int n, int k, int z, mcomplex * f, mcomplex * g)
{
    extern double *Qy, **Rw, *Kz;
    extern double dt, re;
    extern int dimR, qpts, dimQ;
    static double e[3] = { 1. / 3., 1. / 2, 1. };
    extern int Nz, Nx;
    int i, j, x;


    for (x = 0; x < Nx / 2; x++) {
        for (i = 0; i < dimR; ++i) {
            Re(f[x * dimR + i]) = 0;
            Im(f[x * dimR + i]) = 0;
            if (x == 0 && (z == 1 || z == Nz - 1)) {
                for (j = 0; j < qpts; j++) {
                    Im(f[x * dimR + i]) +=
                        Rw[i][j] * Kz[z] * (-2 * M_PI * (1 - Qy[j]) *
                                            cos(2 * M_PI *
                                                (n * dt - dt * e[k] -
                                                 0.01))
                                            + re * Kz[z] * Kz[z] * (1 -
                                                                    Qy[j])
                                            * sin(2 * M_PI *
                                                  (n * dt - dt * e[k] -
                                                   0.01)));
                }
            }
        }
    }

    for (x = 0; x < Nx / 2; x++) {
        for (i = 0; i < dimQ; ++i) {
            Re(g[x * dimQ + i]) = 0.;
            Im(g[x * dimQ + i]) = 0.;
        }
    }
}


/**************forces used for manufacture solution 2 ***************************/
void adjforce_2(int n, int k, int z, mcomplex * f, mcomplex * g)
{
    extern double *Qy, **Rw, *Kz;
    extern double dt, re;
    extern int dimR, qpts, dimQ;
    static double e[3] = { 1. / 3., 1. / 2, 1. };
    extern int Nz, Nx;
    int i, j, x;


    for (x = 0; x < Nx / 2; x++) {
        for (i = 0; i < dimR; ++i) {
            Re(f[x * dimR + i]) = 0;
            Im(f[x * dimR + i]) = 0;
            if (x == 1 && z == 1) {
                for (j = 0; j < qpts; j++) {
                    Im(f[x * dimR + i]) +=
                        Rw[i][j] * 2 * Kz[z] * (-(1 - Qy[j] * Qy[j]) *
                                                cos(2 * M_PI *
                                                    (n * dt - dt * e[k] -
                                                     0.01)) * 2 * M_PI +
                                                re * (2 +
                                                      2 * Kz[z] * Kz[z] *
                                                      (1 -
                                                       Qy[j] * Qy[j])) *
                                                sin(2 * M_PI *
                                                    (n * dt - dt * e[k] -
                                                     0.01)));
                }
            }
        }
    }

    for (x = 0; x < Nx / 2; x++) {
        for (i = 0; i < dimQ; ++i) {
            Re(g[x * dimQ + i]) = 0.;
            Im(g[x * dimQ + i]) = 0.;
        }
    }
}

  /*for (i = 0; i < dimQ; ++i)
     {
     for (j = 0; j < qpts; ++j)
     {
     if(z==1)
     {
     Re(Fa[i][1])+=Qw[i][j]*k11*k11*2*2*Qy[j]*(1-Qy[j]*Qy[j])*sin(2*M_PI*(n*dt-dt*e[k]-0.01))
     *sin(2*M_PI*(n*dt-dt*e[k]));
     }
     }
     } */

    /* for (i = 0; i < dimQ; ++i)
       {
       for (j = 0; j < qpts; ++j)
       {
       if(z==1)
       {
       Im(Fa[i][1])+=Qw[i][j]*k11*(-(12*Qy[j]*Qy[j]-4-k11*k11*2*(1-Qy[j]*Qy[j])*(1-Qy[j]*Qy[j]))
       *cos(2*M_PI*(n*dt-dt*e[k]-0.01))*2*M_PI
       -re*(24-2*k11*k11*2*(12*Qy[j]*Qy[j]-4)+4*k11*k11*k11*k11*(1-Qy[j]*Qy[j])*(1-Qy[j]*Qy[j]))
       *sin(2*M_PI*(n*dt-dt*e[k]-0.01)));
       Re(Fa[i][1])+=-Qw[i][j]*(4*(1-Qy[j]*Qy[j])-32*Qy[j]*Qy[j]+2*k11*k11*(1-Qy[j]*Qy[j])*(1-Qy[j]*Qy[j])
       -4*k11*Qy[j])*k11*k11*(1-Qy[j]*Qy[j])*sin(2*M_PI*(n*dt-dt*e[k]-0.01))
       *sin(2*M_PI*(n*dt-dt*e[k]));
       }
       }
       } */
void increadjforce0_2(int n, int k, int flag, mcomplex * f, mcomplex * g)
{
    extern double *Qy, **Rw;
    extern double dt, re;
    extern int dimR, qpts;
    static double e[3] = { 1. / 3., 1. / 2, 1. };
    int i, j;
    for (i = 0; i < dimR; ++i) {
        Re(f[i]) = 0;
        Im(f[i]) = 0;
        Re(g[i]) = 0;
        Im(g[i]) = 0;
        for (j = 0; j < qpts; ++j) {
            Re(f[i]) +=
                -Rw[i][j] * 4 * M_PI * (1 -
                                        Qy[j]) * cos(2 * M_PI * (n * dt -
                                                                 dt *
                                                                 e[k] -
                                                                 0.01));
            Re(g[i]) +=
                Rw[i][j] * 4 * M_PI * (1 -
                                       Qy[j]) * cos(2 * M_PI * (n * dt -
                                                                dt * e[k] -
                                                                0.01));
        }
    }

}

void increadjforce_2(int n, int k, int z, mcomplex * f, mcomplex * g)
{
    extern double *Qy, **Rw, *Kz, **Qw;
    extern double dt, re;
    extern int dimR, qpts, dimQ;
    static double e[3] = { 1. / 3., 1. / 2, 1. };
    extern int Nz, Nx;
    int i, j, x;


    for (x = 0; x < Nx / 2; x++) {
        for (i = 0; i < dimR; ++i) {
            Re(f[x * dimR + i]) = 0;
            Im(f[x * dimR + i]) = 0;
            if (x == 2 && z == 2) {
                for (j = 0; j < qpts; j++) {
                    Im(f[x * dimR + i]) +=
                        Rw[i][j] * Kz[z] * (-2 * M_PI *
                                            (2 * (1 - Qy[j]) +
                                             (1 -
                                              Qy[j] * Qy[j]) * 4 * Qy[j]) *
                                            cos(2 * M_PI *
                                                (n * dt - dt * e[k] -
                                                 0.01))
                                            + re * (24 * Qy[j] +
                                                    2 * Kz[z] * Kz[z] *
                                                    (2 * (1 - Qy[j]) +
                                                     (1 -
                                                      Qy[j] * Qy[j]) * 4 *
                                                     Qy[j])) * sin(2 *
                                                                   M_PI *
                                                                   (n *
                                                                    dt -
                                                                    dt *
                                                                    e[k] -
                                                                    0.01)));
                }
            }
        }
    }

    for (x = 0; x < Nx / 2; x++) {
        for (i = 0; i < dimQ; ++i) {
            Re(g[x * dimQ + i]) = 0.;
            Im(g[x * dimQ + i]) = 0.;
            if (x == 2 && z == 2) {
                for (j = 0; j < qpts; j++) {
                    Im(g[x * dimQ + i]) +=
                        Qw[i][j] * Kz[z] *
                        (-
                         (12 * Qy[j] * Qy[j] - 4 -
                          2 * Kz[z] * Kz[z] * (1 - Qy[j] * Qy[j]) * (1 -
                                                                     Qy[j]
                                                                     *
                                                                     Qy
                                                                     [j]))
                         * 2 * M_PI * cos(2 * M_PI *
                                          (n * dt - dt * e[k] - 0.01)) -
                         re * (24 -
                               2 * 2 * Kz[z] * Kz[z] * (12 * Qy[j] *
                                                        Qy[j] - 4) +
                               4 * Kz[z] * Kz[z] * Kz[z] * Kz[z] * (1 -
                                                                    Qy[j] *
                                                                    Qy[j])
                               * (1 -
                                  Qy[j] * Qy[j])) * sin(2 * M_PI * (n *
                                                                    dt -
                                                                    dt *
                                                                    e[k] -
                                                                    0.01)));
                }
            }
        }
    }
}

/**************forces used for manufacture solution 3 ***************************/
void adjforce0_3(int n, int k, int flag, mcomplex * f, mcomplex * g)
{
    extern double *Qy, **Rw;
    extern double dt, re;
    extern int dimR, qpts;
    static double e[3] = { 1. / 3., 1. / 2, 1. };

    int i, j;
    for (i = 0; i < dimR; ++i) {
        Re(f[i]) = 0;
        Im(f[i]) = 0;
        Re(g[i]) = 0;
        Im(g[i]) = 0;
        if (flag == 0) {
            for (j = 0; j < qpts; ++j) {
                Re(f[i]) += Rw[i][j] * 2 * re
                    * sin(2 * M_PI * (n * dt - dt * e[k]));
            }
        }
    }
}

void adjforce_3(int n, int k, int z, mcomplex * f, mcomplex * g)
{
    extern double *Qy, **Rw, *Kz, **Qw;
    extern double dt, re;
    extern int dimR, qpts, dimQ;
    static double e[3] = { 1. / 3., 1. / 2, 1. };
    extern int Nz, Nx;
    int i, j, x;


    for (x = 0; x < Nx / 2; x++) {
        for (i = 0; i < dimR; ++i) {
            Re(f[x * dimR + i]) = 0;
            Im(f[x * dimR + i]) = 0;
            if (x == 1 && z == 1) {
                for (j = 0; j < qpts; j++) {
                    Im(f[x * dimR + i]) +=
                        Rw[i][j] * 2 * Kz[z] * (-(1 - Qy[j] * Qy[j]) *
                                                cos(2 * M_PI *
                                                    (n * dt - dt * e[k] -
                                                     0.01)) * 2 * M_PI +
                                                re * (2 +
                                                      2 * Kz[z] * Kz[z] *
                                                      (1 -
                                                       Qy[j] * Qy[j])) *
                                                sin(2 * M_PI *
                                                    (n * dt - dt * e[k] -
                                                     0.01)));
                    Re(f[x * dimR + i]) +=
                        Rw[i][j] * 2 * Kz[z] * Kz[z] * (1 -
                                                        Qy[j] * Qy[j]) *
                        (1 -
                         Qy[j] * Qy[j]) * sin(2 * M_PI * (n * dt -
                                                          dt * e[k] -
                                                          0.01))
                        * sin(2 * M_PI * (n * dt - dt * e[k]));
                }
            }
        }
    }

    for (x = 0; x < Nx / 2; x++) {
        for (i = 0; i < dimQ; ++i) {
            Re(g[x * dimQ + i]) = 0.;
            Im(g[x * dimQ + i]) = 0.;
            if (x == 1 && z == 1) {
                for (j = 0; j < qpts; j++) {
                    Re(g[x * dimQ + i]) +=
                        Qw[i][j] * Kz[z] * Kz[z] * 2 * 2 * Qy[j] * (1 -
                                                                    Qy[j] *
                                                                    Qy[j])
                        * sin(2 * M_PI * (n * dt - dt * e[k] - 0.01))
                        * sin(2 * M_PI * (n * dt - dt * e[k]));
                }
            }
        }
    }
}

void increadjforce0_3(int n, int k, int flag, mcomplex * f, mcomplex * g)
{
    extern double *Qy, **Rw;
    extern double dt, re;
    extern int dimR, qpts;
    static double e[3] = { 1. / 3., 1. / 2, 1. };
    int i, j;
    for (i = 0; i < dimR; ++i) {
        Re(f[i]) = 0;
        Im(f[i]) = 0;
        Re(g[i]) = 0;
        Im(g[i]) = 0;
        for (j = 0; j < qpts; ++j) {
            Re(f[i]) +=
                -Rw[i][j] * 4 * M_PI * (1 -
                                        Qy[j]) * cos(2 * M_PI * (n * dt -
                                                                 dt *
                                                                 e[k] -
                                                                 0.01));
            Re(g[i]) +=
                Rw[i][j] * 4 * M_PI * (1 -
                                       Qy[j]) * cos(2 * M_PI * (n * dt -
                                                                dt * e[k] -
                                                                0.01));
        }
    }

}


void increadjforce_3(int n, int k, int z, mcomplex * f, mcomplex * g)
{
    extern double *Qy, **Rw, *Kz, **Qw;
    extern double dt, re;
    extern int dimR, qpts, dimQ;
    static double e[3] = { 1. / 3., 1. / 2, 1. };
    extern int Nz, Nx;
    int i, j, x;


    for (x = 0; x < Nx / 2; x++) {
        for (i = 0; i < dimR; ++i) {
            Re(f[x * dimR + i]) = 0;
            Im(f[x * dimR + i]) = 0;
            if (x == 2 && z == 2) {
                for (j = 0; j < qpts; j++) {
                    Im(f[x * dimR + i]) +=
                        Rw[i][j] * Kz[z] * (-2 * M_PI *
                                            (2 * (1 - Qy[j]) +
                                             (1 -
                                              Qy[j] * Qy[j]) * 4 * Qy[j]) *
                                            cos(2 * M_PI *
                                                (n * dt - dt * e[k] -
                                                 0.01))
                                            + re * (24 * Qy[j] +
                                                    2 * Kz[z] * Kz[z] *
                                                    (2 * (1 - Qy[j]) +
                                                     (1 -
                                                      Qy[j] * Qy[j]) * 4 *
                                                     Qy[j])) * sin(2 *
                                                                   M_PI *
                                                                   (n *
                                                                    dt -
                                                                    dt *
                                                                    e[k] -
                                                                    0.01)));
                    Re(f[x * dimR + i]) +=
                        Rw[i][j] * Kz[z] * Kz[z] * (2 * (1 - Qy[j]) +
                                                    (1 -
                                                     Qy[j] * Qy[j]) * 4 *
                                                    Qy[j]) * (1 -
                                                              Qy[j] *
                                                              Qy[j])
                        * sin(2 * M_PI * (n * dt - dt * e[k] - 0.01))
                        * sin(2 * M_PI * (n * dt - dt * e[k]));
                }
            }
        }
    }

    for (x = 0; x < Nx / 2; x++) {
        for (i = 0; i < dimQ; ++i) {
            Re(g[x * dimQ + i]) = 0.;
            Im(g[x * dimQ + i]) = 0.;
            if (x == 2 && z == 2) {
                for (j = 0; j < qpts; j++) {
                    Im(g[x * dimQ + i]) +=
                        Qw[i][j] * Kz[z] *
                        (-
                         (12 * Qy[j] * Qy[j] - 4 -
                          2 * Kz[z] * Kz[z] * (1 - Qy[j] * Qy[j]) * (1 -
                                                                     Qy[j]
                                                                     *
                                                                     Qy
                                                                     [j]))
                         * 2 * M_PI * cos(2 * M_PI *
                                          (n * dt - dt * e[k] - 0.01)) -
                         re * (24 -
                               2 * 2 * Kz[z] * Kz[z] * (12 * Qy[j] *
                                                        Qy[j] - 4) +
                               4 * Kz[z] * Kz[z] * Kz[z] * Kz[z] * (1 -
                                                                    Qy[j] *
                                                                    Qy[j])
                               * (1 -
                                  Qy[j] * Qy[j])) * sin(2 * M_PI * (n *
                                                                    dt -
                                                                    dt *
                                                                    e[k] -
                                                                    0.01)));
                    Re(g[x * dimQ + i]) -=
                        Qw[i][j] * (Kz[z] * Kz[z] *
                                    ((1 -
                                      Qy[j] * Qy[j]) * (-16 * Qy[j] *
                                                        Qy[j]) + 4 * (1 -
                                                                      Qy[j]
                                                                      *
                                                                      Qy
                                                                      [j])
                                     * (1 - Qy[j] * Qy[j]))
                                    -
                                    2 * Kz[z] * Kz[z] * (-Kz[z] * Kz[z] *
                                                         (1 -
                                                          Qy[j] * Qy[j]) *
                                                         (1 -
                                                          Qy[j] * Qy[j]) *
                                                         (1 -
                                                          Qy[j] * Qy[j])
                                                         + 2 * Qy[j] * (1 -
                                                                        Qy
                                                                        [j]
                                                                        +
                                                                        (1
                                                                         -
                                                                         Qy
                                                                         [j]
                                                                         *
                                                                         Qy
                                                                         [j])
                                                                        *
                                                                        4 *
                                                                        Qy
                                                                        [j])))
                        * sin(2 * M_PI * (n * dt - dt * e[k] - 0.01))
                        * sin(2 * M_PI * (n * dt - dt * e[k]));
                }
            }
        }
    }
}


/**************forces used for manufacture solution 4 ***************************/
void adjforce0_4(int n, int k, int flag, mcomplex * f, mcomplex * g)
{
    extern double *Qy, **Rw;
    extern double dt, re;
    extern int dimR, qpts;

    int i, j;
    for (i = 0; i < dimR; ++i) {
        Re(f[i]) = 0;
        Im(f[i]) = 0;
        Re(g[i]) = 0;
        Im(g[i]) = 0;
        if (flag == 0) {
            for (j = 0; j < qpts; ++j) {
                Re(f[i]) += Rw[i][j] * 2 * re;
                //*sin(2*M_PI*(n*dt-dt*e[k] ));
            }
        }
    }
}

void adjforce_4(int n, int k, int z, mcomplex * f, mcomplex * g)
{
    extern double *Qy, **Rw, *Kz, **Qw;
    extern double dt, re;
    extern int dimR, qpts, dimQ;
    static double e[3] = { 1. / 3., 1. / 2, 1. };
    extern int Nz, Nx;
    int i, j, x;


    for (x = 0; x < Nx / 2; x++) {
        for (i = 0; i < dimR; ++i) {
            Re(f[x * dimR + i]) = 0;
            Im(f[x * dimR + i]) = 0;
            if (x == 1 && z == 1) {
                for (j = 0; j < qpts; j++) {
                    Im(f[x * dimR + i]) +=
                        Rw[i][j] * 2 * Kz[z] * (-(1 - Qy[j] * Qy[j]) *
                                                cos(2 * M_PI *
                                                    (n * dt - dt * e[k] -
                                                     0.01)) * 2 * M_PI +
                                                re * (2 +
                                                      2 * Kz[z] * Kz[z] *
                                                      (1 -
                                                       Qy[j] * Qy[j])) *
                                                sin(2 * M_PI *
                                                    (n * dt - dt * e[k] -
                                                     0.01)));
                    Re(f[x * dimR + i]) +=
                        Rw[i][j] * 2 * Kz[z] * Kz[z] * (1 -
                                                        Qy[j] * Qy[j]) *
                        (1 -
                         Qy[j] * Qy[j]) * sin(2 * M_PI * (n * dt -
                                                          dt * e[k] -
                                                          0.01));
                    // *sin(2*M_PI*(n*dt-dt*e[k])) ;
                }
            }
        }
    }

    for (x = 0; x < Nx / 2; x++) {
        for (i = 0; i < dimQ; ++i) {
            Re(g[x * dimQ + i]) = 0.;
            Im(g[x * dimQ + i]) = 0.;
            if (x == 1 && z == 1) {
                for (j = 0; j < qpts; j++) {
                    Re(g[x * dimQ + i]) +=
                        Qw[i][j] * Kz[z] * Kz[z] * 2 * 2 * Qy[j] * (1 -
                                                                    Qy[j] *
                                                                    Qy[j])
                        * sin(2 * M_PI * (n * dt - dt * e[k] - 0.01));
                    //  *sin(2*M_PI*(n*dt-dt*e[k]));
                }
            }
        }
    }
}

void increadjforce0_4(int n, int k, int flag, mcomplex * f, mcomplex * g)
{
    extern double *Qy, **Rw;
    extern double dt, re;
    extern int dimR, qpts;
    int i;
    for (i = 0; i < dimR; ++i) {
        Re(f[i]) = 0;
        Im(f[i]) = 0;
        Re(g[i]) = 0;
        Im(g[i]) = 0;
    }

}


void increadjforce_4(int n, int k, int z, mcomplex * f, mcomplex * g)
{
    extern double *Qy, **Rw, *Kz, **Qw;
    extern double dt, re;
    extern int dimR, qpts, dimQ;
    static double e[3] = { 1. / 3., 1. / 2, 1. };
    extern int Nz, Nx;
    int i, j, x;


    for (x = 0; x < Nx / 2; x++) {
        for (i = 0; i < dimR; ++i) {
            Re(f[x * dimR + i]) = 0;
            Im(f[x * dimR + i]) = 0;
            if (x == 1 && z == 1) {
                for (j = 0; j < qpts; j++) {
                    Im(f[x * dimR + i]) +=
                        Rw[i][j] * Kz[z] * (-2 * M_PI * 2 * (1 - Qy[j]) *
                                            cos(2 * M_PI *
                                                (n * dt - dt * e[k] -
                                                 0.01))
                                            +
                                            re * (2 * Kz[z] * Kz[z] *
                                                  (2 * (1 - Qy[j]))) *
                                            sin(2 * M_PI *
                                                (n * dt - dt * e[k] -
                                                 0.01)));
                    //  Re(f[x*dimR+i])+= Rw[i][j]*2*Kz[z]*Kz[z]*(1-Qy[j])*(1-Qy[j]*Qy[j])*sin(2*M_PI*(n*dt-dt*e[k]-0.01));
                    Re(f[x * dimR + i]) +=
                        Rw[i][j] * 2 * 2 * Kz[z] * Kz[z] * (1 -
                                                            Qy[j]) * (1 -
                                                                      Qy[j]
                                                                      *
                                                                      Qy
                                                                      [j])
                        * sin(2 * M_PI * (n * dt - dt * e[k] - 0.01));
                    //*sin(2*M_PI*(n*dt-dt*e[k]));
                }
            }
        }
    }

    for (x = 0; x < Nx / 2; x++) {
        for (i = 0; i < dimQ; ++i) {
            Re(g[x * dimQ + i]) = 0.;
            Im(g[x * dimQ + i]) = 0.;
            if (x == 1 && z == 1) {
                for (j = 0; j < qpts; j++) {
                    Re(g[x * dimQ + i]) +=
                        Qw[i][j] * 2 * Kz[z] * Kz[z] * (2 * Qy[j] *
                                                        (1 - Qy[j]) + 1 -
                                                        Qy[j] * Qy[j]) *
                        sin(2 * M_PI * (n * dt - dt * e[k] - 0.01));
                    //*sin(2*M_PI*(n*dt-dt*e[k])); 
                }
            }
        }
    }
}


/**************forces used for manufacture solution 5 ***************************/
void adjforce0_5(int n, int k, int flag, mcomplex * f, mcomplex * g)
{
    extern double *Qy, **Rw;
    extern double dt, re;
    extern int dimR, qpts;
    static double e[3] = { 1. / 3., 1. / 2, 1. };
    int i, j;
    for (i = 0; i < dimR; ++i) {
        Re(f[i]) = 0;
        Im(f[i]) = 0;
        Re(g[i]) = 0;
        Im(g[i]) = 0;
        if (flag == 0) {
            for (j = 0; j < qpts; ++j) {
                Re(f[i]) +=
                    Rw[i][j] * (-2 * M_PI * (1 - Qy[j] * Qy[j]) *
                                cos(2 * M_PI *
                                    (n * dt - dt * e[k] - 0.01)))
                    +
                    Rw[i][j] * 2 * re * sin(2 * M_PI *
                                            (n * dt - dt * e[k] - 0.01)) +
                    Rw[i][j] * re * 2;
                Re(g[i]) -=
                    Rw[i][j] * (-2 * M_PI * (1 - Qy[j] * Qy[j]) *
                                cos(2 * M_PI *
                                    (n * dt - dt * e[k] - 0.01)))
                    +
                    Rw[i][j] * 2 * re * sin(2 * M_PI *
                                            (n * dt - dt * e[k] - 0.01));

            }
        }
    }
}

void increadjforce0_5(int n, int k, int flag, mcomplex * f, mcomplex * g)
{
    extern double *Qy, **Rw;
    extern double dt, re;
    extern int dimR, qpts;
    static double e[3] = { 1. / 3., 1. / 2, 1. };
    int i, j;
    for (i = 0; i < dimR; ++i) {
        Re(f[i]) = 0;
        Im(f[i]) = 0;
        Re(g[i]) = 0;
        Im(g[i]) = 0;
        for (j = 0; j < qpts; ++j) {
            Re(f[i]) +=
                -Rw[i][j] * 2 * M_PI * (1 -
                                        Qy[j]) * cos(2 * M_PI * (n * dt -
                                                                 dt *
                                                                 e[k] -
                                                                 0.01));
            Re(g[i]) +=
                Rw[i][j] * 2 * M_PI * (1 -
                                       Qy[j]) * cos(2 * M_PI * (n * dt -
                                                                dt * e[k] -
                                                                0.01));
        }
    }

}

/**************forces used for manufacture solution 6 ***************************/
void adjforce0_6(int n, int k, int flag, mcomplex * f, mcomplex * g)
{
    extern double *Qy, **Rw;
    extern double dt, re;
    extern int dimR, qpts;
    static double e[3] = { 1. / 3., 1. / 2, 1. };

    int i, j;
    for (i = 0; i < dimR; ++i) {
        Re(f[i]) = 0;
        Im(f[i]) = 0;
        Re(g[i]) = 0;
        Im(g[i]) = 0;
        if (flag == 0) {
            for (j = 0; j < qpts; ++j) {
                Re(f[i]) +=
                    Rw[i][j] * (-(1 - Qy[j] * Qy[j]) * 2 * M_PI *
                                cos(2 * M_PI *
                                    (n * dt - dt * e[k] - 0.01)) +
                                2 * re * sin(2 * M_PI *
                                             (n * dt - dt * e[k] - 0.01)) +
                                re * 2);
                Re(g[i]) +=
                    -Rw[i][j] * (-(1 - Qy[j] * Qy[j]) * 2 * M_PI *
                                 cos(2 * M_PI *
                                     (n * dt - dt * e[k] - 0.01)) +
                                 2 * re * sin(2 * M_PI *
                                              (n * dt - dt * e[k] -
                                               0.01)) + re * 2);
            }
        }
    }
}


void increadjforce_6(int n, int k, int z, mcomplex * f, mcomplex * g)
{
    extern double *Qy, **Rw, *Kz, **Qw;
    extern double dt, re;
    extern int dimR, qpts, dimQ;
    static double e[3] = { 1. / 3., 1. / 2, 1. };
    extern int Nz, Nx;
    int i, j, x;


    for (x = 0; x < Nx / 2; x++) {
        for (i = 0; i < dimR; ++i) {
            Re(f[x * dimR + i]) = 0;
            Im(f[x * dimR + i]) = 0;
            if (x == 1 && z == 1) {
                for (j = 0; j < qpts; j++) {
                    Im(f[x * dimR + i]) +=
                        Rw[i][j] * Kz[z] * (-2 * M_PI *
                                            (2 * (1 - Qy[j]) +
                                             (1 -
                                              Qy[j] * Qy[j]) * 4 * Qy[j]) *
                                            cos(2 * M_PI *
                                                (n * dt - dt * e[k] -
                                                 0.01))
                                            + re * (24 * Qy[j] +
                                                    2 * Kz[z] * Kz[z] *
                                                    (2 * (1 - Qy[j]) +
                                                     (1 -
                                                      Qy[j] * Qy[j]) * 4 *
                                                     Qy[j])) * sin(2 *
                                                                   M_PI *
                                                                   (n *
                                                                    dt -
                                                                    dt *
                                                                    e[k] -
                                                                    0.01))
                                            +
                                            re * 2 * Kz[z] * Kz[z] * 2 *
                                            (1 - Qy[j]));
                }
            }
        }
    }

    for (x = 0; x < Nx / 2; x++) {
        for (i = 0; i < dimQ; ++i) {
            Re(g[x * dimQ + i]) = 0.;
            Im(g[x * dimQ + i]) = 0.;
            if (x == 1 && z == 1) {
                for (j = 0; j < qpts; j++) {
                    Im(g[x * dimQ + i]) +=
                        Qw[i][j] * Kz[z] *
                        (-
                         (12 * Qy[j] * Qy[j] - 4 -
                          2 * Kz[z] * Kz[z] * (1 - Qy[j] * Qy[j]) * (1 -
                                                                     Qy[j]
                                                                     *
                                                                     Qy
                                                                     [j]))
                         * 2 * M_PI * cos(2 * M_PI *
                                          (n * dt - dt * e[k] - 0.01)) -
                         re * (24 -
                               2 * 2 * Kz[z] * Kz[z] * (12 * Qy[j] *
                                                        Qy[j] - 4) +
                               4 * Kz[z] * Kz[z] * Kz[z] * Kz[z] * (1 -
                                                                    Qy[j] *
                                                                    Qy[j])
                               * (1 -
                                  Qy[j] * Qy[j])) * sin(2 * M_PI * (n *
                                                                    dt -
                                                                    dt *
                                                                    e[k] -
                                                                    0.01)));
                    Re(g[x * dimQ + i]) -=
                        Qw[i][j] * 16 * Kz[z] * Kz[z] * (Qy[j] * Qy[j] *
                                                         Qy[j] * Qy[j] -
                                                         Qy[j] * Qy[j]) *
                        sin(2 * M_PI * (n * dt - dt * e[k] - 0.01));
                }
            }
        }
    }
}

/**************forces used for manufacture solution 7 ***************************/
void force0_7(int n, int k, int flag, mcomplex * f, mcomplex * g)
{
    extern double *Qy, **Rw;
    extern double dt, mpg, re;
    extern int dimR, qpts;
    static double e[3] = { 0., 1. / 2., 2. / 3. };
    extern int Nz, Nx;
    int i, j;

    for (i = 0; i < dimR; ++i) {
        Re(f[i]) = 0;
        Re(g[i]) = 0;
        Im(f[i]) = 0;
        Im(g[i]) = 0;
        for (j = 0; j < qpts; ++j) {
            Re(f[i]) +=
                Rw[i][j] * ((1 - Qy[j] * Qy[j]) *
                            cos(2 * M_PI * (n * dt + e[k] * dt)) * 2 *
                            M_PI +
                            2. * re * sin(2 * M_PI *
                                          (n * dt + dt * e[k])) + mpg);
            Re(g[i]) -=
                Rw[i][j] * ((1 - Qy[j] * Qy[j]) *
                            cos(2 * M_PI * (n * dt + e[k] * dt)) * 2 *
                            M_PI +
                            2. * re * sin(2 * M_PI *
                                          (n * dt + dt * e[k])));
        }
    }

}




void increforce_7(int n, int k, int z, mcomplex * f, mcomplex * g)
{
    extern double *Qy, **Rw, **Qw, *Kz;
    extern double dt, re;
    extern int dimR, qpts, dimQ;
    static double e[3] = { 0., 1. / 2., 2. / 3 };
    extern int Nz, Nx;
    int i, j, x;

    for (x = 0; x < Nx / 2; x++) {
        for (i = 0; i < dimR; ++i) {
            Re(f[x * dimR + i]) = 0;
            Im(f[x * dimR + i]) = 0;
            if (x == 1 && z == 1) {
                for (j = 0; j < qpts; j++) {
                    Im(f[x * dimR + i]) +=
                        Rw[i][j] * Kz[z] * (2 * M_PI *
                                            (2 * (1 - Qy[j]) +
                                             (1 -
                                              Qy[j] * Qy[j]) * 4 * Qy[j]) *
                                            cos(2 * M_PI *
                                                (n * dt + dt * e[k]))
                                            + re * (24 * Qy[j] +
                                                    2 * Kz[z] * Kz[z] *
                                                    (2 * (1 - Qy[j]) +
                                                     (1 -
                                                      Qy[j] * Qy[j]) * 4 *
                                                     Qy[j])) * sin(2 *
                                                                   M_PI *
                                                                   (n *
                                                                    dt +
                                                                    dt *
                                                                    e
                                                                    [k])));
                    Re(f[x * dimR + i]) +=
                        Rw[i][j] * Kz[z] * Kz[z] * (1 -
                                                    Qy[j] * Qy[j]) * (1 -
                                                                      Qy[j]
                                                                      *
                                                                      Qy
                                                                      [j])
                        * 4 * Qy[j] * sin(2 * M_PI *
                                          (n * dt +
                                           dt * e[k])) * sin(2 * M_PI *
                                                             (n * dt +
                                                              dt * e[k]));
                }
            }
        }
    }
    for (x = 0; x < Nx / 2; x++) {
        for (i = 0; i < dimQ; ++i) {
            Re(g[x * dimQ + i]) = 0.;
            Im(g[x * dimQ + i]) = 0.;
            if (x == 1 && z == 1) {
                for (j = 0; j < qpts; j++) {
                    Im(g[x * dimQ + i]) +=
                        Qw[i][j] * Kz[z] *
                        ((12 * Qy[j] * Qy[j] - 4 -
                          2 * Kz[z] * Kz[z] * (1 - Qy[j] * Qy[j]) * (1 -
                                                                     Qy[j]
                                                                     *
                                                                     Qy
                                                                     [j]))
                         * 2 * M_PI * cos(2 * M_PI *
                                          (n * dt + dt * e[k])) -
                         re * (24 -
                               2 * 2 * Kz[z] * Kz[z] * (12 * Qy[j] *
                                                        Qy[j] - 4) +
                               4 * Kz[z] * Kz[z] * Kz[z] * Kz[z] * (1 -
                                                                    Qy[j] *
                                                                    Qy[j])
                               * (1 -
                                  Qy[j] * Qy[j])) * sin(2 * M_PI * (n *
                                                                    dt +
                                                                    dt *
                                                                    e
                                                                    [k])));
                }
            }
        }
    }
}

void adjforce0_7(int n, int k, int flag, mcomplex * f, mcomplex * g)
{
    extern double *Qy, **Rw;
    extern double dt, re;
    extern int dimR, qpts;
    static double e[3] = { 1. / 3., 1. / 2, 1. };

    int i, j;
    for (i = 0; i < dimR; ++i) {
        Re(f[i]) = 0;
        Im(f[i]) = 0;
        Re(g[i]) = 0;
        Im(g[i]) = 0;
        if (flag == 0) {
            for (j = 0; j < qpts; ++j) {
                Re(f[i]) +=
                    Rw[i][j] * (-(1 - Qy[j] * Qy[j]) * 2 * M_PI *
                                cos(2 * M_PI * (n * dt - dt * e[k] - FT)) +
                                2 * re * sin(2 * M_PI *
                                             (n * dt - dt * e[k] - FT)) +
                                re * 2 * sin(2 * M_PI *
                                             (n * dt - dt * e[k])));
                Re(g[i]) +=
                    -Rw[i][j] * (-(1 - Qy[j] * Qy[j]) * 2 * M_PI *
                                 cos(2 * M_PI *
                                     (n * dt - dt * e[k] - FT)) +
                                 2 * re * sin(2 * M_PI *
                                              (n * dt - dt * e[k] - FT)) +
                                 re * 2 * sin(2 * M_PI *
                                              (n * dt - dt * e[k])));
            }
        }
    }
}


void increadjforce_7(int n, int k, int z, mcomplex * f, mcomplex * g)
{
    extern double *Qy, **Rw, *Kz, **Qw;
    extern double dt, re;
    extern int dimR, qpts, dimQ;
    static double e[3] = { 1. / 3., 1. / 2, 1. };
    extern int Nz, Nx;
    int i, j, x;


    for (x = 0; x < Nx / 2; x++) {
        for (i = 0; i < dimR; ++i) {
            Re(f[x * dimR + i]) = 0;
            Im(f[x * dimR + i]) = 0;
            if (x == 1 && z == 1) {
                for (j = 0; j < qpts; j++) {
                    Im(f[x * dimR + i]) +=
                        Rw[i][j] * Kz[z] * (-2 * M_PI *
                                            (2 * (1 - Qy[j]) +
                                             (1 -
                                              Qy[j] * Qy[j]) * 4 * Qy[j]) *
                                            cos(2 * M_PI *
                                                (n * dt - dt * e[k] - FT))
                                            + re * (24 * Qy[j] +
                                                    2 * Kz[z] * Kz[z] *
                                                    (2 * (1 - Qy[j]) +
                                                     (1 -
                                                      Qy[j] * Qy[j]) * 4 *
                                                     Qy[j])) * sin(2 *
                                                                   M_PI *
                                                                   (n *
                                                                    dt -
                                                                    dt *
                                                                    e[k] -
                                                                    FT))
                                            + re * (24 * Qy[j] +
                                                    2 * Kz[z] * Kz[z] *
                                                    (2 * (1 - Qy[j]) +
                                                     (1 -
                                                      Qy[j] * Qy[j]) * 4 *
                                                     Qy[j])) * sin(2 *
                                                                   M_PI *
                                                                   (n *
                                                                    dt -
                                                                    dt *
                                                                    e
                                                                    [k])));
                    Re(f[x * dimR + i]) -=
                        Rw[i][j] * Kz[z] * Kz[z] * (1 -
                                                    Qy[j] * Qy[j]) * (1 -
                                                                      Qy[j]
                                                                      *
                                                                      Qy
                                                                      [j])
                        * 4 * Qy[j] * sin(2 * M_PI *
                                          (n * dt - dt * e[k] -
                                           FT)) * sin(2 * M_PI * (n * dt -
                                                                  dt *
                                                                  e[k]));
                }
            }
        }
    }

    for (x = 0; x < Nx / 2; x++) {
        for (i = 0; i < dimQ; ++i) {
            Re(g[x * dimQ + i]) = 0.;
            Im(g[x * dimQ + i]) = 0.;
            if (x == 1 && z == 1) {
                for (j = 0; j < qpts; j++) {
                    Im(g[x * dimQ + i]) +=
                        Qw[i][j] * Kz[z] *
                        (-
                         (12 * Qy[j] * Qy[j] - 4 -
                          2 * Kz[z] * Kz[z] * (1 - Qy[j] * Qy[j]) * (1 -
                                                                     Qy[j]
                                                                     *
                                                                     Qy
                                                                     [j]))
                         * 2 * M_PI * cos(2 * M_PI *
                                          (n * dt - dt * e[k] - FT)) -
                         re * (24 -
                               2 * 2 * Kz[z] * Kz[z] * (12 * Qy[j] *
                                                        Qy[j] - 4) +
                               4 * Kz[z] * Kz[z] * Kz[z] * Kz[z] * (1 -
                                                                    Qy[j] *
                                                                    Qy[j])
                               * (1 -
                                  Qy[j] * Qy[j])) * sin(2 * M_PI * (n *
                                                                    dt -
                                                                    dt *
                                                                    e[k] -
                                                                    FT)) -
                         re * (24 -
                               2 * 2 * Kz[z] * Kz[z] * (12 * Qy[j] *
                                                        Qy[j] - 4) +
                               4 * Kz[z] * Kz[z] * Kz[z] * Kz[z] * (1 -
                                                                    Qy[j] *
                                                                    Qy[j])
                               * (1 -
                                  Qy[j] * Qy[j])) * sin(2 * M_PI * (n *
                                                                    dt -
                                                                    dt *
                                                                    e
                                                                    [k])));

                }
            }
        }
    }
}


/**************forces used for manufacture solution 8 ***************************/
void force0_8(int n, int k, int flag, mcomplex * f, mcomplex * g)
{
    extern double *Qy, **Rw;
    extern double dt, mpg, re;
    extern int dimR, qpts;
    static double e[3] = { 0., 1. / 2., 2. / 3. };
    extern int Nz, Nx;
    int i, j;

    for (i = 0; i < dimR; ++i) {
        Re(f[i]) = 0;
        Re(g[i]) = 0;
        Im(f[i]) = 0;
        Im(g[i]) = 0;
        for (j = 0; j < qpts; ++j) {
            Re(f[i]) += Rw[i][j] * mpg;
        }
    }

}

void force_8(int n, int k, int z, mcomplex * f, mcomplex * g)
{
    extern double *Qy, **Rw, **Qw, *Kz, *Kx;
    extern double dt, re;
    extern int dimR, qpts, dimQ;
    static double e[3] = { 0., 1. / 2., 2. / 3 };
    extern int Nz, Nx;
    int i, j, x;

    for (x = 0; x < Nx / 2; x++) {
        for (i = 0; i < dimR; ++i) {
            Re(f[x * dimR + i]) = 0;
            Im(f[x * dimR + i]) = 0;
            if (x == 1 && z == 1) {
                for (j = 0; j < qpts; j++) {
                    Im(f[x * dimR + i]) +=
                        Rw[i][j] * 2 * Kz[z] * sin(M_PI * Qy[j] * 2) * (2 *
                                                                        M_PI
                                                                        *
                                                                        cos
                                                                        (2
                                                                         *
                                                                         M_PI
                                                                         *
                                                                         (n
                                                                          *
                                                                          dt
                                                                          +
                                                                          dt
                                                                          *
                                                                          e
                                                                          [k]))
                                                                        +
                                                                        re
                                                                        *
                                                                        (4
                                                                         *
                                                                         M_PI
                                                                         *
                                                                         M_PI
                                                                         +
                                                                         Kz
                                                                         [z]
                                                                         *
                                                                         Kz
                                                                         [z]
                                                                         +
                                                                         Kx
                                                                         [x]
                                                                         *
                                                                         Kx
                                                                         [x])
                                                                        *
                                                                        sin
                                                                        (2
                                                                         *
                                                                         M_PI
                                                                         *
                                                                         (n
                                                                          *
                                                                          dt
                                                                          +
                                                                          dt
                                                                          *
                                                                          e
                                                                          [k])));
                }
            }
        }
    }
    for (x = 0; x < Nx / 2; x++) {
        for (i = 0; i < dimQ; ++i) {
            Re(g[x * dimQ + i]) = 0.;
            Im(g[x * dimQ + i]) = 0.;
        }
    }
}

void increforce0_8(int n, int k, int flag, mcomplex * f, mcomplex * g)
{
    extern double *Qy, **Rw;
    extern double dt, re;
    extern int dimR, qpts;
    static double e[3] = { 0., 1. / 2., 2. / 3. };
    extern int Nz, Nx;
    int i, j;

    for (i = 0; i < dimR; ++i) {
        Re(f[i]) = 0;
        Re(g[i]) = 0;
        Im(f[i]) = 0;
        Im(g[i]) = 0;
        for (j = 0; j < qpts; ++j) {
            Re(f[i]) +=
                Rw[i][j] * 4 * M_PI * sin(M_PI * Qy[j] / 4 +
                                          3 * M_PI / 4) * (2 * M_PI *
                                                           cos(2 * M_PI *
                                                               (n * dt +
                                                                dt *
                                                                e[k])) +
                                                           re * M_PI *
                                                           M_PI / 16 *
                                                           sin(2 * M_PI *
                                                               (n * dt +
                                                                dt *
                                                                e[k])));
            Re(g[i]) -=
                Rw[i][j] * 4 * M_PI * sin(M_PI * Qy[j] / 4 +
                                          3 * M_PI / 4) * (2 * M_PI *
                                                           cos(2 * M_PI *
                                                               (n * dt +
                                                                dt *
                                                                e[k])) +
                                                           re * M_PI *
                                                           M_PI / 16 *
                                                           sin(2 * M_PI *
                                                               (n * dt +
                                                                dt *
                                                                e[k])));
        }
    }

}


void increforce_8(int n, int k, int z, mcomplex * f, mcomplex * g)
{
    extern double *Qy, **Rw, **Qw, *Kz, *Kx;
    extern double dt, re;
    extern int dimR, qpts, dimQ;
    static double e[3] = { 0., 1. / 2., 2. / 3 };
    extern int Nz, Nx;
    int i, j, x;

    for (x = 0; x < Nx / 2; x++) {
        for (i = 0; i < dimR; ++i) {
            Re(f[x * dimR + i]) = 0;
            Im(f[x * dimR + i]) = 0;
            if (x == 2 && z == 2) {
                for (j = 0; j < qpts; j++) {
                    Im(f[x * dimR + i]) +=
                        Rw[i][j] * Kz[z] * 4 * M_PI * sin(M_PI * Qy[j] /
                                                          4 +
                                                          3. * M_PI / 4)
                        * (2 * M_PI *
                           cos(2 * M_PI * (n * dt + dt * e[k])) +
                           re * (M_PI * M_PI / 16 + Kz[z] * Kz[z] +
                                 Kx[x] * Kx[x]) * sin(2 * M_PI * (n * dt +
                                                                  dt *
                                                                  e[k])));
                }
            }
        }
    }
    for (x = 0; x < Nx / 2; x++) {
        for (i = 0; i < dimQ; ++i) {
            Re(g[x * dimQ + i]) = 0.;
            Im(g[x * dimQ + i]) = 0.;
        }
    }
}

void adjforce_8(int n, int k, int z, mcomplex * f, mcomplex * g)
{
    extern double *Qy, **Rw, *Kz, *Kx;
    extern double dt, re;
    extern int dimR, qpts, dimQ;
    static double e[3] = { 1. / 3., 1. / 2, 1. };
    extern int Nz, Nx;
    int i, j, x;

    for (x = 0; x < Nx / 2; x++) {
        for (i = 0; i < dimR; ++i) {
            Re(f[x * dimR + i]) = 0;
            Im(f[x * dimR + i]) = 0;
            if (x == 1 && z == 1) {
                for (j = 0; j < qpts; j++) {
                    Im(f[x * dimR + i]) +=
                        Rw[i][j] * 2 * Kz[z] * sin(2 * M_PI * Qy[j]) *
                        (-2 * M_PI *
                         cos(2 * M_PI * (n * dt - dt * e[k] - FT))
                         + re * (4 * M_PI * M_PI + Kz[z] * Kz[z] +
                                 Kx[x] * Kx[x])
                         * (sin(2 * M_PI * (n * dt - dt * e[k])) +
                            sin(2 * M_PI * (n * dt - dt * e[k] - FT))));
                }
            }
        }
    }
    for (x = 0; x < Nx / 2; x++) {
        for (i = 0; i < dimQ; ++i) {
            Re(g[x * dimQ + i]) = 0.;
            Im(g[x * dimQ + i]) = 0.;
        }
    }
}

void increadjforce0_8(int n, int k, int flag, mcomplex * f, mcomplex * g)
{
    extern double *Qy, **Rw;
    extern double dt, re;
    extern int dimR, qpts;
    static double e[3] = { 1. / 3., 1. / 2, 1. };
    int i, j;
    for (i = 0; i < dimR; ++i) {
        Re(f[i]) = 0;
        Im(f[i]) = 0;
        Re(g[i]) = 0;
        Im(g[i]) = 0;
        for (j = 0; j < qpts; ++j) {
            Re(f[i]) +=
                Rw[i][j] * 4 * M_PI * sin(M_PI * Qy[j] / 4 + 3. * M_PI / 4)
                * (-2 * M_PI * cos(2 * M_PI * (n * dt - dt * e[k] - FT)) +
                   re * M_PI * M_PI / 16 *
                   (sin(2 * M_PI * (n * dt - dt * e[k] - FT)) +
                    sin(2 * M_PI * (n * dt - dt * e[k]))));
            Re(g[i]) -=
                Rw[i][j] * 4 * M_PI * sin(M_PI * Qy[j] / 4 + 3. * M_PI / 4)
                * (-2 * M_PI * cos(2 * M_PI * (n * dt - dt * e[k] - FT)) +
                   re * M_PI * M_PI / 16 *
                   (sin(2 * M_PI * (n * dt - dt * e[k] - FT)) +
                    sin(2 * M_PI * (n * dt - dt * e[k]))));
        }
    }

}

void increadjforce_8(int n, int k, int z, mcomplex * f, mcomplex * g)
{
    extern double *Qy, **Rw, *Kz, *Kx, **Qw;
    extern double dt, re;
    extern int dimR, qpts, dimQ;
    static double e[3] = { 1. / 3., 1. / 2, 1. };
    extern int Nz, Nx;
    int i, j, x;


    for (x = 0; x < Nx / 2; x++) {
        for (i = 0; i < dimR; ++i) {
            Re(f[x * dimR + i]) = 0;
            Im(f[x * dimR + i]) = 0;
            if (x == 2 && z == 2) {
                for (j = 0; j < qpts; j++) {
                    Im(f[x * dimR + i]) +=
                        Rw[i][j] * 4 * M_PI * Kz[z] * sin(M_PI * Qy[j] /
                                                          4 +
                                                          3. * M_PI / 4)
                        * (-2 * M_PI *
                           cos(2 * M_PI * (n * dt - dt * e[k] - FT)) +
                           re * (M_PI * M_PI / 16 + Kz[z] * Kz[z] +
                                 Kx[x] * Kx[x]) * (sin(2 * M_PI *
                                                       (n * dt -
                                                        dt * e[k] - FT)) +
                                                   sin(2 * M_PI *
                                                       (n * dt -
                                                        dt * e[k]))));

                }
            }
        }
    }

    for (x = 0; x < Nx / 2; x++) {
        for (i = 0; i < dimQ; ++i) {
            Re(g[x * dimQ + i]) = 0.;
            Im(g[x * dimQ + i]) = 0.;
        }
    }
}

/**************forces used for manufacture solution 10 ***************************/
void force0_10(int n, int k, int flag, mcomplex * f, mcomplex * g)
{
    extern double *Qy, **Rw;
    extern double dt, re;
    extern int dimR, qpts;
    static double e[3] = { 0., 1. / 2., 2. / 3. };
    extern int Nz, Nx;
    int i, j;

    for (i = 0; i < dimR; ++i) {
        Re(f[i]) = 0;
        Re(g[i]) = 0;
        Im(f[i]) = 0;
        Im(g[i]) = 0;
        for (j = 0; j < qpts; ++j) {
            Re(f[i]) +=
                Rw[i][j] *
                ((1 - Qy[j] * Qy[j] +
                  (1 - Qy[j] * Qy[j]) * (1 -
                                         Qy[j] * Qy[j]) * (-5. / 4)) *
                 cos(2 * M_PI * (n * dt + e[k] * dt)) * 2 * M_PI -
                 re * (-2 * sin(2 * M_PI * (n * dt + dt * e[k])) +
                       (12 * Qy[j] * Qy[j] - 4) * (15. / 16 -
                                                   5. / 4 * sin(2 * M_PI *
                                                                (n * dt +
                                                                 dt *
                                                                 e[k])))));
        }
    }

}

void increforce0_10(int n, int k, int flag, mcomplex * f, mcomplex * g)
{
    extern double *Qy, **Rw;
    extern double dt, re;
    extern int dimR, qpts;
    static double e[3] = { 0., 1. / 2., 2. / 3. };
    extern int Nz, Nx;
    int i, j;

    for (i = 0; i < dimR; ++i) {
        Re(f[i]) = 0;
        Re(g[i]) = 0;
        Im(f[i]) = 0;
        Im(g[i]) = 0;
        for (j = 0; j < qpts; ++j) {
            Re(f[i]) +=
                Rw[i][j] * ((1 - Qy[j] - (1 - Qy[j] * Qy[j]) * 3. / 2) *
                            cos(2 * M_PI * (n * dt + e[k] * dt)) * 2 *
                            M_PI -
                            re * 3 * sin(2 * M_PI * (n * dt + dt * e[k])));
        }
    }

}

void adjforce0_10(int n, int k, int flag, mcomplex * f, mcomplex * g)
{
    extern double *Qy, **Rw;
    extern double dt, re;
    extern int dimR, qpts;
    static double e[3] = { 1. / 3., 1. / 2, 1. };
    extern int Nz, Nx;
    int i, j;

    for (i = 0; i < dimR; ++i) {
        Re(f[i]) = 0;
        Re(g[i]) = 0;
        Im(f[i]) = 0;
        Im(g[i]) = 0;
        for (j = 0; j < qpts; ++j) {
            Re(f[i]) +=
                Rw[i][j] *
                (-
                 (1 - Qy[j] * Qy[j] -
                  5. / 4. * (1 - Qy[j] * Qy[j]) * (1 -
                                                   Qy[j] * Qy[j])) * 2 *
                 M_PI * cos(2 * M_PI * (n * dt - e[k] * dt - FT))
- re * (-2 -
        5 * (3 * Qy[j] * Qy[j] - 1)) * sin(2 * M_PI * (n * dt - dt * e[k] -
                                                       FT))
- re * (-2 * sin(2 * M_PI * (n * dt - dt * e[k])) +
        (12 * Qy[j] * Qy[j] - 4) * (15. / 16 -
                                    5. / 4 * sin(2 * M_PI *
                                                 (n * dt - dt * e[k])))));
        }
    }

}

void increadjforce0_10(int n, int k, int flag, mcomplex * f, mcomplex * g)
{
    extern double *Qy, **Rw;
    extern double dt, re;
    extern int dimR, qpts;
    static double e[3] = { 1. / 3., 1. / 2, 1. };
    extern int Nz, Nx;
    int i, j;

    for (i = 0; i < dimR; ++i) {
        Re(f[i]) = 0;
        Re(g[i]) = 0;
        Im(f[i]) = 0;
        Im(g[i]) = 0;
        for (j = 0; j < qpts; ++j) {
            Re(f[i]) +=
                Rw[i][j] * (-(1 - Qy[j] - (1 - Qy[j] * Qy[j]) * 3. / 2) *
                            cos(2 * M_PI * (n * dt - e[k] * dt - FT)) * 2 *
                            M_PI -
                            re * 3 * sin(2 * M_PI *
                                         (n * dt - dt * e[k] - FT))
                            -
                            re * 3 * sin(2 * M_PI * (n * dt - dt * e[k])));
        }
    }

}

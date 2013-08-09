/* update by shan: add the part to compute H=u x I omega + Iu x omega */

/* Compute H = u x omega for forward, H=u x Iomega+Iu x omega for incremental state */

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include "fftw.h"
#include "rfftw.h"
#include "minChnl.h"

int pass1(int dctr, int n)
{
    /* External Variables */
    extern int qpts, Nx, Nz;
    extern double dt, cfl1, cfl3, *cfl2;
    extern double *Kx, *Kz;
    extern fftw_complex ***CT, ***ICT;  /* 6-by-(3Nz/2)-by-(3*Nx/4+1) */
    extern mcomplex ****U, ****IU;
    extern fftw_plan pf1, pf2;
    extern rfftwnd_plan pr1, pr2;

    extern fftw_plan Ipf1, Ipf2;

    /* Local variables */
    int x, y, z, i, j, idx;
    double uy, norm, t, CFL, Iuy;
    fftw_real *RT, *IRT;        /* real to complex transform */
    fftw_complex *fout = NULL;
    fftw_real *rout = NULL;

    idx = (3 * Nz / 2) * (3 * Nx / 2 + 2);
    RT = (fftw_real *) CT[0][0];
    IRT = (fftw_real *) ICT[0][0];
    norm = 1.0 / ((3. * Nx / 2.) * (3. * Nz / 2.));

    for (y = 0; y < qpts; ++y) {

        /* first compute omega= \nabla x u */

        /* Clear the array CT */

        for (i = 0; i < MAXT; i++) {
            for (j = 0; j < 3 * Nz / 2; j++) {
                memset(CT[i][j], 0,
                       (3 * Nx / 4 + 1) * sizeof(fftw_complex));
            }
        }

        /* Load "first half" of array */
        for (i = 0; i < MAXU; ++i) {
            for (z = 0; z < Nz / 2; ++z) {
                memcpy(CT[i][z], U[z][i][y],
                       (Nx / 2) * sizeof(fftw_complex));
            }
        }
        for (z = 0; z < Nz / 2; ++z) {
            for (x = 0; x < (Nx / 2); ++x) {
                /* omega_x hat */
                Re(CT[WXEL][z][x]) = Re(CT[DZEL][z][x]) +
                    Kz[z] * Im(CT[YEL][z][x]);
                Im(CT[WXEL][z][x]) = Im(CT[DZEL][z][x]) -
                    Kz[z] * Re(CT[YEL][z][x]);
                /* omega_z hat */
                Re(CT[WZEL][z][x]) = -(Re(CT[DXEL][z][x]) +
                                       Kx[x] * Im(CT[YEL][z][x]));
                Im(CT[WZEL][z][x]) = -Im(CT[DXEL][z][x]) +
                    Kx[x] * Re(CT[YEL][z][x]);
                /* omega_y hat */
                Re(CT[WYEL][z][x]) = Kx[x] * Im(CT[ZEL][z][x]) -
                    Kz[z] * Im(CT[XEL][z][x]);
                Im(CT[WYEL][z][x]) = Kz[z] * Re(CT[XEL][z][x]) -
                    Kx[x] * Re(CT[ZEL][z][x]);

            }
        }

        /* Load "second half" of array */
        for (i = 0; i < MAXU; ++i) {
            for (z = Nz / 2 + 1; z < Nz; ++z) {
                memcpy(CT[i][z + Nz / 2], U[z][i][y],
                       (Nx / 2) * sizeof(fftw_complex));
            }
        }
        for (z = Nz + 1; z < (3 * Nz / 2); ++z) {
            for (x = 0; x < (Nx / 2); ++x) {
                /* omega_x hat */
                Re(CT[WXEL][z][x]) = Re(CT[DZEL][z][x]) +
                    Kz[z - Nz / 2] * Im(CT[YEL][z][x]);
                Im(CT[WXEL][z][x]) = Im(CT[DZEL][z][x]) -
                    Kz[z - Nz / 2] * Re(CT[YEL][z][x]);
                /* omega_z hat */
                Re(CT[WZEL][z][x]) = -(Re(CT[DXEL][z][x]) +
                                       Kx[x] * Im(CT[YEL][z][x]));
                Im(CT[WZEL][z][x]) = -Im(CT[DXEL][z][x]) +
                    Kx[x] * Re(CT[YEL][z][x]);
                /* omega_y hat */
                Re(CT[WYEL][z][x]) = Kx[x] * Im(CT[ZEL][z][x]) -
                    Kz[z - Nz / 2] * Im(CT[XEL][z][x]);
                Im(CT[WYEL][z][x]) = Kz[z - Nz / 2] * Re(CT[XEL][z][x]) -
                    Kx[x] * Re(CT[ZEL][z][x]);
            }
        }

        /* inverse Fourier transform */
        for (i = 0; i < MAXT; ++i) {
            /* Each column of CT[i] */
            fftw(pf1, Nx / 2, CT[i][0], 3 * Nx / 4 + 1, 1, fout, -1, -1);

            /* Each row of CT[i] */
            rfftwnd_complex_to_real(pr1, 3 * Nz / 2, CT[i][0], 1,
                                    3 * Nx / 4 + 1, rout, -1, -1);
        }

        /* Check CFL condition */
        if (dctr == 0) {
            CFL = 0.0;
            for (i = 0; i < 3 * Nz / 2; ++i) {
                for (j = 0; j < 3 * Nx / 2; ++j) {
                    t = fabs(RT[XEL * idx + (i * (3 * Nx / 2 + 2) + j)]) *
                        cfl1 +
                        fabs(RT[YEL * idx + (i * (3 * Nx / 2 + 2) + j)]) *
                        cfl2[y] +
                        fabs(RT[ZEL * idx + (i * (3 * Nx / 2 + 2) + j)]) *
                        cfl3;
                    if (t > CFL) {
                        CFL = t;
                    }
                }
            }
            CFL *= dt;
            if (CFL >= CFL_MAX) {
                printf("Pass1:  CFL[%d]:%e.  CFL condition violated.\n",
                       y, CFL);
                return (ERROR);
            }
        }


        /* now compute Iomega= \nabla x Iu */

        /* Clear the array ICT */
        memset(ICT[0][0], 0,
               MAXT * (3 * Nz / 2) * (3 * Nx / 4 +
                                      1) * sizeof(fftw_complex));

        /* Load "first half" of array */
        for (i = 0; i < MAXU; ++i) {
            for (z = 0; z < Nz / 2; ++z) {
                memcpy(ICT[i][z], IU[z][i][y],
                       (Nx / 2) * sizeof(fftw_complex));
            }
        }
        for (z = 0; z < Nz / 2; ++z) {
            for (x = 0; x < (Nx / 2); ++x) {
                /* Iomega_x hat */
                Re(ICT[WXEL][z][x]) = Re(ICT[DZEL][z][x]) +
                    Kz[z] * Im(ICT[YEL][z][x]);
                Im(ICT[WXEL][z][x]) = Im(ICT[DZEL][z][x]) -
                    Kz[z] * Re(ICT[YEL][z][x]);
                /* Iomega_z hat */
                Re(ICT[WZEL][z][x]) = -(Re(ICT[DXEL][z][x]) +
                                        Kx[x] * Im(ICT[YEL][z][x]));
                Im(ICT[WZEL][z][x]) = -Im(ICT[DXEL][z][x]) +
                    Kx[x] * Re(ICT[YEL][z][x]);
                /* Iomega_y hat */
                Re(ICT[WYEL][z][x]) = Kx[x] * Im(ICT[ZEL][z][x]) -
                    Kz[z] * Im(ICT[XEL][z][x]);
                Im(ICT[WYEL][z][x]) = Kz[z] * Re(ICT[XEL][z][x]) -
                    Kx[x] * Re(ICT[ZEL][z][x]);
            }
        }

        /* Load "second half" of array */
        for (i = 0; i < MAXU; ++i) {
            for (z = Nz / 2 + 1; z < Nz; ++z) {
                memcpy(ICT[i][z + Nz / 2], IU[z][i][y],
                       (Nx / 2) * sizeof(fftw_complex));
            }
        }
        for (z = Nz + 1; z < (3 * Nz / 2); ++z) {
            for (x = 0; x < (Nx / 2); ++x) {
                /* omega_x hat */
                Re(ICT[WXEL][z][x]) = Re(ICT[DZEL][z][x]) +
                    Kz[z - Nz / 2] * Im(ICT[YEL][z][x]);
                Im(ICT[WXEL][z][x]) = Im(ICT[DZEL][z][x]) -
                    Kz[z - Nz / 2] * Re(ICT[YEL][z][x]);
                /* omega_z hat */
                Re(ICT[WZEL][z][x]) = -(Re(ICT[DXEL][z][x]) +
                                        Kx[x] * Im(ICT[YEL][z][x]));
                Im(ICT[WZEL][z][x]) = -Im(ICT[DXEL][z][x]) +
                    Kx[x] * Re(ICT[YEL][z][x]);
                /* omega_y hat */
                Re(ICT[WYEL][z][x]) = Kx[x] * Im(ICT[ZEL][z][x]) -
                    Kz[z - Nz / 2] * Im(ICT[XEL][z][x]);
                Im(ICT[WYEL][z][x]) = Kz[z - Nz / 2] * Re(ICT[XEL][z][x]) -
                    Kx[x] * Re(ICT[ZEL][z][x]);
            }
        }

        /* inverse Fourier transform */
        for (i = 0; i < MAXT; ++i) {
            /* Each column of ICT[i] */
            fftw(Ipf1, Nx / 2, ICT[i][0], 3 * Nx / 4 + 1, 1, fout, -1, -1);

            /* Each row of ICT[i] */
            rfftwnd_complex_to_real(pr1, 3 * Nz / 2, ICT[i][0], 1,
                                    3 * Nx / 4 + 1, rout, -1, -1);
        }

        /* Compute H = u x omega for state and IH=u x Iomega+ Iu x omega for incremental */
        for (z = 0; z < (3 * Nz / 2); ++z) {
            for (x = 0; x < 3 * Nx / 2; ++x) {
                Iuy = IRT[YEL * idx + (z * (3 * Nx / 2 + 2) + x)];
                uy = RT[YEL * idx + (z * (3 * Nx / 2 + 2) + x)];

                /*first compute IH= u x Iomega + Iu x omega */
                IRT[HXEL * idx + (z * (3 * Nx / 2 + 2) + x)] =
                    IRT[YEL * idx +
                        (z * (3 * Nx / 2 + 2) + x)] * RT[WZEL * idx +
                                                         (z *
                                                          (3 * Nx / 2 +
                                                           2) + x)] -
                    IRT[ZEL * idx +
                        (z * (3 * Nx / 2 + 2) + x)] * RT[WYEL * idx +
                                                         (z *
                                                          (3 * Nx / 2 +
                                                           2) + x)] +
                    RT[YEL * idx +
                       (z * (3 * Nx / 2 + 2) + x)] * IRT[WZEL * idx +
                                                         (z *
                                                          (3 * Nx / 2 +
                                                           2) + x)] -
                    RT[ZEL * idx +
                       (z * (3 * Nx / 2 + 2) + x)] * IRT[WYEL * idx +
                                                         (z *
                                                          (3 * Nx / 2 +
                                                           2) + x)];

                IRT[HYEL * idx + (z * (3 * Nx / 2 + 2) + x)] =
                    IRT[ZEL * idx +
                        (z * (3 * Nx / 2 + 2) + x)] * RT[WXEL * idx +
                                                         (z *
                                                          (3 * Nx / 2 +
                                                           2) + x)] -
                    IRT[XEL * idx +
                        (z * (3 * Nx / 2 + 2) + x)] * RT[WZEL * idx +
                                                         (z *
                                                          (3 * Nx / 2 +
                                                           2) + x)] +
                    RT[ZEL * idx +
                       (z * (3 * Nx / 2 + 2) + x)] * IRT[WXEL * idx +
                                                         (z *
                                                          (3 * Nx / 2 +
                                                           2) + x)] -
                    RT[XEL * idx +
                       (z * (3 * Nx / 2 + 2) + x)] * IRT[WZEL * idx +
                                                         (z *
                                                          (3 * Nx / 2 +
                                                           2) + x)];

                IRT[HZEL * idx + (z * (3 * Nx / 2 + 2) + x)] =
                    IRT[XEL * idx +
                        (z * (3 * Nx / 2 + 2) + x)] * RT[WYEL * idx +
                                                         (z *
                                                          (3 * Nx / 2 +
                                                           2) + x)] -
                    Iuy * RT[WXEL * idx + (z * (3 * Nx / 2 + 2) + x)] +
                    RT[XEL * idx +
                       (z * (3 * Nx / 2 + 2) + x)] * IRT[WYEL * idx +
                                                         (z *
                                                          (3 * Nx / 2 +
                                                           2) + x)] -
                    uy * IRT[WXEL * idx + (z * (3 * Nx / 2 + 2) + x)];


                /* next compute H=u x omega */
                RT[HXEL * idx + (z * (3 * Nx / 2 + 2) + x)] =
                    RT[YEL * idx +
                       (z * (3 * Nx / 2 + 2) + x)] * RT[WZEL * idx +
                                                        (z *
                                                         (3 * Nx / 2 + 2) +
                                                         x)] -
                    RT[ZEL * idx +
                       (z * (3 * Nx / 2 + 2) + x)] * RT[WYEL * idx +
                                                        (z *
                                                         (3 * Nx / 2 + 2) +
                                                         x)];
                RT[HYEL * idx + (z * (3 * Nx / 2 + 2) + x)] =
                    RT[ZEL * idx +
                       (z * (3 * Nx / 2 + 2) + x)] * RT[WXEL * idx +
                                                        (z *
                                                         (3 * Nx / 2 + 2) +
                                                         x)] -
                    RT[XEL * idx +
                       (z * (3 * Nx / 2 + 2) + x)] * RT[WZEL * idx +
                                                        (z *
                                                         (3 * Nx / 2 + 2) +
                                                         x)];
                RT[HZEL * idx + (z * (3 * Nx / 2 + 2) + x)] =
                    RT[XEL * idx +
                       (z * (3 * Nx / 2 + 2) + x)] * RT[WYEL * idx +
                                                        (z *
                                                         (3 * Nx / 2 + 2) +
                                                         x)] -
                    uy * RT[WXEL * idx + (z * (3 * Nx / 2 + 2) + x)];
            }
        }

        /* Fourier transform to get H hats and IH hats. */
        for (i = 0; i < 3; ++i) {
            /* start of H hats */
            /* Each row of RT[i] */
            rfftwnd_real_to_complex(pr2, 3 * Nz / 2, RT + (i * idx), 1,
                                    3 * Nx / 2 + 2, fout, -1, -1);

            /* Each column of CT[i] */
            fftw(pf2, Nx / 2, CT[i][0], 3 * Nx / 4 + 1, 1, fout, -1, -1);

            /* Put data back in array U */
            for (z = 0; z < Nz / 2; ++z) {
                for (x = 0; x < Nx / 2; ++x) {
                    Re(CT[i][z][x]) = norm * Re(CT[i][z][x]);
                    Im(CT[i][z][x]) = norm * Im(CT[i][z][x]);
                }

                memcpy(U[z][i][y], CT[i][z],
                       Nx / 2 * sizeof(fftw_complex));
            }
            for (z = Nz + 1; z < 3 * Nz / 2; ++z) {
                for (x = 0; x < Nx / 2; ++x) {
                    Re(CT[i][z][x]) = norm * Re(CT[i][z][x]);
                    Im(CT[i][z][x]) = norm * Im(CT[i][z][x]);
                }

                memcpy(U[z - Nz / 2][i][y], CT[i][z],
                       Nx / 2 * sizeof(fftw_complex));
            }
            /* end of H hats */

            /* start of IH hats */

            /* Each row of RT[i] */
            rfftwnd_real_to_complex(pr2, 3 * Nz / 2, IRT + (i * idx), 1,
                                    3 * Nx / 2 + 2, fout, -1, -1);

            /* Each column of CT[i] */
            fftw(Ipf2, Nx / 2, ICT[i][0], 3 * Nx / 4 + 1, 1, fout, -1, -1);

            /* Put data back in array U */
            for (z = 0; z < Nz / 2; ++z) {
                for (x = 0; x < Nx / 2; ++x) {
                    Re(ICT[i][z][x]) = norm * Re(ICT[i][z][x]);
                    Im(ICT[i][z][x]) = norm * Im(ICT[i][z][x]);
                }

                memcpy(IU[z][i][y], ICT[i][z],
                       Nx / 2 * sizeof(fftw_complex));
            }
            for (z = Nz + 1; z < 3 * Nz / 2; ++z) {
                for (x = 0; x < Nx / 2; ++x) {
                    Re(ICT[i][z][x]) = norm * Re(ICT[i][z][x]);
                    Im(ICT[i][z][x]) = norm * Im(ICT[i][z][x]);
                }

                memcpy(IU[z - Nz / 2][i][y], ICT[i][z],
                       Nx / 2 * sizeof(fftw_complex));
            }
            /*End of IH hats */
        }

    }                           /* end for y... */

    return (NO_ERR);
}                               /* end pass1 */

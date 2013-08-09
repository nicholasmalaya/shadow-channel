#include <math.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include "fftw.h"
#include "rfftw.h"
#include "minChnl.h"

int comp_hess(int n, int k)
{
    /* External Variables */
    extern int Nx, Nz;
    extern fftw_complex ***CT, ***ICT;  /* 6-by-(3Nz/2)-by-(3*Nx/4+1) */
    extern mcomplex **hess;
    extern fftw_plan pf1, pf2;
    extern fftw_plan Ipf1, Ipf2;
    extern rfftwnd_plan pr1, pr2;
    extern double re, dt;
    int x, i, z, idx;
    double norm, tmp1;
    fftw_real *RT, *IRT;        /* real to complex transform */
    fftw_complex *fout = NULL;
    fftw_real *rout = NULL;
    extern mcomplex **GUxb, **GUzb;
    extern mcomplex **GIUxb, **GIUzb;
    extern mcomplex **HUxb, **HUzb;
    extern mcomplex **HAUxb, **HAUzb;
    extern mcomplex **IUxb, **IUzb;
    extern mcomplex **IAUxb, **IAUzb;
    idx = (3 * Nz / 2) * (3 * Nx / 2 + 2);
    RT = (fftw_real *) CT[0][0];
    IRT = (fftw_real *) ICT[0][0];
    norm = 1.0 / ((3. * Nx / 2.) * (3. * Nz / 2.));

    memset(CT[0][0], 0,
           9 * (3 * Nz / 2) * (3 * Nx / 4 + 1) * sizeof(fftw_complex));
    memset(ICT[0][0], 0,
           9 * (3 * Nz / 2) * (3 * Nx / 4 + 1) * sizeof(fftw_complex));

    for (z = 0; z < Nz / 2; ++z) {
        memcpy(CT[0][z], GUxb[z], (Nx / 2) * sizeof(fftw_complex));
        memcpy(CT[1][z], GUzb[z], (Nx / 2) * sizeof(fftw_complex));

        memcpy(CT[2][z], GIUxb[z], (Nx / 2) * sizeof(fftw_complex));
        memcpy(CT[3][z], GIUzb[z], (Nx / 2) * sizeof(fftw_complex));

        memcpy(CT[4][z], IUxb[z], (Nx / 2) * sizeof(fftw_complex));
        memcpy(CT[5][z], IUzb[z], (Nx / 2) * sizeof(fftw_complex));

        memcpy(CT[6][z], IAUxb[z], (Nx / 2) * sizeof(fftw_complex));
        memcpy(CT[7][z], IAUzb[z], (Nx / 2) * sizeof(fftw_complex));

        memcpy(ICT[0][z], HUxb[z], (Nx / 2) * sizeof(fftw_complex));
        memcpy(ICT[1][z], HUzb[z], (Nx / 2) * sizeof(fftw_complex));

        memcpy(ICT[2][z], HAUxb[z], (Nx / 2) * sizeof(fftw_complex));
        memcpy(ICT[3][z], HAUzb[z], (Nx / 2) * sizeof(fftw_complex));

    }

    for (z = Nz / 2 + 1; z < Nz; ++z) {
        memcpy(CT[0][z + Nz / 2], GUxb[z],
               (Nx / 2) * sizeof(fftw_complex));
        memcpy(CT[1][z + Nz / 2], GUzb[z],
               (Nx / 2) * sizeof(fftw_complex));
        memcpy(CT[2][z + Nz / 2], GIUxb[z],
               (Nx / 2) * sizeof(fftw_complex));
        memcpy(CT[3][z + Nz / 2], GIUzb[z],
               (Nx / 2) * sizeof(fftw_complex));
        memcpy(CT[4][z + Nz / 2], IUxb[z],
               (Nx / 2) * sizeof(fftw_complex));
        memcpy(CT[5][z + Nz / 2], IUzb[z],
               (Nx / 2) * sizeof(fftw_complex));

        memcpy(CT[6][z + Nz / 2], IAUxb[z],
               (Nx / 2) * sizeof(fftw_complex));
        memcpy(CT[7][z + Nz / 2], IAUzb[z],
               (Nx / 2) * sizeof(fftw_complex));

        memcpy(ICT[0][z + Nz / 2], HUxb[z],
               (Nx / 2) * sizeof(fftw_complex));
        memcpy(ICT[1][z + Nz / 2], HUzb[z],
               (Nx / 2) * sizeof(fftw_complex));

        memcpy(ICT[2][z + Nz / 2], HAUxb[z],
               (Nx / 2) * sizeof(fftw_complex));
        memcpy(ICT[3][z + Nz / 2], HAUzb[z],
               (Nx / 2) * sizeof(fftw_complex));
    }

    Re(ICT[4][1][1]) = 1.;

    /*for( i=0; i< 8; i++){
       for (z=0; z< 3*Nz/2; z++)
       {
       for (x=0; x< 3*Nx/4+1; x++)
       {
       printf(" CT[%d][%d][%d]=%f+i%f\n", i, z, x, Re(CT[i][z][x]), Im(CT[i][z][x]));
       }
       }
       } */
    /* inverse Fourier transform */
    for (i = 0; i < 8; ++i) {
        /* Each column of CT[i] */
        fftw(pf1, Nx / 2, CT[i][0], 3 * Nx / 4 + 1, 1, fout, -1, -1);

        /* Each row of CT[i] */
        rfftwnd_complex_to_real(pr1, 3 * Nz / 2, CT[i][0], 1,
                                3 * Nx / 4 + 1, rout, -1, -1);
    }

    //      for( i=0; i< 8; i++)
    // {
    /*i=0;
       for (z=0; z< 3*Nz/2; z++)
       {
       for (x=0; x< 3*Nx/2; x++)
       {
       printf(" RT[%d][%d][%d]=%f\n", i, z, x, RT[i*idx+z*(3*(Nx/2+2)+x)]);
       }
       } */
    //}

    /* inverse Fourier transform */
    for (i = 0; i < 5; ++i) {
        /* Each column of ICT[i] */
        fftw(Ipf1, Nx / 2, ICT[i][0], 3 * Nx / 4 + 1, 1, fout, -1, -1);

        /* Each row of ICT[i] */
        rfftwnd_complex_to_real(pr1, 3 * Nz / 2, ICT[i][0], 1,
                                3 * Nx / 4 + 1, rout, -1, -1);
    }

    /* compute (dux)*(w.n) and (duz)*(w.n) */
    for (z = 0; z < (3 * Nz / 2); ++z) {
        for (x = 0; x < 3 * Nx / 2; ++x) {
            IRT[idx * 2 + (z * (3 * Nx / 2 + 2) + x)] =
                (RT[(z * (3 * Nx / 2 + 2) + x)] *
                 IRT[2 * idx + (z * (3 * Nx / 2 + 2) + x)] +
                 IRT[(z * (3 * Nx / 2 + 2) + x)] * RT[2 * idx +
                                                      (z *
                                                       (3 * Nx / 2 + 2) +
                                                       x)]) * IRT[4 * idx +
                                                                  (z *
                                                                   (3 *
                                                                    Nx /
                                                                    2 +
                                                                    2) +
                                                                   x)];
            tmp1 = IRT[idx + (z * (3 * Nx / 2 + 2) + x)];
            IRT[idx + (z * (3 * Nx / 2 + 2) + x)] =
                (RT[idx + (z * (3 * Nx / 2 + 2) + x)] *
                 IRT[3 * idx + (z * (3 * Nx / 2 + 2) + x)] + IRT[idx +
                                                                 (z *
                                                                  (3 * Nx /
                                                                   2 + 2) +
                                                                  x)] *
                 RT[3 * idx + (z * (3 * Nx / 2 + 2) + x)]) * IRT[4 * idx +
                                                                 (z *
                                                                  (3 * Nx /
                                                                   2 + 2) +
                                                                  x)];
            IRT[(z * (3 * Nx / 2 + 2) + x)] =
                (RT[(z * (3 * Nx / 2 + 2) + x)] *
                 IRT[(z * (3 * Nx / 2 + 2) + x)] + tmp1 * RT[idx +
                                                             (z *
                                                              (3 * Nx / 2 +
                                                               2) +
                                                              x)]) *
                IRT[4 * idx + (z * (3 * Nx / 2 + 2) + x)];
        }
    }

    for (z = 0; z < (3 * Nz / 2); ++z) {
        for (x = 0; x < 3 * Nx / 2; ++x) {
            RT[2 * idx + z * (3 * Nx / 2 + 2) + x] =
                RT[2 * idx + (z * (3 * Nx / 2 + 2) + x)] * RT[4 * idx +
                                                              (z *
                                                               (3 * Nx /
                                                                2 + 2) +
                                                               x)] +
                RT[3 * idx + (z * (3 * Nx / 2 + 2) + x)] * RT[5 * idx +
                                                              (z *
                                                               (3 * Nx /
                                                                2 + 2) +
                                                               x)];
            tmp1 = RT[idx + (z * (3 * Nx / 2 + 2) + x)];
            RT[idx + (z * (3 * Nx / 2 + 2) + x)] =
                RT[(z * (3 * Nx / 2 + 2) + x)] * RT[6 * idx +
                                                    (z * (3 * Nx / 2 + 2) +
                                                     x)] + RT[idx +
                                                              (z *
                                                               (3 * Nx /
                                                                2 + 2) +
                                                               x)] * RT[7 *
                                                                        idx
                                                                        +
                                                                        (z
                                                                         *
                                                                         (3
                                                                          *
                                                                          Nx
                                                                          /
                                                                          2
                                                                          +
                                                                          2)
                                                                         +
                                                                         x)];
            RT[z * (3 * Nx / 2 + 2) + x] =
                RT[(z * (3 * Nx / 2 + 2) + x)] * RT[4 * idx +
                                                    (z * (3 * Nx / 2 + 2) +
                                                     x)] +
                tmp1 * RT[5 * idx + (z * (3 * Nx / 2 + 2) + x)];
        }
    }

    /* Fourier transform to get Uxb hats and Uzb hats. */
    /* Each row of RT[i] */

    for (i = 0; i < 3; ++i) {
        rfftwnd_real_to_complex(pr2, 3 * Nz / 2, RT + (i * idx), 1,
                                3 * Nx / 2 + 2, fout, -1, -1);

        /* Each column of CT[i] */
        fftw(pf2, Nx / 2, CT[i][0], 3 * Nx / 4 + 1, 1, fout, -1, -1);

        /* constant of FFT */
        for (z = 0; z < Nz / 2; ++z) {
            for (x = 0; x < Nx / 2; ++x) {
                Re(CT[i][z][x]) = norm * Re(CT[i][z][x]);
                Im(CT[i][z][x]) = norm * Im(CT[i][z][x]);
            }
        }
        for (z = Nz + 1; z < 3 * Nz / 2; ++z) {
            for (x = 0; x < Nx / 2; ++x) {
                Re(CT[i][z][x]) = norm * Re(CT[i][z][x]);
                Im(CT[i][z][x]) = norm * Im(CT[i][z][x]);
            }
        }
    }

    for (i = 0; i < 3; ++i) {
        rfftwnd_real_to_complex(pr2, 3 * Nz / 2, IRT + (i * idx), 1,
                                3 * Nx / 2 + 2, fout, -1, -1);

        /* Each column of CT[i] */
        fftw(Ipf2, Nx / 2, ICT[i][0], 3 * Nx / 4 + 1, 1, fout, -1, -1);


        /* constant of FFT */
        for (z = 0; z < Nz / 2; ++z) {
            for (x = 0; x < Nx / 2; ++x) {
                Re(ICT[i][z][x]) = norm * Re(ICT[i][z][x]);
                Im(ICT[i][z][x]) = norm * Im(ICT[i][z][x]);
            }
        }

        for (z = Nz + 1; z < 3 * Nz / 2; ++z) {
            for (x = 0; x < Nx / 2; ++x) {
                Re(ICT[i][z][x]) = norm * Re(ICT[i][z][x]);
                Im(ICT[i][z][x]) = norm * Im(ICT[i][z][x]);
            }
        }
    }

    /*      for (z=0; z< 3*Nz/2; z++)
       {
       for(x=0; x< Nx/2; x++)
       { */

    /*z=1;
       x=1;
       printf("z=1, x=1, CT[%d][%d]=%f+i%f, %f+i%f, %f+i%f, %f, %f, %f\n", z, x, Re(CT[0][z][x]), Im(CT[0][z][x]), 
       Re(CT[1][z][x]), Im(CT[1][z][x]), Re(CT[2][z][x]), Im(CT[2][z][x]),
       -20*sin(2*M_PI*(n*dt-dt))*sin(2*M_PI*(n*dt-dt)),-20*sin(2*M_PI*(n*dt-dt-FT))*sin(2*M_PI*(n*dt-dt)),
       20*sin(2*M_PI*(n*dt-dt-FT))*sin(2*M_PI*(n*dt-dt)));

       printf("ICT: %f+i%f, %f+i%f, %f+i%f, %f, %f, %f\n",  Re(ICT[0][z][x]), Im(ICT[0][z][x]), Re(ICT[1][z][x]),
       Im(ICT[1][z][x]),Re(ICT[2][z][x]), Im(ICT[2][z][x]), 64*sin(2*M_PI*(n*dt-dt))*sin(2*M_PI*(n*dt-dt)),32*sin(2*M_PI*(n*dt-dt-FT))*sin(2*M_PI*(n*dt-dt)) ,
       32*sin(2*M_PI*(n*dt-dt))*sin(2*M_PI*(n*dt-dt))); */

    // }
    // }
    for (z = 0; z < Nz / 2; z++) {
        for (x = 0; x < Nx / 2; x++) {
            Re(hess[z][x]) +=
                re * (-Re(CT[0][z][x]) - Re(CT[1][z][x]) -
                      Re(CT[2][z][x]) - Re(ICT[0][z][x]) -
                      Re(ICT[1][z][x]) - Re(ICT[2][z][x])) * dt;
            Im(hess[z][x]) +=
                2 * re * (-Im(CT[0][z][x]) - Im(CT[1][z][x]) -
                          Im(CT[2][z][x]) - Im(ICT[0][z][x]) -
                          Im(ICT[1][z][x]) - Im(ICT[2][z][x])) * dt;
        }
    }
    for (z = Nz / 2 + 1; z < Nz; z++) {
        for (x = 0; x < Nx / 2; x++) {
            Re(hess[z][x]) +=
                re * (-Re(CT[0][z + Nz / 2][x]) -
                      Re(CT[1][z + Nz / 2][x]) - Re(CT[2][z + Nz / 2][x]) -
                      Re(ICT[0][z + Nz / 2][x]) -
                      Re(ICT[1][z + Nz / 2][x]) -
                      Re(ICT[2][z + Nz / 2][x])) * dt;
            Im(hess[z][x]) +=
                2 * re * (-Im(CT[0][z + Nz / 2][x]) -
                          Im(CT[1][z + Nz / 2][x]) -
                          Im(CT[2][z + Nz / 2][x]) -
                          Im(ICT[0][z + Nz / 2][x]) -
                          Im(ICT[1][z + Nz / 2][x]) -
                          Im(ICT[2][z + Nz / 2][x])) * dt;
        }
    }
    return (NO_ERR);
}

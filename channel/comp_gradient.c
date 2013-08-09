#include <math.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include "fftw.h"
#include "rfftw.h"
#include "minChnl.h"

int comp_gradient(int n, int k)
{
    /* External Variables */
    extern int Nx, Nz;
    extern fftw_complex ***CT;  /* 6-by-(3Nz/2)-by-(3*Nx/4+1) */
    extern mcomplex **grad;
    extern fftw_plan pf1, pf2;
    extern rfftwnd_plan pr1, pr2;
    extern double re, dt;
    int x, i, z, idx;
    double norm;
    fftw_real *RT;              /* real to complex transform */
    fftw_complex *fout = NULL;
    fftw_real *rout = NULL;
    extern mcomplex **GUxb, **GUzb;
    extern mcomplex **GIUxb, **GIUzb;
    idx = (3 * Nz / 2) * (3 * Nx / 2 + 2);
    RT = (fftw_real *) CT[0][0];
    norm = 1.0 / ((3. * Nx / 2.) * (3. * Nz / 2.));

    memset(CT[0][0], 0,
           MAXT * (3 * Nz / 2) * (3 * Nx / 4 + 1) * sizeof(fftw_complex));

    /* store Uxb hat and Uzb hat and w hat on CT for inverse FFT */
    /* for (z=0; z<Nz/2; ++z)
       {
       for( x=0; x< Nx/2; x++)
       {
       printf("AUxb[%d][%d]=%f+i%f, %f+i%f\n", z, x, Re(AUxb[z][x]), Im(AUxb[z][x]),
       Re(AUzb[z][x]), Im(AUzb[z][x]));
       printf("Uxb[%d][%d]=%f+i%f, %f+i%f\n", z, x, Re(Uxb[z][x]), Im(Uxb[z][x]),
       Re(Uzb[z][x]), Im(Uzb[z][x]));
       if(z*z+x*x >0)
       {
       Re(CT[0][z][x])=(Kz[z]*Im(AUxb[z][x])-Kx[x]*Im(AUzb[z][x]))/K2[z][x];
       Im(CT[0][z][x])=(-Kz[z]*Re(AUxb[z][x])+Kx[x]*Re(AUzb[z][x]))/K2[z][x];
       Re(CT[1][z][x])=(-Kx[x]*Im(AUxb[z][x])-Kz[z]*Im(AUzb[z][x]))/K2[z][x];
       Im(CT[1][z][x])=(Kx[x]*Re(AUxb[z][x])+Kz[z]*Re(AUzb[z][x]))/K2[z][x];

       Re(CT[2][z][x])=(Kz[z]*Im(Uxb[z][x])-Kx[x]*Im(Uzb[z][x]))/K2[z][x];
       Im(CT[2][z][x])=(-Kz[z]*Re(Uxb[z][x])+Kx[x]*Re(Uzb[z][x]))/K2[z][x];
       Re(CT[3][z][x])=(-Kx[x]*Im(Uxb[z][x])-Kz[z]*Im(Uzb[z][x]))/K2[z][x];
       Im(CT[3][z][x])=(Kx[x]*Re(Uxb[z][x])+Kz[z]*Re(Uzb[z][x]))/K2[z][x];
       }
       }
       }

       for (z=Nz/2+1; z<Nz; ++z)
       {
       for( x=0; x< Nx/2; x++)
       {
       Re(CT[0][z+Nz/2][x])=(Kz[z]*Im(AUxb[z][x])-Kx[x]*Im(AUzb[z][x]))/K2[z][x];
       Im(CT[0][z+Nz/2][x])=(-Kz[z]*Re(AUxb[z][x])+Kx[x]*Re(AUzb[z][x]))/K2[z][x];
       Re(CT[1][z+Nz/2][x])=(-Kx[x]*Im(AUxb[z][x])-Kz[z]*Im(AUzb[z][x]))/K2[z][x];
       Im(CT[1][z+Nz/2][x])=(Kx[x]*Re(AUxb[z][x])+Kz[z]*Re(AUzb[z][x]))/K2[z][x];

       Re(CT[2][z+Nz/2][x])=(Kz[z]*Im(Uxb[z][x])-Kx[x]*Im(Uzb[z][x]))/K2[z][x];
       Im(CT[2][z+Nz/2][x])=(-Kz[z]*Re(Uxb[z][x])+Kx[x]*Re(Uzb[z][x]))/K2[z][x];
       Re(CT[3][z+Nz/2][x])=(-Kx[x]*Im(Uxb[z][x])-Kz[z]*Im(Uzb[z][x]))/K2[z][x];
       Im(CT[3][z+Nz/2][x])=(Kx[x]*Re(Uxb[z][x])+Kz[z]*Re(Uzb[z][x]))/K2[z][x];
       }
       } */

    for (z = 0; z < Nz / 2; ++z) {
        memcpy(CT[0][z], GUxb[z], (Nx / 2) * sizeof(fftw_complex));
        memcpy(CT[1][z], GUzb[z], (Nx / 2) * sizeof(fftw_complex));
        memcpy(CT[2][z], GIUxb[z], (Nx / 2) * sizeof(fftw_complex));
        memcpy(CT[3][z], GIUzb[z], (Nx / 2) * sizeof(fftw_complex));
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
    }



    /* inverse Fourier transform */
    for (i = 0; i < 4; ++i) {
        /* Each column of CT[i] */
        fftw(pf1, Nx / 2, CT[i][0], 3 * Nx / 4 + 1, 1, fout, -1, -1);

        /* Each row of CT[i] */
        rfftwnd_complex_to_real(pr1, 3 * Nz / 2, CT[i][0], 1,
                                3 * Nx / 4 + 1, rout, -1, -1);
    }

    /*i=0;
       for (z=0; z< 3*Nz/2; z++)
       {
       for (x=0; x< 3*Nx/2; x++)
       {
       printf(" RT[%d][%d][%d]=%f\n", i, z, x, RT[i*idx+z*(3*(Nx/2+2)+x)]);
       }
       } */

    /* compute (dux)*(w.n) and (duz)*(w.n) */
    for (z = 0; z < (3 * Nz / 2); ++z) {
        for (x = 0; x < 3 * Nx / 2; ++x) {
            RT[idx * 2 + (z * (3 * Nx / 2 + 2) + x)] =
                RT[(z * (3 * Nx / 2 + 2) + x)] * RT[2 * idx +
                                                    (z * (3 * Nx / 2 + 2) +
                                                     x)];
            RT[idx * 3 + (z * (3 * Nx / 2 + 2) + x)] =
                RT[idx + (z * (3 * Nx / 2 + 2) + x)] * RT[3 * idx +
                                                          (z *
                                                           (3 * Nx / 2 +
                                                            2) + x)];
            RT[(z * (3 * Nx / 2 + 2) + x)] =
                RT[(z * (3 * Nx / 2 + 2) + x)] * RT[(z * (3 * Nx / 2 + 2) +
                                                     x)];
            RT[idx + (z * (3 * Nx / 2 + 2) + x)] =
                RT[idx + (z * (3 * Nx / 2 + 2) + x)] * RT[idx +
                                                          (z *
                                                           (3 * Nx / 2 +
                                                            2) + x)];
        }
    }

    /*z=0;
       x=0;
       printf(" %f, %f, %f, %f, %f, %f\n",  RT[(z*(3*Nx/2+2)+x)], RT[idx+(z*(3*Nx/2+2)+x)], 
       RT[idx*2+(z*(3*Nx/2+2)+x)], RT[idx*3+(z*(3*Nx/2+2)+x)],2*sin(2*M_PI*(n*dt-dt))*2*sin(2*M_PI*(n*dt-dt)), 
       2*sin(2*M_PI*(n*dt-dt-FT))*2*sin(2*M_PI*(n*dt-dt)) ); */

    /* Fourier transform to get Uxb hats and Uzb hats. */
    for (i = 0; i < 4; ++i) {

        /* Each row of RT[i] */
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
    /*      z=0;
       x=0;
       printf(" %f, %f, %f, %f, %f, %f\n",  Re(CT[0][z][x]), Re(CT[1][z][x]), Re(CT[2][z][x]),
       Re(CT[3][z][x]), 2*sin(2*M_PI*(n*dt-dt))*2*sin(2*M_PI*(n*dt-dt)), 
       2*sin(2*M_PI*(n*dt-dt-FT))*2*sin(2*M_PI*(n*dt-dt)) ); */
    for (z = 0; z < Nz / 2; z++) {
        for (x = 0; x < Nx / 2; x++) {
            Re(grad[z][x]) -=
                (re * (Re(CT[0][z][x]) + Re(CT[1][z][x])) / 2. +
                 re * (Re(CT[2][z][x]) + Re(CT[3][z][x]))) * dt;
            Im(grad[z][x]) -=
                (re * (Im(CT[0][z][x]) + Im(CT[1][z][x])) / 2. +
                 re * (Im(CT[2][z][x]) + Im(CT[3][z][x]))) * dt;
        }
    }
    for (z = Nz / 2 + 1; z < Nz; z++) {
        for (x = 0; x < Nx / 2; x++) {
            Re(grad[z][x]) -=
                (re *
                 (Re(CT[0][z + Nz / 2][x]) +
                  Re(CT[1][z + Nz / 2][x])) / 2. +
                 re * (Re(CT[2][z + Nz / 2][x]) +
                       Re(CT[3][z + Nz / 2][x]))) * dt;
            Im(grad[z][x]) -=
                (re *
                 (Im(CT[0][z + Nz / 2][x]) +
                  Im(CT[1][z + Nz / 2][x])) / 2. +
                 re * (Im(CT[2][z + Nz / 2][x]) +
                       Im(CT[3][z + Nz / 2][x]))) * dt;
        }
    }

    return (NO_ERR);
}

/* tangent forcing for Re variation when (Kx,Kz) != (0,0) */
#include <stdio.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include <assert.h>
#include "arrays.h"
#include "minChnl.h"
#include "mvOps.h"
#include "channel.h"

void tangent_manu(int n, mcomplex **** f)
{
    extern int qpts, dimR, Nx, Nz;
    extern double dt;
    extern double **Rw, **Rs, *Qy;
    extern mcomplex ****IC;
    extern mcomplex *****MIC;
    extern double **MZ;

    static mcomplex ***** MIC_tmp = NULL;
    static mcomplex **** IC_tmp = NULL;
    static mcomplex **** IC_manu = NULL;
    int i, j, z, count;

    if (MIC_tmp == NULL) {
        MIC_tmp = c5Darray(3, Nz, 2, dimR, Nx / 2);
        IC_tmp = c4Darray(Nz, 2, dimR, Nx / 2);
        IC_manu = c4Darray(Nz, 2, dimR, Nx / 2);
        assert (MIC_tmp != NULL);
        assert (IC_tmp != NULL);
        assert (IC_manu != NULL);

        memset(IC_manu[0][0][0], 0, Nz * 2 * dimR * Nx / 2 * sizeof(mcomplex));

        for (i = 0; i < dimR; ++i) {
            for (j = 0; j < T_RSDIAG; ++j) {
                MZ[i][j] = Rs[i][j]; // + re * b[k] * dt * Rps[i][j];
            }
        }
        for (i = 0; i < dimR; ++i) {
            for (j = 0; j < qpts; ++j) {
                Re(IC_manu[0][ALPHA][i][0]) +=
                    Rw[i][j] * (1 - Qy[j]*Qy[j]) * Qy[j];
            }
        }
        bsolve0(MZ, IC_manu[0][ALPHA], RSDIAG - 1, RSDIAG - 1, dimR);
    }

    count = n * 3;
    memmove(MIC_tmp[0][0][0][0], MIC[count][0][0][0],
            3 * (Nz) * 2 * dimR * (Nx / 2) * sizeof(mcomplex));
    memmove(IC_tmp[0][0][0], IC[0][0][0],
            (Nz) * 2 * dimR * (Nx / 2) * sizeof(mcomplex));
    memmove(IC[0][0][0], IC_manu[0][0][0],
            (Nz) * 2 * dimR * (Nx / 2) * sizeof(mcomplex));

    tangent(n, n+1, IC[0][0][0], 0);
    for (z = 0; z < (Nz) * 2 * dimR * (Nx / 2); ++z) {
        Re(f[0][0][0][z]) = (Re(IC[0][0][0][z]) - Re(IC_manu[0][0][0][z])) / dt;
        Im(f[0][0][0][z]) = (Im(IC[0][0][0][z]) - Im(IC_manu[0][0][0][z])) / dt;
    }

    memmove(MIC[count][0][0][0], MIC_tmp[0][0][0][0],
           3 * (Nz) * 2 * dimR * (Nx / 2) * sizeof(mcomplex));
    memmove(IC[0][0][0], IC_tmp[0][0][0],
            (Nz) * 2 * dimR * (Nx / 2) * sizeof(mcomplex));
}


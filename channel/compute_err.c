#include <stdio.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include "minChnl.h"
#include "mvOps.h"

double compute_err(mcomplex **** A, mcomplex **** B)
{
    extern int qpts, Nx, Nz;
    extern double *W;
    int y, z, x, i;
    double err = 0;

    for (y = 0; y < qpts; y++) {
        for (z = 0; z < Nz; ++z) {
            for (x = 0; x < Nx / 2; ++x) {
                for (i = 0; i < 3; i++) {
                    err += (Re(A[z][i][y][x]) - Re(B[z][i][y][x])) *
                        (Re(A[z][i][y][x]) - Re(B[z][i][y][x])) * W[y]
                        + (Im(A[z][i][y][x]) - Im(B[z][i][y][x]))
                        * (Im(A[z][i][y][x]) - Im(B[z][i][y][x])) * W[y];
                }
            }
        }
    }

    return (err);

}

/* tangent forcing for Re variation when (Kx,Kz) != (0,0) */
#include <stdio.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include "minChnl.h"
#include "mvOps.h"


void tangent_forcing0(int n, int k, int flag, mcomplex * fb, mcomplex * fa)
{
    /* External Variables */
    extern int qpts, dimR, dimQ, Nx;
    extern double dt, re;
    extern double *Kx, *Kz, **K2;
    extern mcomplex **Fa, **Fb, **TM;
    extern double **Q, **Qp, **Qpp, **R, **Rp, **Qw, **Qpw, **Rw, **Qs,
        **Qps, **Qpps, **Rs, **Rps, *Rp0;
    extern double **Mz;
    extern mcomplex ****U, ****C;
    extern mcomplex *****MC;
   /* Local variables */
    int i, j, x;
    double s;

    /* ALPHA Forcing */

    /* compute diffusion matrix M = Dv for ALPHAS  */
    if (k == 0) {               /* first step */
        for (i = 0; i < dimR; ++i) {    /* MZ = -(1/RE^2) * D0 */
            for (j = 0; j < T_RSDIAG; ++j) {
                MZ[i][j] = - re * re * Rps[i][j]; 
            }
        }
        /* compute MZ*C[0][ALPHA] and store the result in fa. */
        smMult0(MZ, C[0][ALPHA], fa, RSDIAG - 1, RSDIAG - 1, dimR);
        
    }

    /* diffusion matrix M = -(1/RE)^2 * D0 for ALPHAS */
    for (i = 0; i < dimR; ++i) {
        for (j = 0; j < T_RSDIAG; ++j) {
            MZ[i][j] = - re * re * Rps[i][j];  
        }
    }

    /* Compute alpha forcing  */
    smMult0(MZ, C[0][ALPHA], fa, RSDIAG - 1, RSDIAG - 1, dimR);    

    /* BETA Forcing */
    
    /* Compute beta forcing  */
    smMult0(MZ, C[0][BETA], fb, RSDIAG - 1, RSDIAG - 1, dimR);
    
}

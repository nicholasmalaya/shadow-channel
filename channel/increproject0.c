/***********************************************************
New piece of code to solve the incremental state equation for
 (Kx, Kz)=(0,0). The incremental state equation will be discretized exactly
as the state equation.
***********************************************************************/
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "minChnl.h"
#include "mvOps.h"

/* note for myself: need a new mpg for the incremental state solution */

void increproject0(int count, int k, int flag, func_force_t force)
{

    /* External Variables */
    extern int qpts, dimR;
    extern double dt, re;
    extern mcomplex *Ifa, *Ifb, *Itm;
    extern double **R, **Rp, **Rw, **Rs, **Rps;
    extern double **MZ;
    extern mcomplex ****IU, ****IC;
    extern mcomplex *****MIC;
    /* Local variables */
    int i, j;
    double flux_t;

    /* static double a[3] = {29./96., -3./40.,   1./6.};
       static double b[3] = {37./160., 5./24.,   1./6.};
       static double c[3] = { 8./15.,  5./12.,   3./4.};
       static double d[3] = { 0.,    -17./60.,  -5./12.};
       static  double h[3] = { 0., 8./15., 2./3.}; */

    static double a[3] = { 1. / 3., -1. / 2, 1. / 3. };
    static double b[3] = { 1. / 6, 2. / 3, 0. };
    static double c[3] = { 1. / 2, 1. / 3., 1., };
    static double d[3] = { 0., -1. / 6, -2. / 3 };

    //static double e[4] = { 2. / 3., -2. / 15., -2. / 15., 2. / 35. };
    //static double h[3] = { 1. / 2., 2. / 3., 1. };
    //mcomplex tmp[dimR], tmp2[dimR], add[dimR];
    mcomplex add[dimR];


/* FORCING FOR REYNOLDS NUMBER SENSITIVITY */
//    /* Apply Forcing */
//    if ((force != NULL)&&(k==0)) {
//        memset(Ifa, 0, dimR * sizeof(mcomplex));
//        memset(Ifb, 0, dimR * sizeof(mcomplex));
//        force(count, k, 1, Ifa, Ifb);
//
//        /* Solve for forcing cotribution to da/dt, db/dt */
//        /* MZ = M0 */
//        //memset(MZ[0], 0, dimR * 9 * sizeof(double));
//        for (i = 0; i < dimR; ++i) {
//            for (j = 0; j < T_RSDIAG; ++j) {
//                MZ[i][j] = Rs[i][j];
//            }
//        }
//        
//        bsolve_0(MZ, Ifa, RSDIAG - 1, RSDIAG - 1, dimR);
//        
//        for (i = 0; i < dimR; ++i) {
//            for (j = 0; j < T_RSDIAG; ++j) {
//                MZ[i][j] = Rs[i][j];
//            }
//        }
//
//        bsolve_0(MZ, Ifb, RSDIAG - 1, RSDIAG - 1, dimR);
//
////    f_a_norm = 0.0;
////    f_b_norm = 0.0;
////    for (i = 0; i < dimR; ++i) {
////       f_a_norm += pow(MAGNITUDE(Ifa[i]), 2);
////       f_b_norm += pow(MAGNITUDE(Ifb[i]), 2);
////    }
////    printf("f_a_norm = %e, f_b_norm = %e\n", sqrt(f_a_norm),sqrt(f_b_norm));
//
//
//        for (i = 0; i < dimR; ++i) {
//            Re(IC[0][ALPHA][i][0]) += dt * Re(Ifa[i]); // re * re
//            Im(IC[0][ALPHA][i][0]) += dt * Im(Ifa[i]); // re * re
//            Re(IC[0][BETA][i][0]) += dt * Re(Ifb[i]);  // re * re
//            Im(IC[0][BETA][i][0]) += dt * Im(Ifb[i]);  // re * re
//
//        }
//
//
//    }


    memset(Ifa, 0, dimR * sizeof(mcomplex));
    memset(Ifb, 0, dimR * sizeof(mcomplex));

    /* Create matrices for solving linear system. 
       Right hand side of system:  If this is the first step in the Runge-Kutta
       scheme, compute 
       [M0 + (1/RE)a[k]dt*D0]*IC[0][ALPHA] + dt*c[k]Fu. 
       [M0 + (1/RE)a[k]dt*D0]*C[0][BETA] + dt*c[k]Fw. 
       Otherwise, compute C[0][ALPHA] + dt*c[k]Fu, 
       C[0][BETA] + dt*c[k]Fw, 
       Lhs matrix:  M0 - (1/RE)b[k]dt*D0                          */

    if (k == 0) {               /* first step */
        //memset(MZ[0], 0, dimR * 9 * sizeof(double));
        for (i = 0; i < dimR; ++i) {    /* MZ = M0 + (1/RE)a[k]dt*D0 */
            for (j = 0; j < T_RSDIAG; ++j) {
                MZ[i][j] = Rs[i][j] - re * a[0] * dt * Rps[i][j];
            }
        }

        /* compute [M0 + (1/RE)a[k]dt*D0]*C[0][ALPHA] for x=0 and store the
           result in tm.  Then transfer the result back to C[0][ALPHA].  
           Do the same for the BETAS */
        smMult0(MZ, IC[0][ALPHA], Itm, RSDIAG - 1, RSDIAG - 1, dimR);
        for (i = 0; i < dimR; ++i) {
            Re(IC[0][ALPHA][i][0]) = Re(Itm[i]);
            Im(IC[0][ALPHA][i][0]) = Im(Itm[i]);
        }

        smMult0(MZ, IC[0][BETA], Itm, RSDIAG - 1, RSDIAG - 1, dimR);
        for (i = 0; i < dimR; ++i) {
            Re(IC[0][BETA][i][0]) = Re(Itm[i]);
            Im(IC[0][BETA][i][0]) = Im(Itm[i]);
        }
    }

    /* MZ = M0 - (1/RE)b[k]dt*D0 */
    //memset(MZ[0], 0, dimR * 9 * sizeof(double));
    for (i = 0; i < dimR; ++i) {
        for (j = 0; j < T_RSDIAG; ++j) {
            MZ[i][j] = Rs[i][j]; // + re * b[k] * dt * Rps[i][j];
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
    //memset(MZ[0], 0, dimR * 9 * sizeof(double));
    for (i = 0; i < dimR; ++i) {
        for (j = 0; j < T_RSDIAG; ++j) {
            MZ[i][j] = Rs[i][j] + re * b[k] * dt * Rps[i][j];
        }
    }

    /* Finish computing the right hand size */
    /* array Ifa */
    //memset(Ifa, 0, dimR * sizeof(mcomplex));
    for (i = 0; i < dimR; ++i) {
        for (j = 0; j < qpts; ++j) {
            /* no force */
            Re(Ifa[i]) += Rw[i][j] * (Re(IU[0][HXEL][j][0]));
        }
        //if (force != NULL) {
        //    Re(Ifa[i]) += Re(tmp[i]);
        //    Im(Ifa[i]) += Im(tmp[i]);
        //}
        Re(IC[0][ALPHA][i][0]) += dt * c[k] * Re(Ifa[i]);
    }

    /* Compute a's */
    bsolve0(MZ, IC[0][ALPHA], RSDIAG - 1, RSDIAG - 1, dimR);

    /* MZ = M0 - (1/RE)b[k]dt*D0 */
    //memset(MZ[0], 0, dimR * 9 * sizeof(double));
    for (i = 0; i < dimR; ++i) {
        for (j = 0; j < T_RSDIAG; ++j) {
            MZ[i][j] = Rs[i][j] + re * b[k] * dt * Rps[i][j];
        }
    } 

    /* array Ifb */
    //memset(Ifb, 0, dimR * sizeof(mcomplex));
    for (i = 0; i < dimR; ++i) {
        for (j = 0; j < qpts; ++j) {
            Re(Ifb[i]) += Rw[i][j] * Re(IU[0][HZEL][j][0]);
        }
        //if (force != NULL) {
        //    Re(Ifb[i]) += Re(tmp2[i]);
        //    Im(Ifb[i]) += Im(tmp2[i]);
        //}
        Re(IC[0][BETA][i][0]) += dt * c[k] * Re(Ifb[i]);
    }

    /* Compute b's */
    bsolve0(MZ, IC[0][BETA], RSDIAG - 1, RSDIAG - 1, dimR);


    /* NOW COMPUTE IU HATS */
    /* Iux_hat = R*a+c3/2*(1-y) */
    flux_t = 0;

    for (i = 0; i < dimR; ++i) {
        for (j = 0; j < qpts; ++j) {
            flux_t += Rw[i][j] * Re(IC[0][ALPHA][i][0]);
        }
    }


    /* FORCING FOR MEAN X VELOCITY (meanU) SENSITIVITY */
    if (force != NULL) {
        flux_t = 2.0 - flux_t;
    } else {
        flux_t = - flux_t;
    }


    for (j = 0; j < dimR; ++j)
       {
       Re(IC[0][ALPHA][j][0])= Re(IC[0][ALPHA][j][0])+flux_t*Re(add[j]) * 0.5;
       Im(IC[0][ALPHA][j][0])= Im(IC[0][ALPHA][j][0])+flux_t*Im(add[j]) * 0.5;
       } 

    for (i = 0; i < qpts; ++i) {
        Re(IU[0][XEL][i][0]) = 0.0;
        Im(IU[0][XEL][i][0]) = 0.0;
        for (j = 0; j < dimR; ++j) {
            Re(IU[0][XEL][i][0]) += R[i][j] * Re(IC[0][ALPHA][j][0]);
            Im(IU[0][XEL][i][0]) += R[i][j] * Im(IC[0][ALPHA][j][0]);
        }
    }

    /* dIux_hat = Rp*a-c3/2 */
    for (i = 0; i < qpts; ++i) {
        Re(IU[0][DXEL][i][0]) = 0.0;
        Im(IU[0][DXEL][i][0]) = 0.0;
        for (j = 0; j < dimR; ++j) {
            Re(IU[0][DXEL][i][0]) += Rp[i][j] * Re(IC[0][ALPHA][j][0]);
            Im(IU[0][DXEL][i][0]) += Rp[i][j] * Im(IC[0][ALPHA][j][0]);
        }
    }

    /* Iuy_hat = 0 */
    for (i = 0; i < qpts; ++i) {
        Re(IU[0][YEL][i][0]) = 0.0;
        Im(IU[0][YEL][i][0]) = 0.0;
    }

    // Iuz_hat = R*b +c4/2*(1-y) /
    for (i = 0; i < qpts; ++i) {
        Re(IU[0][ZEL][i][0]) = 0.0;
        Im(IU[0][ZEL][i][0]) = 0.0;
        for (j = 0; j < dimR; ++j) {
            Re(IU[0][ZEL][i][0]) += R[i][j] * Re(IC[0][BETA][j][0]);
            Im(IU[0][ZEL][i][0]) += R[i][j] * Im(IC[0][BETA][j][0]);
        }
        //Re(IU[0][ZEL][i][0]) += Re(Uzb[0][0]) * Uadd[i] / 2.;
        //Im(IU[0][ZEL][i][0]) += Im(Uzb[0][0]) * Uadd[i] / 2.;
    }

    // duz_hat = Rp*b-c4/2 /
    for (i = 0; i < qpts; ++i) {
        Re(IU[0][DZEL][i][0]) = 0.0;
        Im(IU[0][DZEL][i][0]) = 0.0;
        for (j = 0; j < dimR; ++j) {
            Re(IU[0][DZEL][i][0]) += Rp[i][j] * Re(IC[0][BETA][j][0]);
            Im(IU[0][DZEL][i][0]) += Rp[i][j] * Im(IC[0][BETA][j][0]);
        }
        //Re(IU[0][DZEL][i][0]) += -Re(Uzb[0][0]) / 2.;
        //Im(IU[0][DZEL][i][0]) += -Im(Uzb[0][0]) / 2.;
    }

    if (count >= 0) {
        for (i = 0; i < dimR; i++) {
            Re(MIC[count][0][ALPHA][i][0]) = Re(IC[0][ALPHA][i][0]);
            Im(MIC[count][0][ALPHA][i][0]) = Im(IC[0][ALPHA][i][0]);
            Re(MIC[count][0][BETA][i][0]) = Re(IC[0][BETA][i][0]);
            Im(MIC[count][0][BETA][i][0]) = Im(IC[0][BETA][i][0]);
        }
    }

    /* UPDATE RHS FOR NEXT TIME */
    if (k != 2) {               /* not last step */
        /* MZ = M0 + (1/RE)a[k+1]dt*D0 */
        //memset(MZ[0], 0, dimR * 9 * sizeof(double));
        for (i = 0; i < dimR; ++i) {
            for (j = 0; j < T_RSDIAG; ++j) {
                MZ[i][j] = Rs[i][j] - re * a[k + 1] * dt * Rps[i][j];
            }
        }

        /* compute [M0 + (1/RE)a[k+1]dt*D0]*IC[0][ALPHA].  Then update
           C[0][ALPHA].  Do the same for the b's. */
        smMult0(MZ, IC[0][ALPHA], Itm, RSDIAG - 1, RSDIAG - 1, dimR);
        for (i = 0; i < dimR; ++i) {
            Re(IC[0][ALPHA][i][0]) =
                Re(Itm[i]) + dt * d[k + 1] * Re(Ifa[i]);
            Im(IC[0][ALPHA][i][0]) =
                Im(Itm[i]) + dt * d[k + 1] * Im(Ifa[i]);
        }

        /* now b's */
        smMult0(MZ, IC[0][BETA], Itm, RSDIAG - 1, RSDIAG - 1, dimR);
        for (i = 0; i < dimR; ++i) {
            Re(IC[0][BETA][i][0]) =
                Re(Itm[i]) + dt * d[k + 1] * Re(Ifb[i]);
            Im(IC[0][BETA][i][0]) =
                Im(Itm[i]) + dt * d[k + 1] * Im(Ifb[i]);
        }
    }


}

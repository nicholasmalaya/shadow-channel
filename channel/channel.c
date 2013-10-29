/**********************************************************************
main code used to solve State, incremental state, adjoint and incremental adjoint equations.

update: 03/09, working on checking the accuracy of new rk scheme for state equations.
**********************************************************************/

#include <assert.h>
#include <stdlib.h>
#include "main.h"
#include "channel.h"
#include "hdf5.h"
#include "mcomplex.h"

#include "minChnl.h"
#include "mvOps.h"

// convertin main() to init() function
//
// will be called from python to compute primal trajectory
//
// call pattern:
//               * mesh def
//               * params
//               * delta t
//               * run-up time
//               * time chunks
//               * size of time chunks
//
// store all state in memory
//
// nick 
// 9/12/13
//

/***************************************************************
*                                                              *
*                         DESTROY FUNCTION                     *
*                                                              *
****************************************************************/

#define DESTROY_STATUS_FFTW 0x1
#define DESTROY_STATUS_GETMEM 0x2
#define DESTROY_STATUS_WAVENUMS 0x4
#define DESTROY_STATUS_LEGENDRE 0x8
#define DESTROY_STATUS_EVERYTHING 0xffff

void destroy(int status)
{
    /* External Variables.  All external variables are defined in main.h */
    extern double *Kx, *Kz, **K2, *cfl2;
    extern double **Q, **Qp, **Qpp, **R, **Rp, **Qw, **Qpw, **Rw, **Qs,
        **Qps, **Qpps, **Rs, **Rps, *Rp0, **Rpw, **Qppw, *Rpp0;

    extern double *Uadd, *Vadd, *Vpadd;
    extern double *Qy;
    extern double *W;

    extern mcomplex ****U, ****C;    /* state variables */
    extern mcomplex **Fa, **Fb, **TM;
    extern mcomplex *fa, *fb, *tm;
    extern double **MZ;
    extern double ***M;

    extern mcomplex ****IU, ****IC;    /* incremental state variables */
    extern mcomplex **IFa, **IFb, **ITM;
    extern mcomplex *Ifa, *Ifb, *Itm;

    extern mcomplex ****AU, ****AC;    /* adjoint variables and will use
                       the same other variables
                       used in state equations */

    extern mcomplex ****IAU, ****IAC;    /* incremental adjoint variables */

    extern mcomplex **Uxbt, **Uzb;    /* variables used to store dux duz
                       evaluated at y=-1 used for
                       computing boundary conditions for
                       incremental state equations */
    extern mcomplex **Uxb, **Uzb;    /* variables used to store dux duz
                       evaluated at y=-1 from previous 
                       state used for boundary conditions
                       for incremental state equations */
    extern mcomplex **IUxb, **IUzb;
    extern mcomplex **IAUxb, **IAUzb;
    extern mcomplex **AUxb, **AUzb;    /* variables used to store dux duz
                       evaluated at y=-1 used for
                       computing boundary conditions
                       for incremental state equations */
    extern fftw_complex ***CT, ***ICT;    /* variables used in fft */
    extern fftw_plan pf1, pf2;
    extern fftw_plan Ipf1, Ipf2;
    extern rfftwnd_plan pr1, pr2;

    extern mcomplex *****MC, *****MIC;    /* variables used to store state and
                           incremental state solutions
                           between two check points. */

    extern mcomplex ****MU, ****MIU;    /* variables used to store
                           manufacture solutions */
    extern mcomplex ****LU, ****LIU;

    mcomplex tmp;

    if (status & DESTROY_STATUS_FFTW) {
        fftw_destroy_plan(pf1);
        fftw_destroy_plan(pf2);
        rfftwnd_destroy_plan(pr1);
        rfftwnd_destroy_plan(pr2);
        fftw_destroy_plan(Ipf1);
        fftw_destroy_plan(Ipf2);
    }

    if (status & DESTROY_STATUS_GETMEM) {
        freec4Darray(U);
        freec4Darray(C);
        freec3Darray(CT);
        freecMatrix(Fa);
        freecMatrix(Fb);
        freed3Darray(M);
        freecMatrix(TM);
        freedMatrix(MZ);
        freecVector(fa);
        freecVector(fb);
        freecVector(tm);
        freedVector(cfl2);
        freec4Darray(IU);
        freec4Darray(IC);
        freecMatrix(IFa);
        freecMatrix(IFb);
        freecMatrix(ITM);
        freecVector(Ifa);
        freecVector(Ifb);
        freecVector(Itm);
        freec3Darray(ICT);
        freecMatrix(Uxbt);
        freecMatrix(Uzbt);
        freecMatrix(Uxb);
        freecMatrix(Uzb);
        freec4Darray(AU);
        freec4Darray(IAU);
        freec4Darray(AC);
        freec4Darray(IAC);
        freec5Darray(MC);
        freec5Darray(MIC);
        freec4Darray(MU);
        freec4Darray(MIU);
        freec4Darray(LU);
        freec4Darray(LIU);
        freecMatrix(AUxb);
        freecMatrix(AUzb);
        freecMatrix(IUzb);
        freecMatrix(IUxb);
        freecMatrix(IAUxb);
        freecMatrix(IAUzb);
        freecMatrix(grad);
        freecMatrix(GUxb);
        freecMatrix(GUzb);
        freecMatrix(GIUxb);
        freecMatrix(GIUzb);
        freecMatrix(HUxb);
        freecMatrix(HUzb);
        freecMatrix(HAUxb);
        freecMatrix(HAUzb);
        freecMatrix(hess);
    }

    if (status & DESTROY_STATUS_WAVENUMS) {
        freedVector(Kx);
        freedVector(Kz);
        freedMatrix(K2);
    }

    if (status & DESTROY_STATUS_LEGENDRE) {
        freedMatrix(Q);
        freedMatrix(Qp);
        freedMatrix(Qpp);
        freedMatrix(R);
        freedMatrix(Rp);
        freedMatrix(Qw);
        freedMatrix(Qpw);
        freedMatrix(Rw);
        freedMatrix(Qs);
        freedMatrix(Qps);
        freedMatrix(Qpps);
        freedMatrix(Rs);
        freedMatrix(Rps);
        freedVector(Rp0);
        freedVector(W);
        freedVector(Vadd);
        freedVector(Vpadd);
        freedVector(Uadd);
        freedVector(Qy);
        freedVector(Rpp0);
        freedMatrix(Qppw);
        freedMatrix(Rpw);
    }
}

/***************************************************************
*                                                              *
*                       INIT FUNCTION                          *
*                                                              *
****************************************************************/

int
init(int _Nx, int _Ny, int _Nz, double _Lx, double _Lz, double _Re,
     double _flux, double _dt, int _ru_steps,
     int _n_chunk, int _nsteps_chunk, int _restart_flag)
{
    /******************** Definition of all variables ********************/
    /* External Variables.  All external variables are defined in main.h */
    extern int qpts, dimR, dimQ, Nx, Nz;
    extern double dt, re, flux;

    extern double *Kx, *Kz, **K2, *cfl2;
    extern double **Q, **Qp, **Qpp, **R, **Rp, **Qw, **Qpw, **Rw, **Qs,
        **Qps, **Qpps, **Rs, **Rps, *Rp0, **Rpw, **Qppw, *Rpp0;

    extern double *Uadd, *Vadd, *Vpadd;
    extern double *Qy;
    extern double *W;

    extern mcomplex ****U, ****C;    /* state variables */
    extern mcomplex **Fa, **Fb, **TM;
    extern mcomplex *fa, *fb, *tm;
    extern double **MZ;
    extern double ***M;

    extern mcomplex ****IU, ****IC;    /* incremental state variables */
    extern mcomplex **IFa, **IFb, **ITM;
    extern mcomplex *Ifa, *Ifb, *Itm;

    extern mcomplex ****AU, ****AC;    /* adjoint variables and will use
                       the same other variables
                       used in state equations */

    extern mcomplex ****IAU, ****IAC;    /* incremental adjoint variables */

    extern mcomplex **Uxbt, **Uzb;    /* variables used to store dux duz
                       evaluated at y=-1 used for
                       computing boundary conditions for
                       incremental state equations */
    extern mcomplex **Uxb, **Uzb;    /* variables used to store dux duz
                       evaluated at y=-1 from previous 
                       state used for boundary conditions
                       for incremental state equations */
    extern mcomplex **IUxb, **IUzb;
    extern mcomplex **IAUxb, **IAUzb;
    extern mcomplex **AUxb, **AUzb;    /* variables used to store dux duz
                       evaluated at y=-1 used for
                       computing boundary conditions
                       for incremental state equations */
    extern fftw_complex ***CT, ***ICT;    /* variables used in fft */
    extern fftw_plan pf1, pf2;
    extern fftw_plan Ipf1, Ipf2;
    extern rfftwnd_plan pr1, pr2;

    extern mcomplex *****MC, *****MIC;    /* variables used to store state and
                           incremental state solutions
                           between two check points. */

    extern mcomplex ****MU, ****MIU;    /* variables used to store
                           manufacture solutions */
    extern mcomplex ****LU, ****LIU;

    /* Local Variables */
    int n, z, dctr;
    int Ny, sizeRealTransform;
    double Lx, Lz;
    fftw_complex *fout = NULL;
    int restart_flag;
    int count;
    // int checknum, checkstep;
    int ru_steps;
    int i;

    /************************ end of variable definitions ****************/

    Nx = _Nx;
    Ny = _Ny;
    Nz = _Nz;
    Lx = _Lx;
    Lz = _Lz;

    dt = _dt;

    nsteps = _n_chunk * _nsteps_chunk;
    ru_steps = _ru_steps;
    flux = _flux;
    re = _Re;

    restart_flag = _restart_flag;

    printf("Nx,Ny,Nz,Lx  |  Lz,dt,nsteps,ru_steps  |  flux,Re,restart_flag\n"
           "%d %d %d %f  |  %f %f %d %d  |  %f %f %d\n",
           Nx, Ny, Nz, Lx, Lz, dt, nsteps, ru_steps, flux, re, restart_flag);

    re = 1. / re;         /* time step routines assume I pass 1/Re */
    qpts = 3 * Ny / 2;    /* number of quadrature points 
                   (see page 9 of Moser's notes) */
    dimR = Ny - 2;        /* dimR and dimQ denote the number of terms */
    dimQ = Ny - 4;        /* in the truncated expansions for the */

    /* functions v_hat, g_hat, U, W */
    /* (see page 5 of Moser's notes). */
    sizeRealTransform = 3 * Nx / 2;    /* for the FFTs */

    /**************** check input parameters, Nx/4, Nz/2 ***************/
    if (Nx % 4 != 0) {
        printf("Required arguments Ny/4==0\n");
        return (EXIT_FAILURE);
    }
    if (Nz % 2 != 0) {
        printf("Required arguments Nz/2==0\n");
        return (EXIT_FAILURE);
    }
    if (Ny - 4 < 0) {
        printf("Required arguments Nzy>4\n");
        return (EXIT_FAILURE);
    }
    /* Create matrices using Legendre polynomials */
    if (LegendreSetup() != NO_ERR) {
        return (EXIT_FAILURE);
    }

    /*********************end of parameter checking ****************/

    /****************Initialize and allocate all variables ***************/
    /* Compute wave numbers */
    if (waveNums(Nx / 2, Nz, Lx, Lz) != NO_ERR) {
        destroy(DESTROY_STATUS_LEGENDRE);
        return (EXIT_FAILURE);
    }

    /* get memory for 4D arrays and other matrices */
    if (getMem() != NO_ERR) {
        destroy(DESTROY_STATUS_LEGENDRE | DESTROY_STATUS_WAVENUMS);
        return (EXIT_FAILURE);
    }

    /* Create plans for FFTs */
    pf1 = fftw_create_plan_specific(3 * Nz / 2, FFTW_BACKWARD,
                    FFTW_MEASURE | FFTW_IN_PLACE, CT[0][0],
                    3 * Nx / 4 + 1, fout, -1);
    pf2 = fftw_create_plan_specific(3 * Nz / 2, FFTW_FORWARD,
                      FFTW_MEASURE | FFTW_IN_PLACE, CT[0][0],
                      3 * Nx / 4 + 1, fout, -1);
    pr1 = rfftwnd_create_plan(1, &sizeRealTransform, FFTW_COMPLEX_TO_REAL,
                FFTW_MEASURE | FFTW_IN_PLACE);
    pr2 = rfftwnd_create_plan(1, &sizeRealTransform, FFTW_REAL_TO_COMPLEX,
                FFTW_MEASURE | FFTW_IN_PLACE);

    /* Create plans for FFTs */
    Ipf1 = fftw_create_plan_specific(3 * Nz / 2, FFTW_BACKWARD,
                     FFTW_MEASURE | FFTW_IN_PLACE,
                     ICT[0][0], 3 * Nx / 4 + 1, fout, -1);
    Ipf2 = fftw_create_plan_specific(3 * Nz / 2, FFTW_FORWARD,
                      FFTW_MEASURE | FFTW_IN_PLACE, ICT[0][0],
                      3 * Nx / 4 + 1, fout, -1);

    /* set variables for checking CFL condition */
    if (cflVars(Lx, Lz) != 0) {
        printf("Error creating CF variables\n");
        destroy(DESTROY_STATUS_LEGENDRE | DESTROY_STATUS_WAVENUMS
            | DESTROY_STATUS_FFTW);
        return (EXIT_FAILURE);
    }

    /* initalize part */
    memset(C[0][0][0], 0,
           (Nz) * 2 * (Ny - 2) * (Nx / 2) * sizeof(mcomplex));
    memset(IC[0][0][0], 0,
           (Nz) * 2 * (Ny - 2) * (Nx / 2) * sizeof(mcomplex));
    memset(AC[0][0][0], 0,
           (Nz) * 2 * (Ny - 2) * (Nx / 2) * sizeof(mcomplex));
    memset(IAC[0][0][0], 0,
           (Nz) * 2 * (Ny - 2) * (Nx / 2) * sizeof(mcomplex));

    memset(MC[0][0][0][0], 0,
           (nsteps * 3 + 1) * (Nz) * 2 * (Ny - 2) * (Nx / 2)
           * sizeof(mcomplex));
    memset(MIC[0][0][0][0], 0,
           (nsteps * 3 + 1) * (Nz) * 2 * (Ny - 2) * (Nx / 2)
           * sizeof(mcomplex));

    memset(U[0][0][0], 0, (Nz) * 5 * qpts * (Nx / 2) * sizeof(mcomplex));
    memset(AU[0][0][0], 0, (Nz) * 5 * qpts * (Nx / 2) * sizeof(mcomplex));
    memset(IU[0][0][0], 0, (Nz) * 5 * qpts * (Nx / 2) * sizeof(mcomplex));
    memset(IAU[0][0][0], 0, (Nz) * 5 * qpts * (Nx / 2) * sizeof(mcomplex));

    memset(LU[0][0][0], 0, (Nz) * 5 * qpts * (Nx / 2) * sizeof(mcomplex));
    memset(LIU[0][0][0], 0, (Nz) * 5 * qpts * (Nx / 2) * sizeof(mcomplex));

    memset(Uxb[0], 0, Nz * (Nx / 2) * sizeof(mcomplex));
    memset(Uzb[0], 0, Nz * (Nx / 2) * sizeof(mcomplex));

    /******************************end of initialization part ***************/

    /******************restart check ******************************/
    if (restart_flag != 0) {
        restart2(restart_flag);
    } else            // provide laminar solution with perturbations
    {
        memset(U[0][0][0], 0, Nz * 5 * qpts * Nx / 2 * sizeof(mcomplex));
        for (i = 0; i < qpts; i++) {
            Re(U[0][XEL][i][0]) = (1.0 - Qy[i] * Qy[i]) * flux * 3./4.;
            Im(U[0][XEL][i][0]) = 0.0;
        }
        // perturbation
        for (i = 0; i < qpts; i++) {
            Re(U[1][XEL][i][1]) = (rand() / (double)RAND_MAX - 0.5);
            Im(U[1][XEL][i][1]) = (rand() / (double)RAND_MAX - 0.5);
        }
        initAlphaBeta();
    }

    /**********************************end of restart **********************/

    /* store current Fourier coefficient into MC and MIC
       used for solving the adjoint system. */
    // destination // source
    memcpy(MC[0][0][0][0], C[0][0][0],
           (Nz) * 2 * (Ny - 2) * (Nx / 2) * sizeof(mcomplex));
    memcpy(MIC[0][0][0][0], IC[0][0][0],
           (Nz) * 2 * (Ny - 2) * (Nx / 2) * sizeof(mcomplex));

    /***************solving state and incremental state equations ****/

    /* time step for forward problem */
    for (n = restart_flag; n < restart_flag + nsteps + ru_steps; ++n) {
        if (n % 100 == 0) {
            printf("Step %d/%d\n", n, restart_flag + nsteps + ru_steps);
        }
        for (dctr = 0; dctr < 3; ++dctr) {    /* RK steps */

            /* copy the result to Uxbt, Uzbt. Uxb and Uzb will be used
               later for boundary condition of current time stage */
            memcpy(Uxbt[0], Uxb[0],
                   (Nz) * (Nx / 2) * sizeof(fftw_complex));
            memcpy(Uzbt[0], Uzb[0],
                   (Nz) * (Nx / 2) * sizeof(fftw_complex));
            memset(Uxb[0], 0, Nz * (Nx / 2) * sizeof(fftw_complex));
            memset(Uzb[0], 0, Nz * (Nx / 2) * sizeof(fftw_complex));

            /* do FFTs to get H_hats.  After this we have for each (Kx,y,Kz)
               Hx_hat   -->  U[z][HXEL][y][x]
               Hy_hat   -->  U[z][HYEL][y][x]
               Hz_hat   -->  U[z][HZEL][y][x]
             */
            if (pass1(dctr, n) != NO_ERR) {
                printf("Pass1 failure\n");
                n = restart_flag + nsteps + ru_steps;
                break;
            }

            project0(dctr, n, NULL);    /* (kx, kz)=(0, 0),
                               solve for a, b */
            project(n, dctr, 0, 1, NULL);    /* (kx, kz)!=(0, 0),
                               solve for alpha, beta */
            for (z = 1; z < Nz; ++z) {
                if (z == Nz / 2) {
                    memset(U[z][0][0], 0,
                           5 * qpts * (Nx / 2) * sizeof(mcomplex));
                    continue;
                }

                project(n, dctr, z, 0, NULL);
            }

            /*now update the boundary condition using current time step
               solution of state equation */
            if (increBoundary() != NO_ERR) {
                printf("increBoundary failure\n");
                n = restart_flag + nsteps + ru_steps;
                break;
            }

            increproject0(dctr, n, 1, NULL);
            increproject(dctr, 0, 1, n, NULL);
            for (z = 1; z < Nz; ++z) {
                if (z == Nz / 2) {
                    // SET U[z][XEL,YEL,ZEL,DXEL,DZEL] TO ZEROS 
                    memset(IU[z][0][0], 0,
                           5 * qpts * (Nx / 2) * sizeof(mcomplex));
                    continue;
                }

                increproject(dctr, z, 0, n, NULL);
            }

            count = (n - restart_flag - ru_steps) * 3 + dctr + 1;
            if (count >= 0) {
                memcpy(MC[count][0][0][0], C[0][0][0],
                       (Nz) * 2 * (Ny - 2) * (Nx / 2) * sizeof(mcomplex));
                memcpy(MIC[count][0][0][0], IC[0][0][0],
                       (Nz) * 2 * (Ny - 2) * (Nx / 2) * sizeof(mcomplex));
            }

        }        /* end for dctr... */

        /* now writing the results to HDF file at selected time steps */
        if (((n + 1) % 10000 == 0) ||
                (n + 1 == restart_flag + nsteps + ru_steps)) {
            /* when we need to store the current time step results */
            write_data2(n + 1);
            // checknum = checknum + 1;
            // count = 0;
            // memcpy(MC[count][0][0][0], C[0][0][0],
            //        (Nz) * 2 * (Ny - 2) * (Nx / 2) * sizeof(mcomplex));
            // memcpy(MIC[count][0][0][0], IC[0][0][0],
            //        (Nz) * 2 * (Ny - 2) * (Nx / 2) * sizeof(mcomplex));
        }
    }            /* end for n... */

    return (EXIT_SUCCESS);
}                /* end main */

/***************************************************************
*                                                              *
*                         GETSOLN FUNCTION                     *
*                                                              *
****************************************************************/

// This is supposed to be used with ARGOUTVIEW_ARRAY4
// Look for the c_v function in kuramoto.i
void
getsoln(int i_step, mcomplex ** MC_ptr,
    int *Nz_ptr, int *Nvar_ptr, int *Ny_ptr, int *Nx_ptr)
{
    assert(i_step >= 0 && i_step <= nsteps);
    (*MC_ptr) = MC[i_step * 3][0][0][0];
    (*Nz_ptr) = Nz;
    (*Nvar_ptr) = 2;
    (*Ny_ptr) = dimR;
    (*Nx_ptr) = Nx / 2;
}

/***************************************************************
*                                                              *
*                       STATISTICS FUNCTION                    *
*                                                              *
****************************************************************/

void statistics(mcomplex * C_ptr,
             int Nz_dup, int Nvar_dup, int Ny_dup, int Nx_dup,
             double ** us_ptr, int * nstats_ptr, int * qpts_ptr)
{
    extern double * Qy;
    int z;
    double ** us;

    /* copy C_ptr to C */
    assert(Nx_dup == Nx / 2);
    assert(Ny_dup == dimR);
    assert(Nz_dup == Nz);
    assert(Nvar_dup == 2);
    memcpy(C[0][0][0], C_ptr, Nz * 2 * dimR * (Nx / 2) * sizeof(mcomplex));

    /* Compute U from C */
    initAlphaBeta2();

    /* compute statistics */
    us = dMatrix(21, qpts);     // memory leak here, not much though
    memcpy(us[0], Qy, qpts * sizeof(double));
    comp_stat(us + 1);

    /* set return variables */
    *us_ptr = us[0];
    *nstats_ptr = 21;
    *qpts_ptr = qpts;
}

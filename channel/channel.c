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

// for restart2
#include <stdio.h>
#include <string.h>
#include "arrays.h"

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
     double _flux, double _dt, int _nsteps)
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
    int Ny, sizeRealTransform;
    double Lx, Lz;
    fftw_complex *fout = NULL;

    /************************ end of variable definitions ****************/

    Nx = _Nx;
    Ny = _Ny;
    Nz = _Nz;
    Lx = _Lx;
    Lz = _Lz;

    dt = _dt;

    nsteps = _nsteps;
    flux = _flux;
    re = _Re;

    printf("Nx,Ny,Nz,Lx  |  Lz,dt,nsteps  |  flux,Re\n"
           "%d %d %d %f  |  %f %f %d |  %f %f\n",
           Nx, Ny, Nz, Lx, Lz, dt, nsteps, flux, re);

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
           (Nz) * 2 * dimR * (Nx / 2) * sizeof(mcomplex));
    memset(IC[0][0][0], 0,
           (Nz) * 2 * dimR * (Nx / 2) * sizeof(mcomplex));
    memset(AC[0][0][0], 0,
           (Nz) * 2 * dimR * (Nx / 2) * sizeof(mcomplex));
    memset(IAC[0][0][0], 0,
           (Nz) * 2 * dimR * (Nx / 2) * sizeof(mcomplex));

    memset(MC[0][0][0][0], 0,
           (nsteps * 3 + 1) * (Nz) * 2 * dimR * (Nx / 2)
           * sizeof(mcomplex));
    memset(MIC[0][0][0][0], 0,
           (nsteps * 3 + 1) * (Nz) * 2 * dimR * (Nx / 2)
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
    return (EXIT_SUCCESS);
}

/***************************************************************
*                                                              *
*                        SOLVE THE PRIMAL                      *
*                                                              *
****************************************************************/

void primal(int ru_steps, mcomplex *C_given)
{
   int i, n, z, dctr, count;

    /******************restart check ******************************/
    if (C_given != 0) {
        memmove(C[0][0][0], C_given,
                (Nz) * 2 * dimR * (Nx / 2) * sizeof(mcomplex));
        initAlphaBeta2();
    }
    else            // provide laminar solution with perturbations
    {
        memset(U[0][0][0], 0, Nz * 5 * qpts * Nx / 2 * sizeof(mcomplex));
        for (i = 0; i < qpts; i++) {
            Re(U[0][XEL][i][0]) = (1.0 - Qy[i] * Qy[i]) * flux * 3./4;
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

    /* store current Fourier coefficient into MC
       used for solving the adjoint system. */
    // destination // source
    memset(MC[0][0][0][0], 0, (nsteps * 3 + 1) *
           (Nz) * 2 * dimR * (Nx / 2) * sizeof(mcomplex));
    memmove(MC[0][0][0][0], C[0][0][0],
            (Nz) * 2 * dimR * (Nx / 2) * sizeof(mcomplex));

    /***************solving state and incremental state equations ****/

    /* time step for forward problem */
    for (n = 0; n < nsteps + ru_steps; ++n) {
        if (n % 100 == 0) {
            printf("Step %d/%d\n", n, nsteps + ru_steps);
        }
        for (dctr = 0; dctr < 3; ++dctr) {    /* RK steps */

            /* do FFTs to get H_hats.  After this we have for each (Kx,y,Kz)
               Hx_hat   -->  U[z][HXEL][y][x]
               Hy_hat   -->  U[z][HYEL][y][x]
               Hz_hat   -->  U[z][HZEL][y][x]
             */
            if (pass1(dctr, n) != NO_ERR) {
                printf("Pass1 failure\n");
                n = nsteps + ru_steps;
                break;
            }

            count = (n - ru_steps) * 3 + dctr + 1;

            project0(count, dctr, NULL);    /* (kx, kz)=(0, 0),
                               solve for a, b */
            project(count, dctr, 0, 1, NULL);    /* (kx, kz)!=(0, 0),
                               solve for alpha, beta */
            for (z = 1; z < Nz; ++z) {
                if (z == Nz / 2) {
                    memset(U[z][0][0], 0,
                           5 * qpts * (Nx / 2) * sizeof(mcomplex));
                    memset(C[z][0][0], 0,
                           2 * dimR * (Nx / 2) * sizeof(mcomplex));
                    continue;
                }

                project(count, dctr, z, 0, NULL);
            }
        }        /* end for dctr... */
    }            /* end for n... */

    if (C_given != 0) {  /* copy the final solution back to given buffer */
        memmove(C_given, C[0][0][0],
                (Nz) * 2 * dimR * (Nx / 2) * sizeof(mcomplex));
    }
}

/***************************************************************
*                                                              *
*                        SOLVE THE TANGENT                     *
*                                                              *
****************************************************************/

void tangent(int start_step, int end_step, mcomplex *IC_given, int inhomo)
{
    /* Local variables */
    int n, dctr, count, z;

    /* retrieve solution */
    assert (end_step > start_step);
    assert (start_step >= 0);
    assert (end_step <= nsteps);
	assert (IC_given != 0);

    memmove(C[0][0][0], MC[start_step * 3][0][0][0],
            (Nz) * 2 * dimR * (Nx / 2) * sizeof(mcomplex));
    memmove(IC[0][0][0], IC_given,
            (Nz) * 2 * dimR * (Nx / 2) * sizeof(mcomplex));

    // memset(U[0][0][0], 0, (Nz) * 5 * qpts * (Nx / 2) * sizeof(mcomplex));
    // memset(AU[0][0][0], 0, (Nz) * 5 * qpts * (Nx / 2) * sizeof(mcomplex));
    // memset(IU[0][0][0], 0, (Nz) * 5 * qpts * (Nx / 2) * sizeof(mcomplex));
    // memset(IAU[0][0][0], 0, (Nz) * 5 * qpts * (Nx / 2) * sizeof(mcomplex));

    // memset(LU[0][0][0], 0, (Nz) * 5 * qpts * (Nx / 2) * sizeof(mcomplex));
    // memset(LIU[0][0][0], 0, (Nz) * 5 * qpts * (Nx / 2) * sizeof(mcomplex));

    // memset(Uxb[0], 0, Nz * (Nx / 2) * sizeof(mcomplex));
    // memset(Uzb[0], 0, Nz * (Nx / 2) * sizeof(mcomplex));

    initAlphaBeta2();
    incre_initAlphaBeta2();

    /* forcing function */
    func_force_t forcing0 = NULL;
    func_force_tt forcing = NULL;
    if (inhomo != 0) {
        forcing0 = tangent_forcing0;
        //forcing = tangent_forcing;
    } 

    /***************solving state and incremental state equations ****/
    count = start_step * 3;
    memmove(MIC[count][0][0][0], IC[0][0][0],
            (Nz) * 2 * dimR * (Nx / 2) * sizeof(mcomplex));

    /* time step for forward problem */
    for (n = start_step; n < end_step; ++n) {
        for (dctr = 0; dctr < 3; ++dctr) {    /* RK steps */

            /* do FFTs to get H_hats.  After this we have for each (Kx,y,Kz)
               Hx_hat   -->  U[z][HXEL][y][x]
               Hy_hat   -->  U[z][HYEL][y][x]
               Hz_hat   -->  U[z][HZEL][y][x]
             */
            if (pass1(dctr, n) != NO_ERR) {
                printf("Pass1 failure\n");
                n = nsteps;
                break;
            }

            count = n * 3 + dctr + 1;

            increproject0(count, dctr, 1, forcing0);
            increproject(count, dctr, 0, 1, forcing);
            for (z = 1; z < Nz; ++z) {
                if (z == Nz / 2) {
                    // SET U[z][XEL,YEL,ZEL,DXEL,DZEL] TO ZEROS 
                    memset(IU[z][0][0], 0,
                           5 * qpts * (Nx / 2) * sizeof(mcomplex));
                    memset(IC[z][0][0], 0,
                           2 * dimR * (Nx / 2) * sizeof(mcomplex));
                    continue;
                }

                increproject(count, dctr, z, 0, forcing);
            }

            project0(-1, dctr, NULL);    /* (kx, kz)=(0, 0), solve for a, b */
            project(-1, dctr, 0, 1, NULL);  /* (kx, kz)!=(0, 0),
                               solve for alpha, beta */
            for (z = 1; z < Nz; ++z) {
                if (z == Nz / 2) {
                    memset(U[z][0][0], 0,
                           5 * qpts * (Nx / 2) * sizeof(mcomplex));
                    memset(C[z][0][0], 0,
                           2 * dimR * (Nx / 2) * sizeof(mcomplex));
                    continue;
                }
                project(-1, dctr, z, 0, NULL);
            }
        }        /* end for dctr... */

        //if (inhomo != 0) {
        //    tangent_manu(n, IAC);
        //    for (z = 0; z < (Nz) * 2 * dimR * (Nx / 2); ++z) {
        //        Re(IC[0][0][0][z]) += dt * Re(IAC[0][0][0][z]);
        //        Im(IC[0][0][0][z]) += dt * Im(IAC[0][0][0][z]);
        //    }
        //}
    }            /* end for n... */
 
    memmove(IC_given, IC[0][0][0],
            (Nz) * 2 * dimR * (Nx / 2) * sizeof(mcomplex));
}


/***************************************************************
*                                                              *
*                        SOLVE THE ADJOINT                     *
*                                                              *
****************************************************************/

void adjoint(int start_step, int end_step, mcomplex *AC_given, int inhomo,
             double strength)
{
    int n, dctr, count, z;
    mcomplex *ICtmp;
    int len1 = (Nz * 2 * dimR * (Nx / 2));
    ICtmp = cVector(len1); 
    /* retrieve solution */
    assert (end_step < start_step);
    assert (end_step >= 0);
    assert (start_step <= nsteps);
	assert (AC_given != 0);

    memmove(C[0][0][0], MC[start_step * 3][0][0][0],
            (Nz) * 2 * dimR * (Nx / 2) * sizeof(mcomplex));
    memmove(AC[0][0][0], AC_given,
            (Nz) * 2 * dimR * (Nx / 2) * sizeof(mcomplex));


    
    /***************solving adjoint equations ****/

    adj_initAlphaBeta2();
    strength *= dt;

    /* time step for adjoint problem */
    for (n = start_step; n > end_step; --n) {

        count = n * 3;
        memmove(ICtmp, MIC[count][0][0][0],
                (Nz) * 2 * dimR * (Nx / 2) * sizeof(mcomplex));
        ddt_project(n, ICtmp);
        for (z = 0; z < (Nz) * 2 * dimR * (Nx / 2); ++ z) {
            Re(AC[0][0][0][z]) -= 0.5 * strength * Re(ICtmp[z]);
            Im(AC[0][0][0][z]) -= 0.5 * strength * Im(ICtmp[z]);
        }

        for (dctr = 0; dctr < 3; ++dctr) {    /* RK steps */

            count = n * 3 - dctr - 1;

            /*read data from memery */
            if (dctr < 2) {
                memmove(C[0][0][0], MC[count - 1][0][0][0],
                        (Nz) * 2 * dimR * (Nx / 2) * sizeof(mcomplex));

                /* reconstruct the state and incremental state solution u,
                   iu from alpha and beta */
                initAlphaBeta2();

                if (pass2(dctr, n) != NO_ERR) {
                    printf("Pass2 failure\n");
                    n = -1;
                    break;
                }
                memmove(LU[0][0][0], U[0][0][0],
                        (Nz) * 5 * qpts * (Nx / 2) * sizeof(mcomplex));
            }
            /*read data from memery */
            memmove(C[0][0][0], MC[count][0][0][0],
                    (Nz) * 2 * dimR * (Nx / 2) * sizeof(mcomplex));

            /*reconstruct the state and incremental state solution u,
              iu from alpha and beta */
            initAlphaBeta2();

            if (pass2(dctr, n) != NO_ERR) {
                printf("Pass2 failure\n");
                n = -1;
                break;
            }

            adjproject0(count, dctr, NULL);
            adjproject(count, dctr, 0, 1, NULL);
            for (z = 1; z < Nz; ++z) {
                if (z == Nz / 2) {
                    memset(AU[z][0][0], 0,
                           5 * qpts * (Nx / 2) * sizeof(mcomplex));
                    memset(AC[z][0][0], 0,
                           2 * dimR * (Nx / 2) * sizeof(mcomplex));
                    continue;
                }
                adjproject(count, dctr, z, 0, NULL);
            }
        }

        count = (n - 1) * 3;
        memmove(ICtmp, MIC[count][0][0][0],
                (Nz) * 2 * dimR * (Nx / 2) * sizeof(mcomplex));
        ddt_project(n-1, ICtmp);
        for (z = 0; z < (Nz) * 2 * dimR * (Nx / 2); ++ z) {
            Re(AC[0][0][0][z]) -= 0.5 * strength * Re(ICtmp[z]);
            Im(AC[0][0][0][z]) -= 0.5 * strength * Im(ICtmp[z]);
        }

	}

    memmove(AC_given, AC[0][0][0],
            (Nz) * 2 * dimR * (Nx / 2) * sizeof(mcomplex));
    freecVector( ICtmp );
}

/***************************************************************
*                                                              *
*                  PRIMAL DERIVATIVE PROJECTION                *
*                                                              *
****************************************************************/

double ddt_project(int i_step, mcomplex * IC_given)
{

    mcomplex *dCdt;
    double *dudt; 
    double *v; 
    double num, den;
    int x, y, z, i;
    int XFLOW, YFLOW, ZFLOW;
    int len1 = (Nz * 2 * dimR * (Nx / 2));
    int len2 = (3 * 3 * (Nx / 2) * qpts * 3 * (Nz /2));
    int i_plus, i_minus;

    // Allocate memory
    dCdt = cVector(len1);
    dudt = dVector(len2);
    v = dVector(len2);
    memset(dCdt, 0,
           (Nz) * 2 * dimR * (Nx / 2) * sizeof(mcomplex));
    memset(dudt, 0,
           3 * (3 * Nx / 2) * qpts * (3 * Nz / 2) * sizeof(double));
    memset(v, 0,
           3 * (3 * Nx / 2) * qpts * (3 * Nz / 2) * sizeof(double));

    // compute time derivative with finite difference
    if (i_step == 0) {
        i_plus = 3*(i_step + 1);
        i_minus = 3*i_step;
    } else {
        i_plus = 3*i_step;
        i_minus = 3*(i_step - 1);   
    }

    for (z = 0; z < (Nz) * 2 * dimR * (Nx / 2); ++ z) {
        Re(dCdt[z]) = (1.0/dt) * (Re(MC[i_plus][0][0][0][z])-Re(MC[i_minus][0][0][0][z]));
        Im(dCdt[z]) = (1.0/dt) * (Im(MC[i_plus][0][0][0][z])-Im(MC[i_minus][0][0][0][z]));
    }
    // find physical derivative, sensitivity
    spec2phys(dCdt,dudt);
    spec2phys(IC_given,v);
    
    //compute numerator and denominator of projection
    num = 0.0;
    den = 0.0;
    XFLOW = 0 * (3 * Nx / 2) * qpts * (3 * Nz / 2);
    YFLOW = 1 * (3 * Nx / 2) * qpts * (3 * Nz / 2);
    ZFLOW = 2 * (3 * Nx / 2) * qpts * (3 * Nz / 2);
    
    for (y = 0; y < qpts; ++y){
        for (z = 0; z < (3 * Nz / 2); ++z){
            for (x = 0; x < (3 * Nx / 2); ++x){
                i = (x * qpts + y) * (3 * Nz / 2) + z;
                num += dudt[i + XFLOW] * v[i + XFLOW] * W[y];
                den += dudt[i + XFLOW] * dudt[i + XFLOW] * W[y];
                num += dudt[i + YFLOW] * v[i + YFLOW] * W[y];
                den += dudt[i + YFLOW] * dudt[i + YFLOW] * W[y];
                num += dudt[i + ZFLOW] * v[i + ZFLOW] * W[y];
                den += dudt[i + ZFLOW] * dudt[i + ZFLOW] * W[y];
            }
        }
    }
    
    // subtract component of IC_given that is parallel to dCdt


    for (z = 0; z < (Nz) * 2 * dimR * (Nx / 2); ++ z) {
        Re(IC_given[z]) -= (num/den) * Re(dCdt[z]);
        Im(IC_given[z]) -= (num/den) * Im(dCdt[z]);
    }


    // free memory
    freecVector(dCdt);
    freedVector(dudt);
    freedVector(v);
	return (num/den);
}



/***************************************************************
*                                                              *
*                           SOLUTION IO                        *
*                                                              *
****************************************************************/

// expose restart2 to python
void read_solution(char * filename, mcomplex * C_ptr)
{
    extern int qpts, dimR, dimQ, Nx, Nz;
    extern mcomplex ****U, ****C;
    extern mcomplex ****IU, ****IC;
    int x, y, z;
    hid_t file_id1, dataset_a, dataset_b;       /* file identifier */
    hid_t dataset_ia, dataset_ib;
    hid_t complex_id;

    /* define compound datatype for the complex number */
    typedef struct {
        double re;              /*real part */
        double im;              /*imaginary part */
    } complex_t;

    complex_id = H5Tcreate(H5T_COMPOUND, sizeof(complex_t));
    H5Tinsert(complex_id, "real", HOFFSET(complex_t, re),
              H5T_NATIVE_DOUBLE);
    H5Tinsert(complex_id, "imaginary", HOFFSET(complex_t, im),
              H5T_NATIVE_DOUBLE);

    /* define some temporal matrix to store the data to the hdf file */
    complex_t Matrix1[dimR][Nz][Nx / 2];
    complex_t Matrix2[dimR][Nz][Nx / 2];

    complex_t IMatrix1[dimR][Nz][Nx / 2];
    complex_t IMatrix2[dimR][Nz][Nx / 2];

    // open the file and dataset 
    file_id1 = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    dataset_a = H5Dopen1(file_id1, "/data_alpha");
    dataset_b = H5Dopen1(file_id1, "/data_beta");

    dataset_ia = H5Dopen1(file_id1, "/data_ialpha");
    dataset_ib = H5Dopen1(file_id1, "/data_ibeta");

    assert(EXIT_SUCCESS ==
        H5Dread(dataset_a, complex_id, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                Matrix1));
    assert(EXIT_SUCCESS ==
        H5Dread(dataset_b, complex_id, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                Matrix2));

    assert(EXIT_SUCCESS ==
        H5Dread(dataset_ia, complex_id, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                IMatrix1));
    assert(EXIT_SUCCESS ==
        H5Dread(dataset_ib, complex_id, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                IMatrix2));


    assert(EXIT_SUCCESS == H5Dclose(dataset_a));
    assert(EXIT_SUCCESS == H5Dclose(dataset_b));

    assert(EXIT_SUCCESS == H5Dclose(dataset_ia));
    assert(EXIT_SUCCESS == H5Dclose(dataset_ib));

    assert(EXIT_SUCCESS == H5Fclose(file_id1));

    for (y = 0; y < dimR; y++) {
        for (z = 0; z < Nz; ++z) {
            for (x = 0; x < Nx / 2; ++x) {
                Re(C[z][ALPHA][y][x]) = Matrix1[y][z][x].re;
                Im(C[z][ALPHA][y][x]) = Matrix1[y][z][x].im;
                Re(C[z][BETA][y][x]) = Matrix2[y][z][x].re;
                Im(C[z][BETA][y][x]) = Matrix2[y][z][x].im;

                Re(IC[z][ALPHA][y][x]) = IMatrix1[y][z][x].re;
                Im(IC[z][ALPHA][y][x]) = IMatrix1[y][z][x].im;
                Re(IC[z][BETA][y][x]) = IMatrix2[y][z][x].re;
                Im(IC[z][BETA][y][x]) = IMatrix2[y][z][x].im;
            }
        }
    }

    /* compute  ux hat, uy hat, uz hat given alpha, beta,. */
    initAlphaBeta2();

    /* compute the boundary condition for previous stage or time step, 
       result is stored in Uxb and Uzb */
    if (increBoundary() != NO_ERR) {
        printf("increBoundary failure\n");
    }


    /* compute iux, iuy, iuz given i-alpha, i-beta */

    incre_initAlphaBeta2();

    memset(U[Nz / 2][0][0], 0, 5 * qpts * (Nx / 2) * sizeof(mcomplex));
    memset(IU[Nz / 2][0][0], 0, 5 * qpts * (Nx / 2) * sizeof(mcomplex));

    memmove(C_ptr, C[0][0][0], Nz * 2 * dimR * (Nx / 2) * sizeof(mcomplex));
}


// expose write_data2 to python
void save_solution(char * filename, mcomplex * C_ptr)
{
    hid_t file_id1, dataset_a, dataset_b;       /* file identifier */
    hid_t fid1, dataset_ia, dataset_ib;
    hid_t fid2, dataset_Nx, dataset_Ny, dataset_Nz, dataset_dt, dataset_Re,
        dataset_mpg;
    hid_t complex_id;
    hsize_t fdim[] = { dimR, Nz, Nx / 2 };
    hsize_t fdim2[] = { 1 };

    extern mcomplex ****C, ****IC;
    extern int qpts, Nx, Nz, dimR, dimQ;
    extern double dt, re, mpg;
    int x, y, z;
    int Ny = qpts * 2 / 3;
    /* define compound datatype for the complex number */
    typedef struct {
        double re;              /*real part */
        double im;              /*imaginary part */
    } complex_t;

    memmove(C[0][0][0], C_ptr, Nz * 2 * dimR * (Nx / 2) * sizeof(mcomplex));

    /* Compute U from C */
    initAlphaBeta2();

    complex_id = H5Tcreate(H5T_COMPOUND, sizeof(complex_t));
    H5Tinsert(complex_id, "real", HOFFSET(complex_t, re),
              H5T_NATIVE_DOUBLE);
    H5Tinsert(complex_id, "imaginary", HOFFSET(complex_t, im),
              H5T_NATIVE_DOUBLE);

    /* define some temporal matrix to store the data to the hdf file */
    complex_t Matrix1[dimR][Nz][Nx / 2];
    complex_t Matrix2[dimR][Nz][Nx / 2];

    complex_t IMatrix1[dimR][Nz][Nx / 2];
    complex_t IMatrix2[dimR][Nz][Nx / 2];

    for (y = 0; y < dimR; y++) {
        for (z = 0; z < Nz; ++z) {
            for (x = 0; x < Nx / 2; ++x) {
                Matrix1[y][z][x].re = Re(C[z][ALPHA][y][x]);    /* storing solved solutions to Matrix */
                Matrix1[y][z][x].im = Im(C[z][ALPHA][y][x]);
                Matrix2[y][z][x].re = Re(C[z][BETA][y][x]);
                Matrix2[y][z][x].im = Im(C[z][BETA][y][x]);

                IMatrix1[y][z][x].re = Re(IC[z][ALPHA][y][x]);  /* storing solved solutions to Matrix */
                IMatrix1[y][z][x].im = Im(IC[z][ALPHA][y][x]);
                IMatrix2[y][z][x].re = Re(IC[z][BETA][y][x]);
                IMatrix2[y][z][x].im = Im(IC[z][BETA][y][x]);
            }
        }
    }

    /* create File data.h5 for three data sets */
    file_id1 =
        H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    /* create the dataspace */
    fid1 = H5Screate_simple(3, fdim, NULL);
    fid2 = H5Screate_simple(1, fdim2, NULL);

    /* create datasets with name u v w */
    dataset_a =
        H5Dcreate1(file_id1, "/data_alpha", complex_id, fid1, H5P_DEFAULT);
    dataset_b =
        H5Dcreate1(file_id1, "/data_beta", complex_id, fid1, H5P_DEFAULT);

    dataset_ia =
        H5Dcreate1(file_id1, "/data_ialpha", complex_id, fid1,
                   H5P_DEFAULT);
    dataset_ib =
        H5Dcreate1(file_id1, "/data_ibeta", complex_id, fid1, H5P_DEFAULT);


    dataset_Nx =
        H5Dcreate1(file_id1, "data_Nx", H5T_STD_I32LE, fid2, H5P_DEFAULT);
    dataset_Ny =
        H5Dcreate1(file_id1, "data_Ny", H5T_STD_I32LE, fid2, H5P_DEFAULT);
    dataset_Nz =
        H5Dcreate1(file_id1, "data_Nz", H5T_STD_I32LE, fid2, H5P_DEFAULT);

    dataset_dt =
        H5Dcreate1(file_id1, "data_dt", H5T_IEEE_F64LE, fid2, H5P_DEFAULT);
    dataset_Re =
        H5Dcreate1(file_id1, "data_Re", H5T_IEEE_F64LE, fid2, H5P_DEFAULT);
    dataset_mpg =
        H5Dcreate1(file_id1, "data_mpg", H5T_IEEE_F64LE, fid2,
                   H5P_DEFAULT);

    /* write data to corresponding datasets */
    assert(EXIT_SUCCESS ==
        H5Dwrite(dataset_a, complex_id, H5S_ALL, fid1, H5P_DEFAULT,
                 Matrix1));
    assert(EXIT_SUCCESS ==
        H5Dwrite(dataset_b, complex_id, H5S_ALL, fid1, H5P_DEFAULT,
                 Matrix2));
    assert(EXIT_SUCCESS ==
        H5Dwrite(dataset_ia, complex_id, H5S_ALL, fid1, H5P_DEFAULT,
                 IMatrix1));
    assert(EXIT_SUCCESS ==
        H5Dwrite(dataset_ib, complex_id, H5S_ALL, fid1, H5P_DEFAULT,
                 IMatrix2));

    assert(EXIT_SUCCESS ==
        H5Dwrite(dataset_Nx, H5T_NATIVE_INT, H5S_ALL, fid2, H5P_DEFAULT,
                 &Nx));
    assert(EXIT_SUCCESS ==
        H5Dwrite(dataset_Ny, H5T_NATIVE_INT, H5S_ALL, fid2, H5P_DEFAULT,
                 &Ny));
    assert(EXIT_SUCCESS ==
        H5Dwrite(dataset_Nz, H5T_NATIVE_INT, H5S_ALL, fid2, H5P_DEFAULT,
                 &Nz));
    assert(EXIT_SUCCESS ==
        H5Dwrite(dataset_dt, H5T_IEEE_F64LE, H5S_ALL, fid2, H5P_DEFAULT,
                 &dt));
    assert(EXIT_SUCCESS ==
        H5Dwrite(dataset_Re, H5T_IEEE_F64LE, H5S_ALL, fid2, H5P_DEFAULT,
                 &re));
    assert(EXIT_SUCCESS ==
        H5Dwrite(dataset_mpg, H5T_IEEE_F64LE, H5S_ALL, fid2, H5P_DEFAULT,
                 &mpg));

    /* close datasets and file */
    assert(EXIT_SUCCESS == H5Dclose(dataset_a));
    assert(EXIT_SUCCESS == H5Dclose(dataset_b));
    assert(EXIT_SUCCESS == H5Dclose(dataset_ia));
    assert(EXIT_SUCCESS == H5Dclose(dataset_ib));
    assert(EXIT_SUCCESS == H5Dclose(dataset_Nx));
    assert(EXIT_SUCCESS == H5Dclose(dataset_Ny));
    assert(EXIT_SUCCESS == H5Dclose(dataset_Nz));
    assert(EXIT_SUCCESS == H5Dclose(dataset_dt));
    assert(EXIT_SUCCESS == H5Dclose(dataset_Re));
    assert(EXIT_SUCCESS == H5Dclose(dataset_mpg));

    assert(EXIT_SUCCESS == H5Sclose(fid1));
    assert(EXIT_SUCCESS == H5Sclose(fid2));
    assert(EXIT_SUCCESS == H5Fclose(file_id1));
}

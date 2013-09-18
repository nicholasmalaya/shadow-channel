/**********************************************************************
main code used to solve State, incremental state, adjoint and incremental adjoint equations.

update: 03/09, working on checking the accuracy of new rk scheme for state equations.
**********************************************************************/

#include "main.h"
#include "hdf5.h"

int MAXSTEP = 3;

int main(int argc, char **argv)
{
  /*****************************Definition of all variables ****************************/
    /* External Variables.  All external variables are defined in main.h */
    extern int qpts, dimR, dimQ, Nx, Nz;
    extern double dt, re, mpg;

    extern double *Kx, *Kz, **K2, *cfl2;
    extern double **Q, **Qp, **Qpp, **R, **Rp, **Qw, **Qpw, **Rw, **Qs,
        **Qps, **Qpps, **Rs, **Rps, *Rp0, **Rpw, **Qppw, *Rpp0;

    extern double *Uadd, *Vadd, *Vpadd;
    extern double *Qy;
    extern double *W;

    extern mcomplex ****U, ****C;       /* state variables */
    extern mcomplex **Fa, **Fb, **TM;
    extern mcomplex *fa, *fb, *tm;
    extern double **MZ;
    extern double ***M;

    extern mcomplex ****IU, ****IC;     /* incremental state variables */
    extern mcomplex **IFa, **IFb, **ITM;
    extern mcomplex *Ifa, *Ifb, *Itm;

    extern mcomplex ****AU, ****AC;     /* adjoint variables and will use the same other variables
                                           used in state equations */

    extern mcomplex ****IAU, ****IAC;   /* incremental adjoint variables */

    extern mcomplex **Uxbt, **Uzb;      /* variables used to store dux duz evaluated at y=-1 used for
                                           computing boundary conditions for incremental state equations */
    extern mcomplex **Uxb, **Uzb;       /* variables used to store dux duz evaluated at y=-1 from previous 
                                           state used for boundary conditions for incremental state equations */
    extern mcomplex **IUxb, **IUzb;
    extern mcomplex **IAUxb, **IAUzb;
    extern mcomplex **AUxb, **AUzb;     /* variables used to store dux duz evaluated at y=-1 used for
                                           computing boundary conditions for incremental state equations */
    extern fftw_complex ***CT, ***ICT;  /* variables used in fft */
    extern fftw_plan pf1, pf2;
    extern fftw_plan Ipf1, Ipf2;
    extern rfftwnd_plan pr1, pr2;

    extern mcomplex *****MC, *****MIC;  /* variables used to store state and incremental state
                                           solutions between two check points. */

    extern mcomplex ****MU, ****MIU;    /* variables used to store manufacture solutions */
    extern mcomplex ****LU, ****LIU;
    /* Local Variables */
    int n, z, dctr, tsteps;
    int Ny, sizeRealTransform;
    double Lx, Lz;
    fftw_complex *fout = NULL;
    int restart_flag;
    int count;
    int checknum, checkstep;
    /********************************************end of variable definitions ***********/


    /**************************read in and define parameters *********************/
    if (argc != 11) {
        printf
            ("Required arguments are Nx,Ny,Nz,Lx,Lz,dt,tsteps,mpg,Re, restart_flag.\n");
        return (EXIT_FAILURE);
    }

    Nx = atoi(argv[1]);         /* number of terms in the truncated expansion
                                   in x */
    Ny = atoi(argv[2]);         /* number of terms in the truncated expansion
                                   in y */
    Nz = atoi(argv[3]);         /* number of terms in the truncated expansion
                                   in z */
    Lx = atof(argv[4]);         /* length of the interval in x */
    Lz = atof(argv[5]);         /* length of the interval in z */
    dt = atof(argv[6]);         /* value of time step */
    tsteps = atoi(argv[7]);     /* number of time steps to take */
    mpg = atof(argv[8]);        /* mean pressure gradient */
    re = atof(argv[9]);         /* Reynolds number */

    restart_flag = atoi(argv[10]);      /* whether to restart, 0, starting from initial conditions, 
                                           nonzero, reading in stored results on hdf5 and restart from this time step */

    printf("Nx,Ny,Nz,Lx,Lz,dt,tsteps,mpg,Re,restart_flag\n"
           "%d %d %d %f %f %f %d %f %f %d\n",
           Nx,Ny,Nz,Lx,Lz,dt,tsteps,mpg,re,restart_flag);

    re = 1. / re;               /* time step routines assume I pass 1/Re */
    qpts = 3 * Ny / 2;          /* number of quadrature points 
                                   (see page 9 of Moser's notes) */
    dimR = Ny - 2;              /* dimR and dimQ denote the number of terms */
    dimQ = Ny - 4;              /* in the truncated expansions for the */
    /* functions v_hat, g_hat, U, W */
    /* (see page 5 of Moser's notes). */
    sizeRealTransform = 3 * Nx / 2;     /* for the FFTs */

    /********************************** check input parameters, Nx/4, Nz/2 ***************/
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

    /**********************************end of parameter checking ****************/


    /***********************Initialize and allocate all variables ***************/
    /* Compute wave numbers */
    if (waveNums(Nx / 2, Nz, Lx, Lz) != NO_ERR) {
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
        freedVector(Vadd);
        freedVector(Vpadd);
        freedVector(Uadd);
        freedMatrix(Rpw);
        freedMatrix(Qppw);
        freedVector(Qy);
        freedVector(W);
        freedVector(Rpp0);
        return (EXIT_FAILURE);
    }

    /* get memory for 4D arrays and other matrices */
    if (getMem() != NO_ERR) {
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
        freedVector(Vadd);
        freedVector(Vpadd);
        freedVector(Uadd);
        freedMatrix(Rpw);
        freedMatrix(Qppw);
        freedVector(Qy);
        freedVector(W);
        freedVector(Kx);
        freedVector(Kz);
        freedMatrix(K2);
        freedVector(Rpp0);
        return (EXIT_FAILURE);
    }

    /* Create plans for FFTs */
    pf1 = fftw_create_plan_specific(3 * Nz / 2, FFTW_BACKWARD,
                                    FFTW_MEASURE | FFTW_IN_PLACE, CT[0][0],
                                    3 * Nx / 4 + 1, fout, -1);
    pf2 =
        fftw_create_plan_specific(3 * Nz / 2, FFTW_FORWARD,
                                  FFTW_MEASURE | FFTW_IN_PLACE, CT[0][0],
                                  3 * Nx / 4 + 1, fout, -1);
    pr1 =
        rfftwnd_create_plan(1, &sizeRealTransform, FFTW_COMPLEX_TO_REAL,
                            FFTW_MEASURE | FFTW_IN_PLACE);
    pr2 =
        rfftwnd_create_plan(1, &sizeRealTransform, FFTW_REAL_TO_COMPLEX,
                            FFTW_MEASURE | FFTW_IN_PLACE);


    /* Create plans for FFTs */
    Ipf1 = fftw_create_plan_specific(3 * Nz / 2, FFTW_BACKWARD,
                                     FFTW_MEASURE | FFTW_IN_PLACE,
                                     ICT[0][0], 3 * Nx / 4 + 1, fout, -1);
    Ipf2 =
        fftw_create_plan_specific(3 * Nz / 2, FFTW_FORWARD,
                                  FFTW_MEASURE | FFTW_IN_PLACE, ICT[0][0],
                                  3 * Nx / 4 + 1, fout, -1);


    /* set variables for checking CFL condition */
    if (cflVars(Lx, Lz) != 0) {
        printf("Error creating CF variables\n");
        fftw_destroy_plan(pf1);
        fftw_destroy_plan(pf2);
        rfftwnd_destroy_plan(pr1);
        rfftwnd_destroy_plan(pr2);
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
        freedVector(Vadd);
        freedVector(Vpadd);
        freedVector(Uadd);
        freedMatrix(Rpw);
        freedMatrix(Qppw);
        freedVector(Qy);
        freedVector(W);
        freedVector(Kx);
        freedVector(Kz);
        freedMatrix(K2);
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
        freedVector(Rpp0);
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
           (MAXSTEP * 3 + 1) * (Nz) * 2 * (Ny -
                                           2) * (Nx / 2) *
           sizeof(mcomplex));
    memset(MIC[0][0][0][0], 0,
           (MAXSTEP * 3 + 1) * (Nz) * 2 * (Ny -
                                           2) * (Nx / 2) *
           sizeof(mcomplex));

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
    }

    /**********************************end of restart **********************/

    /* store current Fourier coefficient into MC and MIC  used for solving the adjoint system. */
    memcpy(MC[0][0][0][0], C[0][0][0],
           (Nz) * 2 * (Ny - 2) * (Nx / 2) * sizeof(mcomplex));
    memcpy(MIC[0][0][0][0], IC[0][0][0],
           (Nz) * 2 * (Ny - 2) * (Nx / 2) * sizeof(mcomplex));
    count = 0;


    /************************solving state and incremental state equations ****/

    /* time step for forward problem */
    for (n = restart_flag; n < tsteps; ++n) {   /* loop for each timestep */
        printf("Step %d/%d\n", n, tsteps);
        for (dctr = 0; dctr < 3; ++dctr) {      /* RK steps */

            count = count + 1;

            /* copy the result to Uxbt, Uzbt. Uxb and Uzb will be used later for boundary condition
               of current time stage */
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
                n = tsteps;
                break;
            }

            project0(count, dctr, n, NULL);     /* (kx, kz)=(0, 0), solve for a, b */
            project(count, n, dctr, 0, 1, NULL);        /* (kx, kz)!=(0, 0), solve for alpha, beta */
            for (z = 1; z < Nz; ++z) {
                if (z == Nz / 2) {
                    memset(U[z][0][0], 0,
                           5 * qpts * (Nx / 2) * sizeof(mcomplex));
                    continue;
                }

                project(count, n, dctr, z, 0, NULL);
            }

            /*now update the boundary condition using current time step solution of state equation */
            if (increBoundary() != NO_ERR) {
                printf("increBoundary failure\n");
                n = tsteps;
                break;
            }

            increproject0(count, dctr, n, 1, NULL);
            increproject(count, dctr, 0, 1, n, NULL);
            for (z = 1; z < Nz; ++z) {
                if (z == Nz / 2) {
                    // SET U[z][XEL,YEL,ZEL,DXEL,DZEL] TO ZEROS 
                    memset(IU[z][0][0], 0,
                           5 * qpts * (Nx / 2) * sizeof(mcomplex));
                    continue;
                }

                increproject(count, dctr, z, 0, n, NULL);
            }
            memcpy(MC[count][0][0][0], C[0][0][0],
                   (Nz) * 2 * (Ny - 2) * (Nx / 2) * sizeof(mcomplex));
            memcpy(MIC[count][0][0][0], IC[0][0][0],
                   (Nz) * 2 * (Ny - 2) * (Nx / 2) * sizeof(mcomplex));

        }                       /* end for dctr... */

        /* now writing the results to HDF file at selected time steps */
        if (((n + 1) % 100 == 0) && (n + 1 < tsteps)) { /* when we need to store the current time step results */
            write_data2(n);
            checknum = checknum + 1;
            count = 0;
            memcpy(MC[count][0][0][0], C[0][0][0],
                   (Nz) * 2 * (Ny - 2) * (Nx / 2) * sizeof(mcomplex));
            memcpy(MIC[count][0][0][0], IC[0][0][0],
                   (Nz) * 2 * (Ny - 2) * (Nx / 2) * sizeof(mcomplex));
        }


    }                           /* end for n... */

    /****************************end of solving forward systems. *******************/

    return (EXIT_SUCCESS);


    /*******************************start of solving backward systems *************/
    for (checkstep = checknum; checkstep >= 0; checkstep--) {

        restart_flag = checkstep * MAXSTEP;
        printf("\n");
        printf(" restarting from time step: %d\n", restart_flag);
        restart2(restart_flag);
        count = 0;


        memcpy(MC[count][0][0][0], C[0][0][0],
               (Nz) * 2 * (Ny - 2) * (Nx / 2) * sizeof(mcomplex));
        memcpy(MIC[count][0][0][0], IC[0][0][0],
               (Nz) * 2 * (Ny - 2) * (Nx / 2) * sizeof(mcomplex));
        for (n = restart_flag; n < restart_flag + MAXSTEP; ++n) {       /* loop for each timestep */
            for (dctr = 0; dctr < 3; ++dctr) {  /* RK steps */

                count = count + 1;

                /* copy the result to Uxbt, Uzbt. Uxb and Uzb will be used later for boundary condition
                   of current time stage */
                memcpy(Uxbt[0], Uxb[0],
                       (Nz) * (Nx / 2) * sizeof(fftw_complex));
                memcpy(Uzbt[0], Uzb[0],
                       (Nz) * (Nx / 2) * sizeof(fftw_complex));
                memset(Uxb[0], 0, Nz * (Nx / 2) * sizeof(fftw_complex));
                memset(Uzb[0], 0, Nz * (Nx / 2) * sizeof(fftw_complex));

                if (pass1(dctr, n) != NO_ERR) {
                    printf("Pass1 failure\n");
                    checkstep = -1;
                    break;
                }

                project0(count, dctr, n, NULL); /* (kx, kz)=(0, 0), solve for a, b */
                project(count, n, dctr, 0, 1, NULL);    /* (kx, kz)!=(0, 0), solve for alpha, beta */
                for (z = 1; z < Nz; ++z) {
                    if (z == Nz / 2) {
                        memset(U[z][0][0], 0,
                               5 * qpts * (Nx / 2) * sizeof(mcomplex));
                        continue;
                    }

                    project(count, n, dctr, z, 0, NULL);
                }

                /*now update the boundary condition using current time step solution of state equation */
                if (increBoundary() != NO_ERR) {
                    printf("increBoundary failure\n");
                    checkstep = -1;
                    break;
                }

                increproject0(count, dctr, n, 1, NULL);
                increproject(count, dctr, 0, 1, n, NULL);
                for (z = 1; z < Nz; ++z) {
                    if (z == Nz / 2) {
                        // SET U[z][XEL,YEL,ZEL,DXEL,DZEL] TO ZEROS 
                        memset(IU[z][0][0], 0,
                               5 * qpts * (Nx / 2) * sizeof(mcomplex));
                        continue;
                    }

                    increproject(count, dctr, z, 0, n, NULL);
                }

            }                   /* end for dctr... */
        }

        memcpy(Uxb[0], AUxb[0], (Nz) * (Nx / 2) * sizeof(fftw_complex));
        memcpy(Uzb[0], AUzb[0], (Nz) * (Nx / 2) * sizeof(fftw_complex));

        /* time step for backward equations */
        for (n = restart_flag + MAXSTEP; n > restart_flag; --n) {       /* loop for each timestep */
            for (dctr = 0; dctr < 3; ++dctr) {  /* RK steps */

                /* copy the result to Uxbt, Uzbt. Uxb and Uzb will be used later for boundary condition
                   of current time stage */
                memcpy(Uxbt[0], Uxb[0],
                       (Nz) * (Nx / 2) * sizeof(fftw_complex));
                memcpy(Uzbt[0], Uzb[0],
                       (Nz) * (Nx / 2) * sizeof(fftw_complex));
                memset(Uxb[0], 0, Nz * (Nx / 2) * sizeof(fftw_complex));
                memset(Uzb[0], 0, Nz * (Nx / 2) * sizeof(fftw_complex));
                memset(AUxb[0], 0, Nz * (Nx / 2) * sizeof(fftw_complex));
                memset(AUzb[0], 0, Nz * (Nx / 2) * sizeof(fftw_complex));

                count = count - 1;

                /*read data from memery */
                if (dctr < 2) {
                    memcpy(C[0][0][0], MC[count - 1][0][0][0],
                           (Nz) * 2 * (Ny -
                                       2) * (Nx / 2) * sizeof(mcomplex));
                    memcpy(IC[0][0][0], MIC[count - 1][0][0][0],
                           (Nz) * 2 * (Ny -
                                       2) * (Nx / 2) * sizeof(mcomplex));


                    /*reconstruct the state and incremental state solution u, iu from alpha and beta */
                    initAlphaBeta2();
                    if (increBoundary() != NO_ERR) {
                        printf("increBoundary failure\n");
                    }
                    incre_initAlphaBeta2();

                    if (pass2(dctr, n) != NO_ERR) {
                        printf("Pass2 failure\n");
                        n = restart_flag - 1;
                        break;
                    }
                    memcpy(LU[0][0][0], U[0][0][0],
                           (Nz) * 5 * qpts * (Nx / 2) * sizeof(mcomplex));
                    memcpy(LIU[0][0][0], IU[0][0][0],
                           (Nz) * 5 * qpts * (Nx / 2) * sizeof(mcomplex));

                }
                /*read data from memery */
                memcpy(C[0][0][0], MC[count][0][0][0],
                       (Nz) * 2 * (Ny - 2) * (Nx / 2) * sizeof(mcomplex));
                memcpy(IC[0][0][0], MIC[count][0][0][0],
                       (Nz) * 2 * (Ny - 2) * (Nx / 2) * sizeof(mcomplex));


                /*reconstruct the state and incremental state solution u, iu from alpha and beta */
                initAlphaBeta2();
                if (increBoundary() != NO_ERR) {
                    printf("increBoundary failure\n");
                }
                incre_initAlphaBeta2();

                memcpy(AUxb[0], Uxb[0],
                       (Nz) * (Nx / 2) * sizeof(fftw_complex));
                memcpy(AUzb[0], Uzb[0],
                       (Nz) * (Nx / 2) * sizeof(fftw_complex));

                if (pass2(dctr, n) != NO_ERR) {
                    printf("Pass2 failure\n");
                    n = restart_flag - 1;
                    break;
                }

                memset(Uxb[0], 0, Nz * (Nx / 2) * sizeof(mcomplex));
                memset(Uzb[0], 0, Nz * (Nx / 2) * sizeof(mcomplex));

                adjproject0(dctr, n, count, NULL);
                adjproject(n, dctr, 0, 1, count, NULL);
                for (z = 1; z < Nz; ++z) {
                    if (z == Nz / 2) {
                        memset(AU[z][0][0], 0,
                               5 * qpts * (Nx / 2) * sizeof(mcomplex));
                        continue;
                    }

                    adjproject(n, dctr, z, 0, count, NULL);
                }

                /*now update the boundary condition using current time step solution of state equation */
                if (increBoundary() != NO_ERR) {
                    printf("increBoundary failure\n");
                    n = restart_flag - 1;
                    break;
                }

                increadjproject0(dctr, n, count, NULL);
                increadjproject(n, dctr, 0, 1, count, NULL);
                for (z = 1; z < Nz; ++z) {
                    if (z == Nz / 2) {
                        // SET U[z][XEL,YEL,ZEL,DXEL,DZEL] TO ZEROS 
                        memset(IAU[z][0][0], 0,
                               5 * qpts * (Nx / 2) * sizeof(mcomplex));
                        continue;
                    }

                    increadjproject(n, dctr, z, 0, count, NULL);
                }
                if (n == restart_flag + 1 && dctr == 2) {
                    memcpy(AUxb[0], Uxb[0],
                           (Nz) * (Nx / 2) * sizeof(fftw_complex));
                    memcpy(AUzb[0], Uzb[0],
                           (Nz) * (Nx / 2) * sizeof(fftw_complex));
                }
            }
        }
    }

    /* clean up... */
    fftw_destroy_plan(pf1);
    fftw_destroy_plan(pf2);
    rfftwnd_destroy_plan(pr1);
    rfftwnd_destroy_plan(pr2);
    fftw_destroy_plan(Ipf1);
    fftw_destroy_plan(Ipf2);

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
    freedVector(Kx);
    freedVector(Kz);
    freedMatrix(K2);
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
    freedVector(Rp0);
    freec3Darray(ICT);
    freec4Darray(IU);
    freec4Darray(IC);
    freecMatrix(IFa);
    freecMatrix(IFb);
    freecMatrix(ITM);
    freecVector(Ifa);
    freecVector(Ifb);
    freecVector(Itm);
    freecMatrix(Uxbt);
    freecMatrix(Uzbt);
    freecMatrix(Uxb);
    freecMatrix(Uzb);
    freedVector(Vadd);
    freedVector(Vpadd);
    freedVector(Uadd);
    freec4Darray(AU);
    freec4Darray(AC);
    freec4Darray(IAU);
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
    freedVector(Rpp0);
    return (EXIT_SUCCESS);

}                               /* end main */

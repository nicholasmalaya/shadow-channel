#include "main.h"
#include "hdf5.h"

int main(int argc, char **argv)
{
    /* External Variables.  All external variables are defined in main.h */
    extern int qpts, dimR, dimQ, Nx, Nz;
    extern double dt, re, mpg;
    extern double *Kx, *Kz, **K2, *cfl2;
    extern double **Q, **Qp, **Qpp, **R, **Rp, **Qw, **Qpw, **Rw, **Qs,
        **Qps, **Qpps, **Rs, **Rps, *Rp0;
    extern double **MZ;
    extern double ***M;
    extern double *Uadd, *Vadd, *Vpadd;
    extern double *Qy;
    extern double *W;

    extern mcomplex **Fa, **Fb, **TM;
    extern mcomplex *fa, *fb, *tm;
    extern double **MZ;
    extern double ***M;

    extern mcomplex ****IU;             /* incremental state variables */
    extern mcomplex **IFa, **IFb, **ITM;
    extern mcomplex *Ifa, *Ifb, *Itm;

    extern mcomplex ****AC;             /* adjoint variables and will use the same other variables
                                           used in state equations */

    extern mcomplex ****IAU, ****IAC;   /* incremental adjoint variables */

    extern mcomplex **Uxbt, **Uzb;      /* variables used to store dux duz evaluated at y=-1 used for
                                           computing boundary conditions for incremental state equations */
    extern mcomplex **Uxb, **Uzb;       /* variables used to store dux duz evaluated at y=-1 from previous 
                                           state used for boundary conditions for incremental state equations */

    extern fftw_complex ***CT, ***ICT;  /* variables used in fft */
    extern fftw_plan pf1, pf2;
    extern fftw_plan Ipf1, Ipf2;
    extern rfftwnd_plan pr1, pr2;

    extern mcomplex *****MIC;  /* variables used to store state and incremental state
                                           solutions between two check points. */

    extern mcomplex ****MU, ****MIU;    /* variables used to store manufacture solutions */


    extern mcomplex ****LU, ****LIU;

    // static double e[3] = { 1. / 2., 2. / 3, 1. };
    static double g[3] = { 0, 1., 2. / 3. };
    static double ee[3] = { 1. / 3., 1. / 2., 1. };
    // static double ss[3] = { 1., 2. / 3., 1. };
    /* Local Variables */
    int n, z, dctr, tsteps, j, i, x, y;
    double k11 = 2 * M_PI * 1.;
    double Lx, Lz, Ny;
    fftw_complex *fout = NULL;
    double err, err2;
    int sizeRealTransform;
    int restart_flag;
    int count;
    if (argc != 11) {
        printf
            ("Required arguments are Nx,Ny,Nz,Lx,Lz,dt,tsteps,mpg,Re.\n");
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

    re = 1. / re;               /* time step routines assume I pass 1/Re */
    qpts = 3 * Ny / 2;          /* number of quadrature points 
                                   (see page 9 of Moser's notes) */
    dimR = Ny - 2;              /* dimR and dimQ denote the number of terms */
    dimQ = Ny - 4;              /* in the truncated expansions for the */
    /* functions v_hat, g_hat, U, W */
    /* (see page 5 of Moser's notes). */
    sizeRealTransform = 3 * Nx / 2;     /* for the FFTs */


    //assert((Nx%4==0)&&(Nz%2==0));

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
        freedVector(Kx);
        freedVector(Kz);
        freedMatrix(K2);
        freedVector(Rp0);
        freedVector(Vadd);
        freedVector(Vpadd);
        freedMatrix(Rpw);
        freedVector(Uadd);
        freedVector(Qy);
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
        freedVector(Kx);
        freedVector(Kz);
        freedMatrix(K2);
        freec4Darray(U);
        freec4Darray(C);
        freec3Darray(CT);
        freec3Darray(ICT);
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
        freedVector(Vadd);
        freedVector(Vpadd);
        freedVector(Uadd);
        freecMatrix(ITM);
        freecMatrix(IFa);
        freecMatrix(IFb);
        freecMatrix(Uxb);
        freecMatrix(Uzb);
        freecMatrix(Uzbt);
        freecMatrix(Uxbt);
        freecVector(Ifa);
        freecVector(Ifb);
        freecVector(Itm);
        freedVector(Qy);
        freec4Darray(IU);
        freec4Darray(IC);
        freec4Darray(AU);
        freec4Darray(AC);
        freec4Darray(IAU);
        freec4Darray(IAC);
        freec5Darray(MC);
        freec5Darray(MIC);
        freedMatrix(Rpw);
        freec4Darray(LU);
        freec4Darray(LIU);
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

    /* for (i=0; i< (MAXSTEP); i++)
       {
       // Re(MC[i][0][XEL][0][0])=1.;
       Re(MC[3*i][0][XEL][0][0])=1.*sin(2*M_PI*(i*dt));
       Re(MC[3*i+1][0][XEL][0][0])=1.*sin(2*M_PI*(i*dt+dt*0.5));
       Re(MC[3*i+2][0][XEL][0][0])=1.*sin(2*M_PI*(i*dt+dt*2./3));
       } */


    count = 600;


    /* time step for backward equations */
    for (n = tsteps; n > 0; --n) {      /* loop for each timestep */
        for (dctr = 0; dctr < 3; ++dctr) {      /* RK steps */

            /* copy the result to Uxbt, Uzbt. Uxb and Uzb will be used later for boundary condition
               of current time stage */
            memcpy(Uxbt[0], Uxb[0],
                   (Nz) * (Nx / 2) * sizeof(fftw_complex));
            memcpy(Uzbt[0], Uzb[0],
                   (Nz) * (Nx / 2) * sizeof(fftw_complex));
            memset(Uxb[0], 0, Nz * (Nx / 2) * sizeof(fftw_complex));
            memset(Uzb[0], 0, Nz * (Nx / 2) * sizeof(fftw_complex));

            count = count - 1;

            /*read data from memery */
            if (dctr < 2) {
                memcpy(C[0][0][0], MC[count - 1][0][0][0],
                       (Nz) * 2 * (Ny - 2) * (Nx / 2) * sizeof(mcomplex));
                memcpy(IC[0][0][0], MIC[count - 1][0][0][0],
                       (Nz) * 2 * (Ny - 2) * (Nx / 2) * sizeof(mcomplex));


                /*reconstruct the state and incremental state solution u, iu from alpha and beta */
                initAlphaBeta2();
                if (increBoundary() != NO_ERR) {
                    printf("increBoundary failure\n");
                }
                incre_initAlphaBeta2();

                if (pass2(dctr, n) != NO_ERR) {
                    printf("Pass2 failure\n");
                    n = 0;
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

            if (pass2(dctr, n) != NO_ERR) {
                printf("Pass2 failure\n");
                n = 0;
                break;
            }

            memset(Uxb[0], 0, Nz * (Nx / 2) * sizeof(mcomplex));
            memset(Uzb[0], 0, Nz * (Nx / 2) * sizeof(mcomplex));

            adjproject0(dctr, n, count, adjforce0);
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
                n = tsteps;
                break;
            }
            /*  for (z = 0; z < Nz; ++z)
               {
               for (x = 0; x < Nx/2; ++x)
               {
               printf("Uxb[%d][%d]=%f+i%f, Uzb=%f+i%f\n",z, x, Re(Uxb[z][x]), Im(Uxb[z][x]), Re(Uzb[z][x]), Im(Uzb[z][x]));
               }
               } */
            increadjproject0(dctr, n, count, increadjforce0);
            increadjproject(n, dctr, 0, 1, count, increadjforce);
            for (z = 1; z < Nz; ++z) {
                if (z == Nz / 2) {
                    // SET U[z][XEL,YEL,ZEL,DXEL,DZEL] TO ZEROS 
                    memset(AU[z][0][0], 0,
                           5 * qpts * (Nx / 2) * sizeof(mcomplex));
                    continue;
                }

                increadjproject(n, dctr, z, 0, count, increadjforce);
            }
        }

        err = 0;
        err2 = 0;
        for (y = 0; y < qpts; y++) {
            for (z = 0; z < Nz; ++z) {
                for (x = 0; x < Nx / 2; ++x) {

                    if (x == 0 && z == 0) {
                        err2 =
                            (Re(AU[z][XEL][y][x]) -
                             (1 -
                              Qy[y] * Qy[y]) * sin(2 * M_PI * (n * dt -
                                                               dt -
                                                               0.01))) *
                            (Re(AU[z][XEL][y][x]) -
                             (1 -
                              Qy[y] * Qy[y]) * sin((n * dt - dt -
                                                    0.01) * 2 * M_PI)) *
                            W[y]
                            +
                            Im(AU[z][XEL][y][x]) * Im(AU[z][XEL][y][x]) *
                            W[y] +
                            Re(AU[z][YEL][y][x]) * Re(AU[z][YEL][y][x]) *
                            W[y]
                            +
                            Im(AU[z][YEL][y][x]) * Im(AU[z][YEL][y][x]) *
                            W[y] +
                            Re(AU[z][ZEL][y][x]) * Re(AU[z][ZEL][y][x]) *
                            W[y]
                            +
                            Im(AU[z][ZEL][y][x]) * Im(AU[z][ZEL][y][x]) *
                            W[y];
                        err = err2 + err;
                    }

                    else {
                        err2 =
                            Re(AU[z][XEL][y][x]) * Re(AU[z][XEL][y][x]) *
                            W[y]
                            +
                            Im(AU[z][XEL][y][x]) * Im(AU[z][XEL][y][x]) *
                            W[y] +
                            Re(AU[z][YEL][y][x]) * Re(AU[z][YEL][y][x]) *
                            W[y]
                            +
                            Im(AU[z][YEL][y][x]) * Im(AU[z][YEL][y][x]) *
                            W[y] +
                            Re(AU[z][ZEL][y][x]) * Re(AU[z][ZEL][y][x]) *
                            W[y]
                            +
                            Im(AU[z][ZEL][y][x]) * Im(AU[z][ZEL][y][x]) *
                            W[y];
                        err = err2 + err;
                    }

                }
            }
        }
        printf("%d,%d,  %25.16e", n, count, sqrt(err));
        err = 0;
        err2 = 0;
        for (y = 0; y < qpts; y++) {
            for (z = 0; z < Nz; ++z) {
                for (x = 0; x < Nx / 2; ++x) {

                    if (x == 0 && z == 0) {
                        err2 =
                            (Re(IAU[z][XEL][y][x]) -
                             (2 - Qy[y] -
                              Qy[y] * Qy[y]) * sin(2 * M_PI * (n * dt -
                                                               dt -
                                                               0.01))) *
                            (Re(IAU[z][XEL][y][x]) -
                             (2 - Qy[y] -
                              Qy[y] * Qy[y]) * sin((n * dt - dt -
                                                    0.01) * 2 * M_PI)) *
                            W[y]
                            +
                            Im(IAU[z][XEL][y][x]) * Im(IAU[z][XEL][y][x]) *
                            W[y] +
                            Re(IAU[z][YEL][y][x]) * Re(IAU[z][YEL][y][x]) *
                            W[y]
                            +
                            Im(IAU[z][YEL][y][x]) * Im(IAU[z][YEL][y][x]) *
                            W[y] +
                            Re(IAU[z][ZEL][y][x]) * Re(IAU[z][ZEL][y][x]) *
                            W[y]
                            +
                            Im(IAU[z][ZEL][y][x]) * Im(IAU[z][ZEL][y][x]) *
                            W[y];
                        err = err2 + err;
                    } else if (x == 0 && z == 1) {
                        err2 =
                            (Re(IAU[z][XEL][y][x]) -
                             (1 -
                              Qy[y]) * sin(2 * M_PI * (n * dt - dt -
                                                       0.01))) *
                            (Re(IAU[z][XEL][y][x]) -
                             (1 -
                              Qy[y]) * sin((n * dt - dt -
                                            0.01) * 2 * M_PI)) * W[y]
                            +
                            Im(IAU[z][XEL][y][x]) * Im(IAU[z][XEL][y][x]) *
                            W[y] +
                            Re(IAU[z][YEL][y][x]) * Re(IAU[z][YEL][y][x]) *
                            W[y]
                            +
                            Im(AU[z][YEL][y][x]) * Im(IAU[z][YEL][y][x]) *
                            W[y] +
                            Re(IAU[z][ZEL][y][x]) * Re(IAU[z][ZEL][y][x]) *
                            W[y]
                            +
                            Im(IAU[z][ZEL][y][x]) * Im(IAU[z][ZEL][y][x]) *
                            W[y];
                        err = err2 + err;
                    } else if (x == 0 && z == Nz - 1) {
                        err2 =
                            (Re(IAU[z][XEL][y][x]) -
                             (1 -
                              Qy[y]) * sin(2 * M_PI * (n * dt - dt -
                                                       0.01))) *
                            (Re(IAU[z][XEL][y][x]) -
                             (1 -
                              Qy[y]) * sin((n * dt - dt -
                                            0.01) * 2 * M_PI)) * W[y]
                            +
                            Im(IAU[z][XEL][y][x]) * Im(IAU[z][XEL][y][x]) *
                            W[y] +
                            Re(IAU[z][YEL][y][x]) * Re(IAU[z][YEL][y][x]) *
                            W[y]
                            +
                            Im(IAU[z][YEL][y][x]) * Im(IAU[z][YEL][y][x]) *
                            W[y] +
                            Re(IAU[z][ZEL][y][x]) * Re(IAU[z][ZEL][y][x]) *
                            W[y]
                            +
                            Im(IAU[z][ZEL][y][x]) * Im(IAU[z][ZEL][y][x]) *
                            W[y];
                        err = err2 + err;
                    } else {
                        err2 =
                            Re(IAU[z][XEL][y][x]) * Re(IAU[z][XEL][y][x]) *
                            W[y]
                            +
                            Im(IAU[z][XEL][y][x]) * Im(IAU[z][XEL][y][x]) *
                            W[y] +
                            Re(IAU[z][YEL][y][x]) * Re(IAU[z][YEL][y][x]) *
                            W[y]
                            +
                            Im(IAU[z][YEL][y][x]) * Im(IAU[z][YEL][y][x]) *
                            W[y] +
                            Re(IAU[z][ZEL][y][x]) * Re(IAU[z][ZEL][y][x]) *
                            W[y]
                            +
                            Im(IAU[z][ZEL][y][x]) * Im(IAU[z][ZEL][y][x]) *
                            W[y];
                        err = err2 + err;
                    }
                    /*if (err2>1e-12)
                       {
                       printf("%d, %d, %d,  %f+i%f, %f+i%f, %f+i%f, %f\n",x,z,y, Re(IAU[z][XEL][y][x]), Im(IAU[z][XEL][y][x]),
                       Re(IAU[z][YEL][y][x]), Im(IAU[z][YEL][y][x]),Re(IAU[z][ZEL][y][x]), Im(IAU[z][ZEL][y][x]),
                       (1-Qy[y])*sin(2*M_PI*(n*dt-dt-0.01)));
                       } */
                }
            }
        }
        //      printf("%d,%d,  %25.16e\n", n, count, sqrt(err));
        printf("  %25.16e\n", sqrt(err));
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
    return (EXIT_SUCCESS);

}                               /* end main */

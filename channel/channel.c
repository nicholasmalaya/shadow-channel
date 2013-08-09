#include<assert.h>
#include"channel.h"
#include"channel_main.h"

int channel_construct(int _Nx, int _Ny, int _Nz, double _Lx, double _Lz,
                      double _dt, double _mpg, double _re, int _restart_flag)
{
    extern int qpts, dimR, dimQ, Nx, Nz;
    extern double dt, re, mpg;

    Nx = _Nx; Nz = _Nz; dt = _dt;
    mpg = _mpg; re = _re;

    re = 1. / re;               /* time step routines assume I pass 1/Re */
    qpts = 3 * _Ny / 2;         /* number of quadrature points 
                                   (see page 9 of Moser's notes) */
    dimR = _Ny - 2;             /* dimR and dimQ denote the number of terms */
    dimQ = _Ny - 4;             /* in the truncated expansions for the */
    /* functions v_hat, g_hat, U, W */
    /* (see page 5 of Moser's notes). */
    int sizeRealTransform = 3 * Nx / 2;     /* for the FFTs */

    /****************** check input parameters, Nx/4, Nz/2 ***************/
    if (Nx % 4 != 0) {
        printf("Required arguments Ny/4==0\n");
        return (EXIT_FAILURE);
    }
    if (Nz % 2 != 0) {
        printf("Required arguments Nz/2==0\n");
        return (EXIT_FAILURE);
    }
    if (_Ny - 4 < 0) {
        printf("Required arguments Nzy>4\n");
        return (EXIT_FAILURE);
    }
    /* Create matrices using Legendre polynomials */
    if (LegendreSetup() != NO_ERR) {
        return (EXIT_FAILURE);
    }

    /*****************end of parameter checking ****************/

    /**************Initialize and allocate all variables ***************/
    extern double **Q, **Qp, **Qpp, **R, **Rp, **Qw, **Qpw, **Rw, **Qs,
        **Qps, **Qpps, **Rs, **Rps, *Rp0, **Rpw, **Qppw, *Rpp0;

    /* Compute wave numbers */
    if (waveNums(Nx / 2, Nz, _Lx, _Lz) != NO_ERR) {
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

    extern double *Kx, *Kz, **K2, *cfl2;

    extern mcomplex *****MC, *****MIC;  /* variables used to store state
                                           and incremental state solutions
                                           between two check points. */

    extern mcomplex ****MU, ****MIU;    /* variables used to store manufacture
                                           solutions */
    extern mcomplex ****LU, ****LIU;
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
    extern fftw_complex ***CT, ***ICT;  /* variables used in fft */
    extern fftw_plan pf1, pf2;
    extern fftw_plan Ipf1, Ipf2;
    extern rfftwnd_plan pr1, pr2;

    fftw_complex *fout = NULL;

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
    if (cflVars(_Lx, _Lz) != 0) {
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
           (Nz) * 2 * (_Ny - 2) * (Nx / 2) * sizeof(mcomplex));
    memset(IC[0][0][0], 0,
           (Nz) * 2 * (_Ny - 2) * (Nx / 2) * sizeof(mcomplex));
    memset(AC[0][0][0], 0,
           (Nz) * 2 * (_Ny - 2) * (Nx / 2) * sizeof(mcomplex));
    memset(IAC[0][0][0], 0,
           (Nz) * 2 * (_Ny - 2) * (Nx / 2) * sizeof(mcomplex));

    memset(MC[0][0][0][0], 0,
           (MAXSTEP * 3 + 1) * (Nz) * 2 * (_Ny - 2) * (Nx / 2) *
           sizeof(mcomplex));
    memset(MIC[0][0][0][0], 0,
           (MAXSTEP * 3 + 1) * (Nz) * 2 * (_Ny - 2) * (Nx / 2) *
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
    if (_restart_flag != 0) {
        restart2(_restart_flag);
    }

    /**********************************end of restart **********************/

    /* store current Fourier coefficient into MC and MIC
       used for solving the adjoint system. */
    memcpy(MC[0][0][0][0], C[0][0][0], channel_size());
    memcpy(MIC[0][0][0][0], IC[0][0][0], channel_size());

    /* double/complex compiler check */
    double c_double[2] = {1.0, 2.0};
    mcomplex * c_mcomplex = (mcomplex *) c_double;
    assert( Im(*c_mcomplex) == 1.0 );
    assert( Re(*c_mcomplex) == 2.0 );

    return 0;
}


int channel_size()
{
    int Ny = dimR + 2;
    return (Nz) * 2 * (Ny - 2) * (Nx / 2) * sizeof(mcomplex);
}


int channel_forward(double * uTan, double T0, double T1)
{
    return 0;
}


int channel_backward(double * uAdj, double T1, double T0)
{
    return 0;
}

int channel_ddt(double * dudt)
{
    return 0;
}

void channel_destroy()
{}

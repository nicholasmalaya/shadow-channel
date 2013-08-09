#include<assert.h>
#include<math.h>
#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<cblas.h>

#include"kuramoto.h"

int N_GRID=0, N_CHUNK=0, N_STEP=0;
double C_CONST=0;
double DT_STEP=0;
double *** SOLN_U=0;
double *** SOLN_V=0;

// Runge Kutta coefficients
double RK[3][3] = {{ 1./2,  0,    0 },
                   {-1./6,  1./3, 0 },
                   { 0,    -2./3, 1.}};



void
ddt(const double * u, double * dudt)
{
    double dx = 128. / (N_GRID + 1);
    for (int i = 0; i < N_GRID; ++i)
    {
        double up = (i < N_GRID - 1) ? u[i+1] : 0;
        double upp = (i < N_GRID - 2) ? u[i+2] : u[i];
        double um = (i > 0) ? u[i-1] : 0;
        double umm = (i > 1) ? u[i-2] : u[i];
        double ux = (up - um) / (2 * dx);
        double u2x = (up*up - um*um) / (2 * dx);
        double uxx = (up + um - 2 * u[i]) / (dx * dx);
        double uxxp = (upp + u[i] - 2 * up) / (dx * dx);
        double uxxm = (umm + u[i] - 2 * um) / (dx * dx);
        double uxxxx = (uxxp + uxxm - 2 * uxx) / (dx * dx);
        dudt[i] = -C_CONST * ux - 0.5 * u2x - uxx - uxxxx;
    }
}


void
ddtTan(const double * u, const double * v, double * dvdt, int inhomo)
{
    double dx = 128. / (N_GRID + 1);
    for (int i = 0; i < N_GRID; ++i)
    {
        double up = (i < N_GRID - 1) ? u[i+1] : 0;
        double um = (i > 0) ? u[i-1] : 0;
        double ux = (up - um) / (2 * dx);
        double vp = (i < N_GRID - 1) ? v[i+1] : 0;
        double vpp = (i < N_GRID - 2) ? v[i+2] : v[i];
        double vm = (i > 0) ? v[i-1] : 0;
        double vmm = (i > 1) ? v[i-2] : v[i];
        double vx = (vp - vm) / (2 * dx);
        double v2x = (up*vp - um*vm) / (2 * dx);
        double vxx = (vp + vm - 2 * v[i]) / (dx * dx);
        double vxxp = (vpp + v[i] - 2 * vp) / (dx * dx);
        double vxxm = (vmm + v[i] - 2 * vm) / (dx * dx);
        double vxxxx = (vxxp + vxxm - 2 * vxx) / (dx * dx);
        dvdt[i] = -C_CONST * vx - 0.5 * v2x - vxx - vxxxx - inhomo * ux;
    }
}


void
ddtAdj(const double * u, const double * w, double * dwdt)
{
    double dx = 128. / (N_GRID + 1);
    for (int i = 0; i < N_GRID; ++i) {
        double wp = (i < N_GRID - 1) ? w[i+1] : 0;
        double wpp = (i < N_GRID - 2) ? w[i+2] : w[i];
        double wm = (i > 0) ? w[i-1] : 0;
        double wmm = (i > 1) ? w[i-2] : w[i];
        double wx = (wp - wm) / (2 * dx);
        double w2x = u[i] * (wp - wm) / (2 * dx);
        double wxx = (wp + wm - 2 * w[i]) / (dx * dx);
        double wxxp = (wpp + w[i] - 2 * wp) / (dx * dx);
        double wxxm = (wmm + w[i] - 2 * wm) / (dx * dx);
        double wxxxx = (wxxp + wxxm - 2 * wxx) / (dx * dx);
        dwdt[i] = -C_CONST * wx - 0.5 * w2x + wxx + wxxxx;
    }
}


void
stepPrimal(const double * u0, double * u, double dt)
{
    double dudt0[N_GRID], dudt1[N_GRID], *dudt2 = dudt0;

    memmove(u, u0, sizeof(double) * N_GRID);

    ddt(u, dudt0);
    cblas_daxpy(N_GRID, dt * RK[0][0], dudt0, 1, u, 1);

    ddt(u, dudt1);
    cblas_daxpy(N_GRID, dt * RK[1][0], dudt0, 1, u, 1);
    cblas_daxpy(N_GRID, dt * RK[1][1], dudt1, 1, u, 1);

    ddt(u, dudt2);
    cblas_daxpy(N_GRID, dt * RK[2][1], dudt1, 1, u, 1);
    cblas_daxpy(N_GRID, dt * RK[2][2], dudt2, 1, u, 1);
}


void
stepTangent(const double * u0, const double * v0, double * v, double dt,
            int inhomo)
{
    double u1[N_GRID], u2[N_GRID], dudt0[N_GRID], dudt1[N_GRID];
    double * dvdt0 = dudt0, * dvdt1 = dudt1, * dvdt2 = dudt0;

    // Primal
    ddt(u0, dudt0);
    memmove(u1, u0, sizeof(double) * N_GRID);
    cblas_daxpy(N_GRID, dt * RK[0][0], dudt0, 1, u1, 1);

    ddt(u1, dudt1);
    memmove(u2, u1, sizeof(double) * N_GRID);
    cblas_daxpy(N_GRID, dt * RK[1][0], dudt0, 1, u2, 1);
    cblas_daxpy(N_GRID, dt * RK[1][1], dudt1, 1, u2, 1);

    // Tagent
    memmove(v, v0, sizeof(double) * N_GRID);
    ddtTan(u0, v, dvdt0, inhomo);
    cblas_daxpy(N_GRID, dt * RK[0][0], dvdt0, 1, v, 1);

    ddtTan(u1, v, dvdt1, inhomo);
    cblas_daxpy(N_GRID, dt * RK[1][0], dvdt0, 1, v, 1); 
    cblas_daxpy(N_GRID, dt * RK[1][1], dvdt1, 1, v, 1);

    ddtTan(u2, v, dvdt2, inhomo);
    cblas_daxpy(N_GRID, dt * RK[2][1], dvdt1, 1, v, 1);
    cblas_daxpy(N_GRID, dt * RK[2][2], dvdt2, 1, v, 1);
}


void
stepAdjoint(const double * u0, const double * v0, double strength,
            const double * w0, double * w, double dt)
{
    double u1[N_GRID], u2[N_GRID], dudt0[N_GRID], dudt1[N_GRID];
    double dwdt_u2w3[N_GRID], dwdt_u1w2[N_GRID], dwdt_u1w3[N_GRID],
           dwdt_u0w1[N_GRID], dwdt_u0w2[N_GRID];

    // Primal
    ddt(u0, dudt0);
    memmove(u1, u0, sizeof(double) * N_GRID);
    cblas_daxpy(N_GRID, dt * RK[0][0], dudt0, 1, u1, 1);

    ddt(u1, dudt1);
    memmove(u2, u1, sizeof(double) * N_GRID);
    cblas_daxpy(N_GRID, dt * RK[1][0], dudt0, 1, u2, 1);
    cblas_daxpy(N_GRID, dt * RK[1][1], dudt1, 1, u2, 1);

    // Adjoint -- w is w3
    memmove(w, w0, sizeof(double) * N_GRID);
    ddtAdj(u2, w, dwdt_u2w3);
    ddtAdj(u1, w, dwdt_u1w3);
    cblas_daxpy(N_GRID, -dt * RK[2][2], dwdt_u2w3, 1, w, 1);

    // w is now w2
    ddtAdj(u1, w, dwdt_u1w2);
    ddtAdj(u0, w, dwdt_u0w2);
    cblas_daxpy(N_GRID, -dt * RK[2][1], dwdt_u1w3, 1, w, 1); 
    cblas_daxpy(N_GRID, -dt * RK[1][1], dwdt_u1w2, 1, w, 1);

    // w is now w1
    ddtAdj(u0, w, dwdt_u0w1);
    cblas_daxpy(N_GRID, -dt * RK[1][0], dwdt_u0w2, 1, w, 1);
    cblas_daxpy(N_GRID, -dt * RK[0][0], dwdt_u0w1, 1, w, 1);

    // Adjoint -- source term
    double v0DotUt = cblas_ddot(N_GRID, v0, 1, dudt0, 1);
    double utDotUt = cblas_ddot(N_GRID, dudt0, 1, dudt0, 1);

    cblas_daxpy(N_GRID, dt * strength, v0, 1, w, 1);
    cblas_daxpy(N_GRID, dt * strength * (-v0DotUt / utDotUt), dudt0, 1, w, 1);
}


double
project_ddt(int i_chunk, double * v)
{
    assert (i_chunk >= 0 && i_chunk <= N_CHUNK);
    double dudt[N_GRID];
    ddt(SOLN_U[i_chunk][0], dudt);

    double vDotUt = cblas_ddot(N_GRID, v, 1, dudt, 1);
    double utDotUt = cblas_ddot(N_GRID, dudt, 1, dudt, 1);

    cblas_daxpy(N_GRID, -vDotUt / utDotUt, dudt, 1, v, 1);
    return vDotUt / utDotUt;
}


void
alloc_space_for_big_arrays_U_and_V()
{
    // allocate space for big arrays
    SOLN_U = (double ***) malloc(sizeof(double **) * (N_CHUNK + 1));
    SOLN_V = (double ***) malloc(sizeof(double **) * N_CHUNK);
    assert (SOLN_U != 0 && SOLN_V != 0);

    SOLN_U[0] = (double **) malloc(sizeof(double*) * (N_CHUNK*N_STEP + 1));
    SOLN_V[0] = (double **) malloc(sizeof(double*) * N_CHUNK*N_STEP);
    assert (SOLN_U[0] != 0 && SOLN_V[0] != 0);

    SOLN_U[0][0] = (double *) malloc(sizeof(double) * 
                                            (N_GRID + N_CHUNK*N_STEP*N_GRID));
    SOLN_V[0][0] = (double *) malloc(sizeof(double) * N_CHUNK*N_STEP*N_GRID);
    assert (SOLN_U[0][0] != 0 && SOLN_V[0][0] != 0);

    for (int i = 0; i < N_CHUNK; ++ i)
    {
        SOLN_U[i] = SOLN_U[0] + i * N_STEP;
        SOLN_V[i] = SOLN_V[0] + i * N_STEP;

        for (int j = 0; j < N_STEP; ++ j)
        {
            SOLN_U[i][j] = SOLN_U[0][0] + (i * N_STEP + j) * N_GRID;
            SOLN_V[i][j] = SOLN_V[0][0] + (i * N_STEP + j) * N_GRID;
        }
    }

    // This stores the last dangling solution
    SOLN_U[N_CHUNK] = SOLN_U[N_CHUNK - 1] + N_STEP;
    SOLN_U[N_CHUNK][0] = SOLN_U[N_CHUNK - 1][N_STEP - 1] + N_GRID;
}


void
run_to_T0_with_random_initial_condition(double * u, double T0, double dt_max)
{
    assert(T0 > 0);
    assert(dt_max > 0);
    int n0_steps = (int) ceil(T0 / dt_max);
    double dt0 = T0 / n0_steps;

    assert(u != 0);
    for (int i = 0; i < N_GRID; ++ i)
    {
        u[i] = rand() / (double) RAND_MAX;
    }

    for (int i = 0; i < n0_steps; ++ i)
    {
        stepPrimal(u, u, dt0);
    }
}


void
init(double c, int n_grid, double T0,
     int n_chunk, double t_chunk, double dt_max)
{
    C_CONST = c;

    assert(n_grid > 0);
    N_GRID = n_grid;

    assert(n_chunk > 0);
    N_CHUNK = n_chunk;
    assert(t_chunk > 0);
    N_STEP = (int) ceil(t_chunk / dt_max);
    DT_STEP = t_chunk / N_STEP;

    alloc_space_for_big_arrays_U_and_V();

    run_to_T0_with_random_initial_condition(SOLN_U[0][0], T0, dt_max);

    // run in each time chunk
    for (int i_chunk = 0; i_chunk < N_CHUNK; ++ i_chunk)
    {
        double ** u = SOLN_U[i_chunk];
        for (int i = 0; i < N_STEP; ++ i)
        {
            // when i == N_STEP - 1, u[i+1] is SOLN_U{i_chunk + 1][0]
            // there's one more memory slot at the last i_chunk, so don't worry
            stepPrimal(u[i], u[i+1], DT_STEP);
        }
    }
}


void
tangent(int i_chunk, double * v0, int inhomo)
{
    assert (i_chunk >= 0 && i_chunk < N_CHUNK);
    double ** u = SOLN_U[i_chunk];
    double ** v = SOLN_V[i_chunk];

    memmove(v[0], v0, sizeof(double) * N_GRID);

    for (int i = 0; i < N_STEP - 1; ++ i)
    {
        stepTangent(u[i], v[i], v[i+1], DT_STEP, inhomo);
    }
    stepTangent(u[N_STEP - 1], v[N_STEP - 1], v0, DT_STEP, inhomo);
}


void
adjoint(int i_chunk, double * w0, double forcing)
{
    assert (i_chunk >= 0 && i_chunk < N_CHUNK);
    double ** u = SOLN_U[i_chunk];
    double ** v = SOLN_V[i_chunk];

    for (int i = N_STEP - 1; i >= 0; -- i)
    {
        stepAdjoint(u[i], v[i], 1, w0, w0, DT_STEP);
    }
}

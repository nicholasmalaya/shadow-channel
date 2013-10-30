#include<string.h>
#include<stdio.h>
#include<stdlib.h>
#include<cblas.h>
#include"kuramoto.h"

int
main()
{
    const int n_grid = 127;
    double u0[n_grid];
    for (int i = 0; i < n_grid; ++ i) u0[i] = rand() / (double)RAND_MAX;
    init(0.5, u0, n_grid, 100, 10, 2, 0.2);


    const int n = 10;
    double u[n][N_GRID];
    double v0[N_GRID], w0[N_GRID];
    double v1[N_GRID], w1[N_GRID];

    for (int i = 0; i < N_GRID; ++ i)
    {
        u[0][i] = rand() / (double)RAND_MAX;
        v0[i] = rand() / (double)RAND_MAX;
        w1[i] = rand() / (double)RAND_MAX;
    }

    for (int i = 0; i < n-1; ++ i)
    {
        stepPrimal(u[i], u[i+1], DT_STEP);
    }

    memmove(v1, v0, sizeof(double) * N_GRID);
    for (int i = 0; i < n-1; ++ i)
    {
        stepTangent(u[i], v1, v1, DT_STEP, 0);
    }

    memmove(w0, w1, sizeof(double) * N_GRID);
    for (int i = n - 2; i >= 0; -- i)
    {
        stepAdjoint(u[i], v0, 0, w0, w0, DT_STEP);
    }

    double vDotW0 = cblas_ddot(N_GRID, v0, 1, w0, 1);
    double vDotW1 = cblas_ddot(N_GRID, v1, 1, w1, 1);

    printf("%f %f %f\n", vDotW0, vDotW1, vDotW0 - vDotW1);

    const double EPS = 0.0000001;
    double up[N_GRID], um[N_GRID];
    memmove(up, u[0], sizeof(double) * N_GRID);
    memmove(um, u[0], sizeof(double) * N_GRID);
    cblas_daxpy(N_GRID, +0.5 * EPS, v0, 1, up, 1);
    cblas_daxpy(N_GRID, -0.5 * EPS, v0, 1, um, 1);

    for (int i = 0; i < n-1; ++ i)
    {
        stepPrimal(up, up, DT_STEP);
        stepPrimal(um, um, DT_STEP);
    }

    double du[N_GRID];
    for (int i = 0; i < N_GRID; ++ i)
    {
        du[i] = (up[i] - um[i]) / EPS;
    }

    double v1Nrm2 = cblas_dnrm2(N_GRID, v1, 1);
    double duNrm2 = cblas_dnrm2(N_GRID, du, 1);

    cblas_daxpy(N_GRID, -1, v1, 1, du, 1);
    printf("%f %f %f\n", v1Nrm2, duNrm2, cblas_dnrm2(N_GRID, du, 1));

	

    return 0;
}

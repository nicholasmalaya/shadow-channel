#ifndef KURAMOTO_H
#define KURAMOTO_H

#include<assert.h>

void
ddt(const double * u, double * dudt);

void
ddtTan(const double * u, const double * v, double * dvdt, int inhomo);

void
ddtAdj(const double * u, const double * w, double * dwdt);

void
stepPrimal(const double * u0, double * u, double dt);

void
stepTangent(const double * u0, const double * v0, double * v, double dt,
            int inhomo);

void
stepAdjoint(const double * u0, const double * v0, double strength,
            const double * w0, double * w, double dt);

double
project_ddt(int i_chunk, int i_step, double * v);

void
init(double c, double * u0, int n_grid, double T0,
     int n_chunk, double t_chunk, double dt_max);

void
tangent(int i_chunk, double * v0, int inhomo);

void
adjoint(int i_chunk, double * w0, double forcing);

extern int N_GRID, N_CHUNK, N_STEP;
extern double C_CONST;
extern double DT_STEP;
extern double *** SOLN_U;
extern double *** SOLN_V;

#endif

#ifndef CHANNEL_CHANNEL_H
#define CHANNEL_CHANNEL_H

#include"mcomplex.h"

int init(int _Nx, int _Ny, int _Nz, double _Lx, double _Lz, double _Re,
         double _flux, double _dt, int _nsteps);

void primal(int ru_steps, mcomplex *C_given);

void tangent(int start_step, int end_step, mcomplex *IC_given, int inhomo);

void adjoint(int start_step, int end_step, mcomplex *AC_given, int inhomo, double strength);

void tangent_forcing0(int n, int k, int flag, mcomplex * f_a, mcomplex * f_b);

void tangent_forcing(int n, int k, int z, mcomplex ** f_a, mcomplex ** f_b);

void destroy(int status);

void getsoln(int i_step, mcomplex ** MC_ptr,
             int * Nz_ptr, int * Nvar_ptr, int * Ny_ptr, int * Nx_ptr);

void statistics(mcomplex * C_ptr, double * us_ptr);

void save_solution(char * filename, mcomplex * C_ptr);

void read_solution(char * filename, mcomplex * C_ptr);

#endif

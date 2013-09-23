#ifndef CHANNEL_CHANNEL_H
#define CHANNEL_CHANNEL_H

#include"mcomplex.h"

int init(int _Nx, int _Ny, int _Nz, double _Lx, double _Lz, double _Re,
         double _mpg, double _dt, int _rut,
         int _n_chunk, int _nsteps_chunk, int _restart_flag);

void destroy(int status);

void getsoln(int i_step, mcomplex ** MC_ptr,
             int * Nz_ptr, int * Nvar_ptr, int * Ny_ptr, int * Nx_ptr);

#endif

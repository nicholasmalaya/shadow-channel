%module channel
%{
#define SWIG_FILE_WITH_INIT
#include "channel.h"
%}


%include "numpy.i"

%init %{
import_array();
%}

int init(int _Nx, int _Ny, int _Nz, double _Lx, double _Lz, double _Re,
         double _mpg, double _dt, int _rut,
         int _n_chunk, int _nsteps_chunk, int _restart_flag);

void destroy(int status);

%numpy_typemaps(mcomplex, NPY_CDOUBLE, int);

%apply (mcomplex ** ARGOUTVIEW_ARRAY4, int * DIM1, int * DIM2, int * DIM3, int * DIM4)
      {(mcomplex ** MC_ptr, int * Nz_ptr, int * Nvar_ptr, int * Ny_ptr, int * Nx_ptr)}


void getsoln(int i_step, mcomplex ** MC_ptr,
             int * Nz_ptr, int * Nvar_ptr, int * Ny_ptr, int * Nx_ptr);

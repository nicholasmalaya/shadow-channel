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
         double _flux, double _dt, int _nstep);

void destroy(int status);

%numpy_typemaps(mcomplex, NPY_CDOUBLE, int);

%apply (mcomplex ** ARGOUTVIEW_ARRAY4, int * DIM1, int * DIM2, int * DIM3, int * DIM4)
      {(mcomplex ** MC_ptr, int * Nz_ptr, int * Nvar_ptr, int * Ny_ptr, int * Nx_ptr)}


void getsoln(int i_step, mcomplex ** MC_ptr,
             int * Nz_ptr, int * Nvar_ptr, int * Ny_ptr, int * Nx_ptr);

%apply (mcomplex * IN_ARRAY4, int DIM1, int DIM2, int DIM3, int DIM4)
      {(mcomplex * IC_given, int Nz_dup, int Nvar_dup, int Ny_dup, int Nx_dup)}

%inline %{
    void
    c_tangent(int start_step, int end_step, mcomplex ****IC_given, int Nz_dup, int Nvar_dup, int Ny_dup, int Nx_dup, int inhomo)
    {
        tangent(start_step, end_step, IC_given, inhomo);
    }
%}


%apply (mcomplex * IN_ARRAY4, int DIM1, int DIM2, int DIM3, int DIM4)
      {(mcomplex * C, int Nz_dup, int Nvar_dup, int Ny_dup, int Nx_dup)}
%apply (double ** ARGOUTVIEW_ARRAY2, int * DIM1, int * DIM2)
      {(double ** us_ptr, int * nstats_ptr, int * qpts_ptr)}

void statistics(mcomplex * C,
             int Nz_dup, int Nvar_dup, int Ny_dup, int Nx_dup,
             double ** us_ptr, int * nstats_ptr, int * qpts_ptr);

void restart2(int restart_flag);

void write_data2(int n);



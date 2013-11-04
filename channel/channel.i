%module c_channel
%{
#define SWIG_FILE_WITH_INIT
#include "main.h"
#include "channel.h"
%}

%include "numpy.i"

%init %{
import_array();
%}

%numpy_typemaps(mcomplex, NPY_CDOUBLE, int);


%inline %{
int c_init(int _Nx, int _Ny, int _Nz, double _Lx, double _Lz, double _Re,
           double _flux, double _dt, int _nstep)
{
    return init(_Nx, _Ny, _Nz, _Lx, _Lz, _Re, _flux, _dt, _nstep);
}
%}


%inline %{
void c_destroy()
{
    destroy(0xffffffff);
}
%}


%apply (mcomplex * INPLACE_ARRAY4, int DIM1, int DIM2, int DIM3, int DIM4)
      {(mcomplex * C_given, int Nz_dup, int Nvar_dup, int Ny_dup, int Nx_dup)}

%inline %{
void c_primal(int ru_steps,
        mcomplex *C_given, int Nz_dup, int Nvar_dup, int Ny_dup, int Nx_dup)
{
    if (Nx_dup == 0 && Ny_dup == 0 && Nz_dup == 0 && Nvar_dup == 0)
        primal(ru_steps, 0);
    else
    {
        assert(Nx_dup == Nx / 2);
        assert(Ny_dup == dimR);
        assert(Nz_dup == Nz);
        assert(Nvar_dup == 2);

        primal(ru_steps, C_given);
    }
}
%}


%apply (mcomplex * INPLACE_ARRAY4, int DIM1, int DIM2, int DIM3, int DIM4)
      {(mcomplex * IC_given, int Nz_dup, int Nvar_dup, int Ny_dup, int Nx_dup)}

%inline %{
void c_tangent(int start_step, int end_step,
        mcomplex *IC_given, int Nz_dup, int Nvar_dup, int Ny_dup, int Nx_dup,
        int inhomo)
{
    assert(Nx_dup == Nx / 2);
    assert(Ny_dup == dimR);
    assert(Nz_dup == Nz);
    assert(Nvar_dup == 2);

    tangent(start_step, end_step, IC_given, inhomo);
}
%}


%apply (mcomplex * INPLACE_ARRAY4, int DIM1, int DIM2, int DIM3, int DIM4)
      {(mcomplex * IC_given, int Nz_dup, int Nvar_dup, int Ny_dup, int Nx_dup)}

%inline %{
void c_adjoint(int start_step, int end_step,
        mcomplex *AC_given, int Nz_dup, int Nvar_dup, int Ny_dup, int Nx_dup,
        int inhomo)
{
    assert(Nx_dup == Nx / 2);
    assert(Ny_dup == dimR);
    assert(Nz_dup == Nz);
    assert(Nvar_dup == 2);

    adjoint(start_step, end_step, AC_given, inhomo);
}
%}


%apply (mcomplex ** ARGOUTVIEW_ARRAY4, int * DIM1, int * DIM2, int * DIM3, int * DIM4)
      {(mcomplex ** MC_ptr, int * Nz_ptr, int * Nvar_ptr, int * Ny_ptr, int * Nx_ptr)}

/* This is supposed to be used with ARGOUTVIEW_ARRAY4
   Look for the c_v function in kuramoto.i */
%inline %{
void c_getsoln(int i_step, mcomplex ** MC_ptr,
               int * Nz_ptr, int * Nvar_ptr, int * Ny_ptr, int * Nx_ptr)
{       
    assert(i_step >= 0 && i_step <= nsteps);
    (*MC_ptr) = MC[i_step * 3][0][0][0];
    (*Nz_ptr) = Nz;
    (*Nvar_ptr) = 2;
    (*Ny_ptr) = dimR;
    (*Nx_ptr) = Nx / 2;
}
%}


%apply (mcomplex * INPLACE_ARRAY4, int DIM1, int DIM2, int DIM3, int DIM4)
      {(mcomplex * C_ptr, int Nz_dup, int Nvar_dup, int Ny_dup, int Nx_dup)}
%apply (double * INPLACE_ARRAY2, int DIM1, int DIM2)
      {(double * us_ptr, int nstats_ptr, int qpts_ptr)}

%inline %{
void c_statistics(mcomplex * C_ptr,
             int Nz_dup, int Nvar_dup, int Ny_dup, int Nx_dup,
             double * us_ptr, int nstats_dup, int qpts_dup)
{
    assert(Nx_dup == Nx / 2);
    assert(Ny_dup == dimR);
    assert(Nz_dup == Nz);
    assert(Nvar_dup == 2);

    assert (nstats_dup == 21);
    assert (qpts_dup == qpts);

    statistics(C_ptr, us_ptr);
}
%}


%apply (mcomplex * INPLACE_ARRAY4, int DIM1, int DIM2, int DIM3, int DIM4)
      {(mcomplex * C_ptr, int Nz_dup, int Nvar_dup, int Ny_dup, int Nx_dup)}

%inline %{
void c_read_solution(char * filename,
        mcomplex * C_ptr, int Nz_dup, int Nvar_dup, int Ny_dup, int Nx_dup)
{
    assert(Nx_dup == Nx / 2);
    assert(Ny_dup == dimR);
    assert(Nz_dup == Nz);
    assert(Nvar_dup == 2);

    read_solution(filename, C_ptr);
}
%}


%apply (mcomplex * INPLACE_ARRAY4, int DIM1, int DIM2, int DIM3, int DIM4)
      {(mcomplex * C_ptr, int Nz_dup, int Nvar_dup, int Ny_dup, int Nx_dup)}

%inline %{
void c_save_solution(char * filename,
        mcomplex * C_ptr, int Nz_dup, int Nvar_dup, int Ny_dup, int Nx_dup)
{
    assert(Nx_dup == Nx / 2);
    assert(Ny_dup == dimR);
    assert(Nz_dup == Nz);
    assert(Nvar_dup == 2);

    save_solution(filename, C_ptr);
}
%}


/* global variables */
%inline %{
int c_qpts() { return qpts; }
int c_dimR() { return dimR; }
int c_dimQ() { return dimQ; }
int c_Nx() { return Nx; }
int c_Nz() { return Nz; }
int c_nsteps() { return nsteps; } 
%}

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
      {(mcomplex * AC_given, int Nz_dup, int Nvar_dup, int Ny_dup, int Nx_dup)}

%inline %{
void c_adjoint(int start_step, int end_step,
        mcomplex *AC_given, int Nz_dup, int Nvar_dup, int Ny_dup, int Nx_dup,
        int inhomo, double strength)
{
    assert(Nx_dup == Nx / 2);
    assert(Ny_dup == dimR);
    assert(Nz_dup == Nz);
    assert(Nvar_dup == 2);

    adjoint(start_step, end_step, AC_given, inhomo, strength);
}
%}


%apply (mcomplex * INPLACE_ARRAY4, int DIM1, int DIM2, int DIM3, int DIM4)
      {(mcomplex * IC_given, int Nz_dup, int Nvar_dup, int Ny_dup, int Nx_dup)}

%inline %{
double c_ddt_project(int i_step,
        mcomplex *IC_given, int Nz_dup, int Nvar_dup, int Ny_dup, int Nx_dup)
{
    assert(Nx_dup == Nx / 2);
    assert(Ny_dup == dimR);
    assert(Nz_dup == Nz);
    assert(Nvar_dup == 2);

    return ddt_project(i_step, IC_given);
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

%apply (mcomplex ** ARGOUTVIEW_ARRAY4, int * DIM1, int * DIM2, int * DIM3, int * DIM4)
      {(mcomplex ** MIC_ptr, int * Nz_ptr, int * Nvar_ptr, int * Ny_ptr, int * Nx_ptr)}

/* This is supposed to be used with ARGOUTVIEW_ARRAY4
   Look for the c_v function in kuramoto.i */
%inline %{
void c_getincresoln(int i_step, mcomplex ** MIC_ptr,
               int * Nz_ptr, int * Nvar_ptr, int * Ny_ptr, int * Nx_ptr)
{       
    assert(i_step >= 0 && i_step <= nsteps);
    (*MIC_ptr) = MIC[i_step * 3][0][0][0];
    (*Nz_ptr) = Nz;
    (*Nvar_ptr) = 2;
    (*Ny_ptr) = dimR;
    (*Nx_ptr) = Nx / 2;
}
%}


%apply (mcomplex * INPLACE_ARRAY4, int DIM1, int DIM2, int DIM3, int DIM4)
      {(mcomplex * C_ptr, int Nz_dup, int Nvar_dup, int Ny_dup, int Nx_dup)}
%apply (double * INPLACE_ARRAY4, int DIM1, int DIM2, int DIM3, int DIM4)
      {(double * flow_ptr, int Nd_flow, int Nx_flow, int Ny_flow, int Nz_flow)}

%inline %{
void c_spec2phys(mcomplex * C_ptr,
                 int Nz_dup, int Nvar_dup, int Ny_dup, int Nx_dup,
                 double * flow_ptr,
                 int Nd_flow, int Nx_flow, int Ny_flow, int Nz_flow)
{
    assert(Nx_dup == Nx / 2);
    assert(Ny_dup == dimR);
    assert(Nz_dup == Nz);
    assert(Nvar_dup == 2);

    assert (Nd_flow == 3);
    assert (Nx_flow == 3 * Nx / 2);
    assert (Ny_flow == qpts);
    assert (Nz_flow == 3 * Nz / 2);

    spec2phys(C_ptr, flow_ptr);
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


%apply (double ** ARGOUTVIEW_ARRAY1, int * DIM1)
      {(double ** y_ptr, int * qpts_ptr1)}
%apply (double ** ARGOUTVIEW_ARRAY1, int * DIM1)
      {(double ** w_ptr, int * qpts_ptr2)}

%inline %{
void c_quad(double ** y_ptr, int * qpts_ptr1, double ** w_ptr, int * qpts_ptr2)
{       
    (*y_ptr) = Qy;
    (*w_ptr) = W;
    (*qpts_ptr1) = (*qpts_ptr2) = qpts;
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

%module channel
%{
#define SWIG_FILE_WITH_INIT
#include "channel.h"
%}


%include "numpy.i"

%init %{
import_array();
%}

%apply (double* INPLACE_ARRAY1, int DIM1) {(double* u0, int n_grid)}

%inline %{
    void
    c_init(double c, double * u0, int n_grid, double T0,
           int n_chunk, double t_chunk, double dt_max)
    {
        init(c, u0, n_grid, T0, n_chunk, t_chunk, dt_max);
    }
%}

%apply (double* INPLACE_ARRAY1, int DIM1) {(double* v, int n)}

%inline %{
    double
    c_project_ddt(int i_chunk, int i_step, double * v, int n)
    {
        assert (n == N_GRID);
        return project_ddt(i_chunk, i_step, v);
    }
%}

%apply (double* INPLACE_ARRAY1, int DIM1) {(double* v0, int n0)}

%inline %{
    void
    c_tangent(int i_chunk, double * v0, int n0, int inhomo)
    {
        assert (n0 == N_GRID);
        tangent(i_chunk, v0, inhomo);
    }
%}

%apply (double* INPLACE_ARRAY1, int DIM1) {(double* w0, int n0)}

%inline %{
    void
    c_adjoint(int i_chunk, double * w0, int n0, double forcing)
    {
        assert (n0 == N_GRID);
        adjoint(i_chunk, w0, forcing);
    }
%}


%apply (double ** ARGOUTVIEW_ARRAY1, int * DIM1)
      {(double ** u_ptr, int * n_ptr)}

%inline %{
    void
    c_u(int i_chunk, int i_step, double ** u_ptr, int * n_ptr)
    {
        assert (i_chunk >= 0 && i_chunk < N_CHUNK);
        assert (i_step >= 0 && i_step <= N_STEP);
        *n_ptr = N_GRID;
        *u_ptr = SOLN_U[i_chunk][i_step];
    }
%}


%apply (double ** ARGOUTVIEW_ARRAY1, int * DIM1)
      {(double ** v_ptr, int * n_ptr)}

%inline %{
    void
    c_v(int i_chunk, int i_step, double ** v_ptr, int * n_ptr)
    {
        assert (i_chunk >= 0 && i_chunk < N_CHUNK);
        assert (i_step >= 0 && i_step < N_STEP);
        *n_ptr = N_GRID;
        *v_ptr = SOLN_V[i_chunk][i_step];
    }
%}


%immutable;
int N_GRID, N_CHUNK, N_STEP;
double C_CONST;
double DT_STEP;
%mutable;

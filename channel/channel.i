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
    c_channel_destroy()
    {
      channel_destroy();
    }
%}


%immutable;
double N_GRID;;
%mutable;

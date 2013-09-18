/* Copyright Qiqi Wang 2013 qiqi@mit.edu
   
   This is the first attempt to create an object oriented interface to
   the UT Austin channel DNS code for the Least Squares Shadowing method

   Original channel code by Vanessa, advised by Bob Moser
   Tangent and adjoint code by Shan, advised by Bob Moser
   This interface is written by Qiqi during 2013 summer J.T.Oden fellowship
*/

#ifndef CHANNEL_H
#define CHANNEL_H

int channel_construct(int _Nx, int _Ny, int _Nz, double _Lx, double _Lz,
                      double _dt, double _mpg, double _re, int _restart_flag);

int channel_size();

int channel_forward(double * uTan, double T0, double T1);

int channel_backward(double * uAdj, double T1, double T0);

int channel_ddt(double * dudt);

void channel_destroy();

extern double N_GRID;

#endif


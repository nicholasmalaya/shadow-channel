#ifndef MYFORCE_H
#define MYFORCE_H

#include <math.h>
#include <stdio.h>
#include <string.h>
#include "fftw.h"
#include "rfftw.h"
#include "minChnl.h"
#include "arrays.h"

void  force0(int n, int k, int flag, double *f, double *g);
void  increforce0(int n, int k, int flag, double *f, double *g);
void increforce(int n, int k, int z, double *f, double *g);
void  force0_2(int n, int k, int flag, double *f, double *g);
void  increforce0_2(int n, int k, int flag, double *f, double *g);
void increforce_2(int n, int k, int z, double *f, double *g);
void  force0_3(int n, int k, int flag, double *f, double *g);
void  force_3(int n, int k, int flag, double *f, double *g);
void  increforce0_3(int n, int k, int flag, double *f, double *g);
void increforce_3(int n, int k, int z, double *tmp, double *tmp2);
void  force0_4(int n, int k, int z, double *f, double *g);
void increforce_4(int n, int k, int z, double *tmp, double *tmp2);
void increforce_5(int n, int k, int z, double *tmp);
void  force0_5(int n, int k, int z, double *f);
#endif

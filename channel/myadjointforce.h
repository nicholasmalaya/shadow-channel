#ifndef MYADJOINTFORCE_H
#define MYADJOINTFORCE_H

#include <math.h>
#include <stdio.h>
#include <string.h>
#include "fftw.h"
#include "rfftw.h"
#include "minChnl.h"
#include "arrays.h"

void adjforce0(int n, int k, int flag, mcomplex * f, mcomplex * g);
void increadjforce0(int n, int k, int flag, mcomplex * f, mcomplex * g);
void increadjforce(int n, int k, int z, mcomplex * f, mcomplex * g);

void adjforce_2(int n, int k, int z, mcomplex * f, mcomplex * g);
void increadjforce0_2(int n, int k, int flag, mcomplex * f, mcomplex * g);
void increadjforce_2(int n, int k, int z, mcomplex * f, mcomplex * g);

void adjforce0_3(int n, int k, int flag, mcomplex * f, mcomplex * g);
void adjforce_3(int n, int k, int z, mcomplex * f, mcomplex * g);
void increadjforce0_3(int n, int k, int flag, mcomplex * f, mcomplex * g);
void increadjforce_3(int n, int k, int z, mcomplex * f, mcomplex * g);


void adjforce0_4(int n, int k, int flag, mcomplex * f, mcomplex * g);
void adjforce_4(int n, int k, int z, mcomplex * f, mcomplex * g);
void increadjforce0_4(int n, int k, int flag, mcomplex * f, mcomplex * g);
void increadjforce_4(int n, int k, int z, mcomplex * f, mcomplex * g);

void adjforce0_5(int n, int k, int flag, mcomplex * f, mcomplex * g);
void increadjforce0_5(int n, int k, int flag, mcomplex * f, mcomplex * g);

void adjforce0_6(int n, int k, int flag, mcomplex * f, mcomplex * g);
void increadjforce_6(int n, int k, int z, mcomplex * f, mcomplex * g);

void adjforce0_7(int n, int k, int flag, mcomplex * f, mcomplex * g);
void increadjforce_7(int n, int k, int z, mcomplex * f, mcomplex * g);
void force0_7(int n, int k, int flag, mcomplex * f, mcomplex * g);
void increforce_7(int n, int k, int z, mcomplex * f, mcomplex * g);


void force0_8(int n, int k, int flag, mcomplex * f, mcomplex * g);
void force_8(int n, int k, int z, mcomplex * f, mcomplex * g);
void increforce0_8(int n, int k, int flag, mcomplex * f, mcomplex * g);
void increforce_8(int n, int k, int z, mcomplex * f, mcomplex * g);
void adjforce_8(int n, int k, int flag, mcomplex * f, mcomplex * g);
void increadjforce0_8(int n, int k, int flag, mcomplex * f, mcomplex * g);
void increadjforce_8(int n, int k, int z, mcomplex * f, mcomplex * g);

void force0_10(int n, int k, int flag, mcomplex * f, mcomplex * g);
void increforce0_10(int n, int k, int flag, mcomplex * f, mcomplex * g);
void adjforce0_10(int n, int k, int flag, mcomplex * f, mcomplex * g);
void increadjforce0_10(int n, int k, int flag, mcomplex * f, mcomplex * g);
#endif

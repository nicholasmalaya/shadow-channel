/**********************************************************************
Author: Vanessa Lopez
	Department of Computer Science
	University of Illinois at Urbana-Champaign
	1304 W. Springfield Ave.
	Urbana, IL 61801-2987
	Email: vlopez@cse.uiuc.edu

Copyright 1999.  This code represents preliminary work toward the
author's thesis.  The code is not to be redistributed without the
author's permission.  Any corrections should be communicated to
the author.  Any modifications or reuse of the code must retain
acknowledgement of the original source.
**********************************************************************/

#ifndef MCOMPLEX_H
#define MCOMPLEX_H

#include "fftw.h"
typedef fftw_complex mcomplex;
#define Re(c)  ((c).re)
#define Im(c)  ((c).im)

#endif

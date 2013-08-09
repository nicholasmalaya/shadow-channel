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

#include <stdlib.h>
#include "fftw.h"
#include "rfftw.h"
#include "minChnl.h"
#include "arrays.h"

int getMem(void)
{
    /* External Variables */
    extern int qpts, dimR, Nx, Nz;
    extern double *cfl2, **MZ;
    extern double ***M;
    extern mcomplex ****U, ****C, ****IU, ****IC;
    extern mcomplex ****AU, ****AC, ****IAU, ****IAC;
    extern mcomplex *****MC, *****MIC;
    extern mcomplex **Fa, **Fb, **TM;
    extern mcomplex **IFa, **IFb, **ITM;
    extern mcomplex **Uxb, **Uzb;
    extern mcomplex **Uxbt, **Uzbt;
    extern mcomplex **AUxb, **AUzb;
    extern mcomplex **IUxb, **IUzb;
    extern mcomplex **IAUxb, **IAUzb;
    extern mcomplex **grad, **hess;
    extern mcomplex *fa, *fb, *tm;
    extern mcomplex *Ifa, *Ifb, *Itm;
    extern mcomplex **GUxb, **GUzb;
    extern mcomplex **GIUxb, **GIUzb;
    extern fftw_complex ***CT, ***ICT;
    extern mcomplex ****MU, ****MIU;
    extern mcomplex ****LU, ****LIU;
    extern mcomplex **HUxb, **HUzb;
    extern mcomplex **HAUxb, **HAUzb;
    if ((U = c4Darray(Nz, 5, qpts, Nx / 2)) == NULL) {
        printf("ERROR: No memory for U array.\n");
        return (NO_MEM);
    }
    if ((C = c4Darray(Nz, 2, dimR, Nx / 2)) == NULL) {
        printf("ERROR: No memory for C array.\n");
        freec4Darray(U);
        return (NO_MEM);
    }

    if ((CT = c3Darray(9, 3 * Nz / 2, 3 * Nx / 4 + 1)) == NULL) {
        printf("ERROR: No memory for CT array.\n");
        freec4Darray(U);
        freec4Darray(C);
        return (NO_MEM);
    }
    if ((Fa = cMatrix(dimR, Nx / 2)) == NULL) {
        printf("ERROR: No memory for Fa array.\n");
        freec4Darray(U);
        freec4Darray(C);
        freec3Darray(CT);
        return (NO_MEM);
    }
    if ((Fb = cMatrix(dimR, Nx / 2)) == NULL) {
        printf("ERROR: No memory for Fb array.\n");
        freec4Darray(U);
        freec4Darray(C);
        freec3Darray(CT);
        freecMatrix(Fa);
        return (NO_MEM);
    }
    if ((M = d3Darray(dimR, 9, Nx / 2)) == NULL) {
        printf("ERROR: No memory for M array.\n");
        freec4Darray(U);
        freec4Darray(C);
        freec3Darray(CT);
        freecMatrix(Fa);
        freecMatrix(Fb);
        return (NO_MEM);
    }
    if ((TM = cMatrix(dimR, Nx / 2)) == NULL) {
        printf("ERROR: No memory for TM array.\n");
        freec4Darray(U);
        freec4Darray(C);
        freec3Darray(CT);
        freecMatrix(Fa);
        freecMatrix(Fb);
        freed3Darray(M);
        return (NO_MEM);
    }
    if ((MZ = dMatrix(dimR, 5)) == NULL) {
        printf("ERROR: No memory for M0 array.\n");
        freec4Darray(U);
        freec4Darray(C);
        freec3Darray(CT);
        freecMatrix(Fa);
        freecMatrix(Fb);
        freed3Darray(M);
        freecMatrix(TM);
        return (NO_MEM);
    }
    if ((fa = cVector(dimR)) == NULL) {
        printf("ERROR: No memory for fa array.\n");
        freec4Darray(U);
        freec4Darray(C);
        freec3Darray(CT);
        freecMatrix(Fa);
        freecMatrix(Fb);
        freed3Darray(M);
        freecMatrix(TM);
        freedMatrix(MZ);
        return (NO_MEM);
    }
    if ((fb = cVector(dimR)) == NULL) {
        printf("ERROR: No memory for fb array.\n");
        freec4Darray(U);
        freec4Darray(C);
        freec3Darray(CT);
        freecMatrix(Fa);
        freecMatrix(Fb);
        freed3Darray(M);
        freecMatrix(TM);
        freedMatrix(MZ);
        freecVector(fa);
        return (NO_MEM);
    }
    if ((tm = cVector(dimR)) == NULL) {
        printf("ERROR: No memory for fb array.\n");
        freec4Darray(U);
        freec4Darray(C);
        freec3Darray(CT);
        freecMatrix(Fa);
        freecMatrix(Fb);
        freed3Darray(M);
        freecMatrix(TM);
        freedMatrix(MZ);
        freecVector(fa);
        freecVector(fb);
        return (NO_MEM);
    }
    if ((cfl2 = dVector(qpts)) == NULL) {
        printf("ERROR: No memory for cfl2 array.\n");
        freec4Darray(U);
        freec4Darray(C);
        freec3Darray(CT);
        freecMatrix(Fa);
        freecMatrix(Fb);
        freed3Darray(M);
        freecMatrix(TM);
        freedMatrix(MZ);
        freecVector(fa);
        freecVector(fb);
        freecVector(tm);
        return (NO_MEM);
    }

    if ((IU = c4Darray(Nz, 5, qpts, Nx / 2)) == NULL) {
        printf("ERROR: No memory for IU array.\n");
        freec4Darray(U);
        freec4Darray(C);
        freec3Darray(CT);
        freecMatrix(Fa);
        freecMatrix(Fb);
        freed3Darray(M);
        freecMatrix(TM);
        freedMatrix(MZ);
        freecVector(fa);
        freecVector(fb);
        freecVector(tm);
        freedVector(cfl2);
        return (NO_MEM);
    }

    if ((IC = c4Darray(Nz, 2, dimR, Nx / 2)) == NULL) {
        printf("ERROR: No memory for IC array.\n");
        freec4Darray(U);
        freec4Darray(C);
        freec3Darray(CT);
        freecMatrix(Fa);
        freecMatrix(Fb);
        freed3Darray(M);
        freecMatrix(TM);
        freedMatrix(MZ);
        freecVector(fa);
        freecVector(fb);
        freecVector(tm);
        freedVector(cfl2);
        freec4Darray(IU);
        return (NO_MEM);
    }

    if ((IFa = cMatrix(dimR, Nx / 2)) == NULL) {
        printf("ERROR: No memory for IFa array.\n");
        freec4Darray(U);
        freec4Darray(C);
        freec3Darray(CT);
        freecMatrix(Fa);
        freecMatrix(Fb);
        freed3Darray(M);
        freecMatrix(TM);
        freedMatrix(MZ);
        freecVector(fa);
        freecVector(fb);
        freecVector(tm);
        freedVector(cfl2);
        freec4Darray(IU);
        freec4Darray(IC);
        return (NO_MEM);
    }


    if ((IFb = cMatrix(dimR, Nx / 2)) == NULL) {
        printf("ERROR: No memory for IFb array.\n");
        freec4Darray(U);
        freec4Darray(C);
        freec3Darray(CT);
        freecMatrix(Fa);
        freecMatrix(Fb);
        freed3Darray(M);
        freecMatrix(TM);
        freedMatrix(MZ);
        freecVector(fa);
        freecVector(fb);
        freecVector(tm);
        freedVector(cfl2);
        freec4Darray(IU);
        freec4Darray(IC);
        freecMatrix(IFa);
        return (NO_MEM);
    }

    if ((ITM = cMatrix(dimR, Nx / 2)) == NULL) {
        printf("ERROR: No memory for ITM array.\n");
        freec4Darray(U);
        freec4Darray(C);
        freec3Darray(CT);
        freecMatrix(Fa);
        freecMatrix(Fb);
        freed3Darray(M);
        freecMatrix(TM);
        freedMatrix(MZ);
        freecVector(fa);
        freecVector(fb);
        freecVector(tm);
        freedVector(cfl2);
        freec4Darray(IU);
        freec4Darray(IC);
        freecMatrix(IFa);
        freecMatrix(IFb);
        return (NO_MEM);
    }

    if ((Ifa = cVector(dimR)) == NULL) {
        printf("ERROR: No memory for Ifa array.\n");
        freec4Darray(U);
        freec4Darray(C);
        freec3Darray(CT);
        freecMatrix(Fa);
        freecMatrix(Fb);
        freed3Darray(M);
        freecMatrix(TM);
        freedMatrix(MZ);
        freecVector(fa);
        freecVector(fb);
        freecVector(tm);
        freedVector(cfl2);
        freec4Darray(IU);
        freec4Darray(IC);
        freecMatrix(IFa);
        freecMatrix(IFb);
        freecMatrix(ITM);
        return (NO_ERR);
    }

    if ((Ifb = cVector(dimR)) == NULL) {
        printf("ERROR: No memory for Ifb array.\n");
        freec4Darray(U);
        freec4Darray(C);
        freec3Darray(CT);
        freecMatrix(Fa);
        freecMatrix(Fb);
        freed3Darray(M);
        freecMatrix(TM);
        freedMatrix(MZ);
        freecVector(fa);
        freecVector(fb);
        freecVector(tm);
        freedVector(cfl2);
        freec4Darray(IU);
        freec4Darray(IC);
        freecMatrix(IFa);
        freecMatrix(IFb);
        freecMatrix(ITM);
        freecVector(Ifa);
        return (NO_ERR);
    }


    if ((Itm = cVector(dimR)) == NULL) {
        printf("ERROR: No memory for Itm array.\n");
        freec4Darray(U);
        freec4Darray(C);
        freec3Darray(CT);
        freecMatrix(Fa);
        freecMatrix(Fb);
        freed3Darray(M);
        freecMatrix(TM);
        freedMatrix(MZ);
        freecVector(fa);
        freecVector(fb);
        freecVector(tm);
        freedVector(cfl2);
        freec4Darray(IU);
        freec4Darray(IC);
        freecMatrix(IFa);
        freecMatrix(IFb);
        freecMatrix(ITM);
        freecVector(Ifa);
        freecVector(Ifb);
        return (NO_ERR);
    }


    if ((ICT = c3Darray(9, 3 * Nz / 2, 3 * Nx / 4 + 1)) == NULL) {
        printf("ERROR: No memory for ICT array.\n");
        freec4Darray(U);
        freec4Darray(C);
        freec3Darray(CT);
        freecMatrix(Fa);
        freecMatrix(Fb);
        freed3Darray(M);
        freecMatrix(TM);
        freedMatrix(MZ);
        freecVector(fa);
        freecVector(fb);
        freecVector(tm);
        freedVector(cfl2);
        freec4Darray(IU);
        freec4Darray(IC);
        freecMatrix(IFa);
        freecMatrix(IFb);
        freecMatrix(ITM);
        freecVector(Ifa);
        freecVector(Ifb);
        freecVector(Itm);
        return (NO_MEM);
    }

    if ((Uxbt = cMatrix(Nz, Nx / 2)) == NULL) {
        printf("ERROR: No memory for Uxbt array.\n");
        freec4Darray(U);
        freec4Darray(C);
        freec3Darray(CT);
        freecMatrix(Fa);
        freecMatrix(Fb);
        freed3Darray(M);
        freecMatrix(TM);
        freedMatrix(MZ);
        freecVector(fa);
        freecVector(fb);
        freecVector(tm);
        freedVector(cfl2);
        freec4Darray(IU);
        freec4Darray(IC);
        freecMatrix(IFa);
        freecMatrix(IFb);
        freecMatrix(ITM);
        freecVector(Ifa);
        freecVector(Ifb);
        freecVector(Itm);
        freec3Darray(ICT);
        return (NO_MEM);
    }

    if ((Uzbt = cMatrix(Nz, Nx / 2)) == NULL) {
        printf("ERROR: No memory for Uzbt array.\n");
        freec4Darray(U);
        freec4Darray(C);
        freec3Darray(CT);
        freecMatrix(Fa);
        freecMatrix(Fb);
        freed3Darray(M);
        freecMatrix(TM);
        freedMatrix(MZ);
        freecVector(fa);
        freecVector(fb);
        freecVector(tm);
        freedVector(cfl2);
        freec4Darray(IU);
        freec4Darray(IC);
        freecMatrix(IFa);
        freecMatrix(IFb);
        freecMatrix(ITM);
        freecVector(Ifa);
        freecVector(Ifb);
        freecVector(Itm);
        freec3Darray(ICT);
        freecMatrix(Uxbt);
        return (NO_MEM);
    }


    if ((Uxb = cMatrix(Nz, Nx / 2)) == NULL) {
        printf("ERROR: No memory for Uxb array.\n");
        freec4Darray(U);
        freec4Darray(C);
        freec3Darray(CT);
        freecMatrix(Fa);
        freecMatrix(Fb);
        freed3Darray(M);
        freecMatrix(TM);
        freedMatrix(MZ);
        freecVector(fa);
        freecVector(fb);
        freecVector(tm);
        freedVector(cfl2);
        freec4Darray(IU);
        freec4Darray(IC);
        freecMatrix(IFa);
        freecMatrix(IFb);
        freecMatrix(ITM);
        freecVector(Ifa);
        freecVector(Ifb);
        freecVector(Itm);
        freec3Darray(ICT);
        freecMatrix(Uxbt);
        freecMatrix(Uzbt);
        return (NO_MEM);
    }

    if ((Uzb = cMatrix(Nz, Nx / 2)) == NULL) {
        printf("ERROR: No memory for Uzb array.\n");
        freec4Darray(U);
        freec4Darray(C);
        freec3Darray(CT);
        freecMatrix(Fa);
        freecMatrix(Fb);
        freed3Darray(M);
        freecMatrix(TM);
        freedMatrix(MZ);
        freecVector(fa);
        freecVector(fb);
        freecVector(tm);
        freedVector(cfl2);
        freec4Darray(IU);
        freec4Darray(IC);
        freecMatrix(IFa);
        freecMatrix(IFb);
        freecMatrix(ITM);
        freecVector(Ifa);
        freecVector(Ifb);
        freecVector(Itm);
        freec3Darray(ICT);
        freecMatrix(Uxbt);
        freecMatrix(Uzbt);
        freecMatrix(Uxb);
        return (NO_MEM);
    }

    if ((AU = c4Darray(Nz, 5, qpts, Nx / 2)) == NULL) {
        printf("ERROR: No memory for AU array.\n");
        freec4Darray(U);
        freec4Darray(C);
        freec3Darray(CT);
        freecMatrix(Fa);
        freecMatrix(Fb);
        freed3Darray(M);
        freecMatrix(TM);
        freedMatrix(MZ);
        freecVector(fa);
        freecVector(fb);
        freecVector(tm);
        freedVector(cfl2);
        freec4Darray(IU);
        freec4Darray(IC);
        freecMatrix(IFa);
        freecMatrix(IFb);
        freecMatrix(ITM);
        freecVector(Ifa);
        freecVector(Ifb);
        freecVector(Itm);
        freec3Darray(ICT);
        freecMatrix(Uxbt);
        freecMatrix(Uzbt);
        freecMatrix(Uxb);
        freecMatrix(Uzb);
        return (NO_MEM);
    }

    if ((IAU = c4Darray(Nz, 5, qpts, Nx / 2)) == NULL) {
        printf("ERROR: No memory for IAU array.\n");
        freec4Darray(U);
        freec4Darray(C);
        freec3Darray(CT);
        freecMatrix(Fa);
        freecMatrix(Fb);
        freed3Darray(M);
        freecMatrix(TM);
        freedMatrix(MZ);
        freecVector(fa);
        freecVector(fb);
        freecVector(tm);
        freedVector(cfl2);
        freec4Darray(IU);
        freec4Darray(IC);
        freecMatrix(IFa);
        freecMatrix(IFb);
        freecMatrix(ITM);
        freecVector(Ifa);
        freecVector(Ifb);
        freecVector(Itm);
        freec3Darray(ICT);
        freecMatrix(Uxbt);
        freecMatrix(Uzbt);
        freecMatrix(Uxb);
        freecMatrix(Uzb);
        freec4Darray(AU);
        return (NO_MEM);
    }


    if ((AC = c4Darray(Nz, 2, dimR, Nx / 2)) == NULL) {
        printf("ERROR: No memory for AC array.\n");
        freec4Darray(U);
        freec4Darray(C);
        freec3Darray(CT);
        freecMatrix(Fa);
        freecMatrix(Fb);
        freed3Darray(M);
        freecMatrix(TM);
        freedMatrix(MZ);
        freecVector(fa);
        freecVector(fb);
        freecVector(tm);
        freedVector(cfl2);
        freec4Darray(IU);
        freec4Darray(IC);
        freecMatrix(IFa);
        freecMatrix(IFb);
        freecMatrix(ITM);
        freecVector(Ifa);
        freecVector(Ifb);
        freecVector(Itm);
        freec3Darray(ICT);
        freecMatrix(Uxbt);
        freecMatrix(Uzbt);
        freecMatrix(Uxb);
        freecMatrix(Uzb);
        freec4Darray(AU);
        freec4Darray(IAU);
        return (NO_MEM);
    }


    if ((IAC = c4Darray(Nz, 2, dimR, Nx / 2)) == NULL) {
        printf("ERROR: No memory for IAC array.\n");
        freec4Darray(U);
        freec4Darray(C);
        freec3Darray(CT);
        freecMatrix(Fa);
        freecMatrix(Fb);
        freed3Darray(M);
        freecMatrix(TM);
        freedMatrix(MZ);
        freecVector(fa);
        freecVector(fb);
        freecVector(tm);
        freedVector(cfl2);
        freec4Darray(IU);
        freec4Darray(IC);
        freecMatrix(IFa);
        freecMatrix(IFb);
        freecMatrix(ITM);
        freecVector(Ifa);
        freecVector(Ifb);
        freecVector(Itm);
        freec3Darray(ICT);
        freecMatrix(Uxbt);
        freecMatrix(Uzbt);
        freecMatrix(Uxb);
        freecMatrix(Uzb);
        freec4Darray(AU);
        freec4Darray(IAU);
        freec4Darray(AC);
        return (NO_MEM);
    }

    if ((MC = c5Darray(3 * MAXSTEP + 1, Nz, 2, dimR, Nx / 2)) == NULL) {
        printf("ERROR: No memory for MC array.\n");
        freec4Darray(U);
        freec4Darray(C);
        freec3Darray(CT);
        freecMatrix(Fa);
        freecMatrix(Fb);
        freed3Darray(M);
        freecMatrix(TM);
        freedMatrix(MZ);
        freecVector(fa);
        freecVector(fb);
        freecVector(tm);
        freedVector(cfl2);
        freec4Darray(IU);
        freec4Darray(IC);
        freecMatrix(IFa);
        freecMatrix(IFb);
        freecMatrix(ITM);
        freecVector(Ifa);
        freecVector(Ifb);
        freecVector(Itm);
        freec3Darray(ICT);
        freecMatrix(Uxbt);
        freecMatrix(Uzbt);
        freecMatrix(Uxb);
        freecMatrix(Uzb);
        freec4Darray(AU);
        freec4Darray(IAU);
        freec4Darray(AC);
        freec4Darray(IAC);
        return (NO_MEM);
    }

    if ((MIC = c5Darray(3 * MAXSTEP + 1, Nz, 2, dimR, Nx / 2)) == NULL) {
        printf("ERROR: No memory for MIC array.\n");
        freec4Darray(U);
        freec4Darray(C);
        freec3Darray(CT);
        freecMatrix(Fa);
        freecMatrix(Fb);
        freed3Darray(M);
        freecMatrix(TM);
        freedMatrix(MZ);
        freecVector(fa);
        freecVector(fb);
        freecVector(tm);
        freedVector(cfl2);
        freec4Darray(IU);
        freec4Darray(IC);
        freecMatrix(IFa);
        freecMatrix(IFb);
        freecMatrix(ITM);
        freecVector(Ifa);
        freecVector(Ifb);
        freecVector(Itm);
        freec3Darray(ICT);
        freecMatrix(Uxbt);
        freecMatrix(Uzbt);
        freecMatrix(Uxb);
        freecMatrix(Uzb);
        freec4Darray(AU);
        freec4Darray(IAU);
        freec4Darray(AC);
        freec4Darray(IAC);
        freec5Darray(MC);
        return (NO_MEM);
    }

    if ((MU = c4Darray(Nz, 5, qpts, Nx / 2)) == NULL) {
        printf("ERROR: No memory for MU array.\n");
        freec4Darray(U);
        freec4Darray(C);
        freec3Darray(CT);
        freecMatrix(Fa);
        freecMatrix(Fb);
        freed3Darray(M);
        freecMatrix(TM);
        freedMatrix(MZ);
        freecVector(fa);
        freecVector(fb);
        freecVector(tm);
        freedVector(cfl2);
        freec4Darray(IU);
        freec4Darray(IC);
        freecMatrix(IFa);
        freecMatrix(IFb);
        freecMatrix(ITM);
        freecVector(Ifa);
        freecVector(Ifb);
        freecVector(Itm);
        freec3Darray(ICT);
        freecMatrix(Uxbt);
        freecMatrix(Uzbt);
        freecMatrix(Uxb);
        freecMatrix(Uzb);
        freec4Darray(AU);
        freec4Darray(IAU);
        freec4Darray(AC);
        freec4Darray(IAC);
        freec5Darray(MC);
        freec5Darray(MIC);
        return (NO_MEM);
    }

    if ((MIU = c4Darray(Nz, 5, qpts, Nx / 2)) == NULL) {
        printf("ERROR: No memory for MIU array.\n");
        freec4Darray(U);
        freec4Darray(C);
        freec3Darray(CT);
        freecMatrix(Fa);
        freecMatrix(Fb);
        freed3Darray(M);
        freecMatrix(TM);
        freedMatrix(MZ);
        freecVector(fa);
        freecVector(fb);
        freecVector(tm);
        freedVector(cfl2);
        freec4Darray(IU);
        freec4Darray(IC);
        freecMatrix(IFa);
        freecMatrix(IFb);
        freecMatrix(ITM);
        freecVector(Ifa);
        freecVector(Ifb);
        freecVector(Itm);
        freec3Darray(ICT);
        freecMatrix(Uxbt);
        freecMatrix(Uzbt);
        freecMatrix(Uxb);
        freecMatrix(Uzb);
        freec4Darray(AU);
        freec4Darray(IAU);
        freec4Darray(AC);
        freec4Darray(IAC);
        freec5Darray(MC);
        freec5Darray(MIC);
        freec4Darray(MU);
        return (NO_MEM);
    }

    if ((LU = c4Darray(Nz, 5, qpts, Nx / 2)) == NULL) {
        printf("ERROR: No memory for LU array.\n");
        freec4Darray(U);
        freec4Darray(C);
        freec3Darray(CT);
        freecMatrix(Fa);
        freecMatrix(Fb);
        freed3Darray(M);
        freecMatrix(TM);
        freedMatrix(MZ);
        freecVector(fa);
        freecVector(fb);
        freecVector(tm);
        freedVector(cfl2);
        freec4Darray(IU);
        freec4Darray(IC);
        freecMatrix(IFa);
        freecMatrix(IFb);
        freecMatrix(ITM);
        freecVector(Ifa);
        freecVector(Ifb);
        freecVector(Itm);
        freec3Darray(ICT);
        freecMatrix(Uxbt);
        freecMatrix(Uzbt);
        freecMatrix(Uxb);
        freecMatrix(Uzb);
        freec4Darray(AU);
        freec4Darray(IAU);
        freec4Darray(AC);
        freec4Darray(IAC);
        freec5Darray(MC);
        freec5Darray(MIC);
        freec4Darray(MU);
        freec4Darray(MIU);
        return (NO_MEM);
    }

    if ((LIU = c4Darray(Nz, 5, qpts, Nx / 2)) == NULL) {
        printf("ERROR: No memory for LIU array.\n");
        freec4Darray(U);
        freec4Darray(C);
        freec3Darray(CT);
        freecMatrix(Fa);
        freecMatrix(Fb);
        freed3Darray(M);
        freecMatrix(TM);
        freedMatrix(MZ);
        freecVector(fa);
        freecVector(fb);
        freecVector(tm);
        freedVector(cfl2);
        freec4Darray(IU);
        freec4Darray(IC);
        freecMatrix(IFa);
        freecMatrix(IFb);
        freecMatrix(ITM);
        freecVector(Ifa);
        freecVector(Ifb);
        freecVector(Itm);
        freec3Darray(ICT);
        freecMatrix(Uxbt);
        freecMatrix(Uzbt);
        freecMatrix(Uxb);
        freecMatrix(Uzb);
        freec4Darray(AU);
        freec4Darray(IAU);
        freec4Darray(AC);
        freec4Darray(IAC);
        freec5Darray(MC);
        freec5Darray(MIC);
        freec4Darray(MU);
        freec4Darray(MIU);
        freec4Darray(LU);
        return (NO_MEM);
    }

    if ((AUxb = cMatrix(Nz, Nx / 2)) == NULL) {
        printf("ERROR: No memory for AUxb array.\n");
        freec4Darray(U);
        freec4Darray(C);
        freec3Darray(CT);
        freecMatrix(Fa);
        freecMatrix(Fb);
        freed3Darray(M);
        freecMatrix(TM);
        freedMatrix(MZ);
        freecVector(fa);
        freecVector(fb);
        freecVector(tm);
        freedVector(cfl2);
        freec4Darray(IU);
        freec4Darray(IC);
        freecMatrix(IFa);
        freecMatrix(IFb);
        freecMatrix(ITM);
        freecVector(Ifa);
        freecVector(Ifb);
        freecVector(Itm);
        freec3Darray(ICT);
        freecMatrix(Uxbt);
        freecMatrix(Uzbt);
        freecMatrix(Uxb);
        freecMatrix(Uzb);
        freec4Darray(AU);
        freec4Darray(IAU);
        freec4Darray(AC);
        freec4Darray(IAC);
        freec5Darray(MC);
        freec5Darray(MIC);
        freec4Darray(MU);
        freec4Darray(MIU);
        freec4Darray(LU);
        freec4Darray(LIU);
        return (NO_MEM);
    }

    if ((AUzb = cMatrix(Nz, Nx / 2)) == NULL) {
        printf("ERROR: No memory for AUzb array.\n");
        freec4Darray(U);
        freec4Darray(C);
        freec3Darray(CT);
        freecMatrix(Fa);
        freecMatrix(Fb);
        freed3Darray(M);
        freecMatrix(TM);
        freedMatrix(MZ);
        freecVector(fa);
        freecVector(fb);
        freecVector(tm);
        freedVector(cfl2);
        freec4Darray(IU);
        freec4Darray(IC);
        freecMatrix(IFa);
        freecMatrix(IFb);
        freecMatrix(ITM);
        freecVector(Ifa);
        freecVector(Ifb);
        freecVector(Itm);
        freec3Darray(ICT);
        freecMatrix(Uxbt);
        freecMatrix(Uzbt);
        freecMatrix(Uxb);
        freecMatrix(Uzb);
        freec4Darray(AU);
        freec4Darray(IAU);
        freec4Darray(AC);
        freec4Darray(IAC);
        freec5Darray(MC);
        freec5Darray(MIC);
        freec4Darray(MU);
        freec4Darray(MIU);
        freec4Darray(LU);
        freec4Darray(LIU);
        freecMatrix(AUxb);
        return (NO_MEM);
    }

    if ((IUxb = cMatrix(Nz, Nx / 2)) == NULL) {
        printf("ERROR: No memory for IUxb array.\n");
        freec4Darray(U);
        freec4Darray(C);
        freec3Darray(CT);
        freecMatrix(Fa);
        freecMatrix(Fb);
        freed3Darray(M);
        freecMatrix(TM);
        freedMatrix(MZ);
        freecVector(fa);
        freecVector(fb);
        freecVector(tm);
        freedVector(cfl2);
        freec4Darray(IU);
        freec4Darray(IC);
        freecMatrix(IFa);
        freecMatrix(IFb);
        freecMatrix(ITM);
        freecVector(Ifa);
        freecVector(Ifb);
        freecVector(Itm);
        freec3Darray(ICT);
        freecMatrix(Uxbt);
        freecMatrix(Uzbt);
        freecMatrix(Uxb);
        freecMatrix(Uzb);
        freec4Darray(AU);
        freec4Darray(IAU);
        freec4Darray(AC);
        freec4Darray(IAC);
        freec5Darray(MC);
        freec5Darray(MIC);
        freec4Darray(MU);
        freec4Darray(MIU);
        freec4Darray(LU);
        freec4Darray(LIU);
        freecMatrix(AUxb);
        freecMatrix(AUzb);
        return (NO_MEM);
    }

    if ((IUzb = cMatrix(Nz, Nx / 2)) == NULL) {
        printf("ERROR: No memory for IUzb array.\n");
        freec4Darray(U);
        freec4Darray(C);
        freec3Darray(CT);
        freecMatrix(Fa);
        freecMatrix(Fb);
        freed3Darray(M);
        freecMatrix(TM);
        freedMatrix(MZ);
        freecVector(fa);
        freecVector(fb);
        freecVector(tm);
        freedVector(cfl2);
        freec4Darray(IU);
        freec4Darray(IC);
        freecMatrix(IFa);
        freecMatrix(IFb);
        freecMatrix(ITM);
        freecVector(Ifa);
        freecVector(Ifb);
        freecVector(Itm);
        freec3Darray(ICT);
        freecMatrix(Uxbt);
        freecMatrix(Uzbt);
        freecMatrix(Uxb);
        freecMatrix(Uzb);
        freec4Darray(AU);
        freec4Darray(IAU);
        freec4Darray(AC);
        freec4Darray(IAC);
        freec5Darray(MC);
        freec5Darray(MIC);
        freec4Darray(MU);
        freec4Darray(MIU);
        freec4Darray(LU);
        freec4Darray(LIU);
        freecMatrix(AUxb);
        freecMatrix(AUzb);
        freecMatrix(IUxb);
        return (NO_MEM);
    }

    if ((IAUxb = cMatrix(Nz, Nx / 2)) == NULL) {
        printf("ERROR: No memory for IAUxb array.\n");
        freec4Darray(U);
        freec4Darray(C);
        freec3Darray(CT);
        freecMatrix(Fa);
        freecMatrix(Fb);
        freed3Darray(M);
        freecMatrix(TM);
        freedMatrix(MZ);
        freecVector(fa);
        freecVector(fb);
        freecVector(tm);
        freedVector(cfl2);
        freec4Darray(IU);
        freec4Darray(IC);
        freecMatrix(IFa);
        freecMatrix(IFb);
        freecMatrix(ITM);
        freecVector(Ifa);
        freecVector(Ifb);
        freecVector(Itm);
        freec3Darray(ICT);
        freecMatrix(Uxbt);
        freecMatrix(Uzbt);
        freecMatrix(Uxb);
        freecMatrix(Uzb);
        freec4Darray(AU);
        freec4Darray(IAU);
        freec4Darray(AC);
        freec4Darray(IAC);
        freec5Darray(MC);
        freec5Darray(MIC);
        freec4Darray(MU);
        freec4Darray(MIU);
        freec4Darray(LU);
        freec4Darray(LIU);
        freecMatrix(AUxb);
        freecMatrix(AUzb);
        freecMatrix(IUzb);
        freecMatrix(IUxb);
        return (NO_MEM);
    }

    if ((IAUzb = cMatrix(Nz, Nx / 2)) == NULL) {
        printf("ERROR: No memory for IUzb array.\n");
        freec4Darray(U);
        freec4Darray(C);
        freec3Darray(CT);
        freecMatrix(Fa);
        freecMatrix(Fb);
        freed3Darray(M);
        freecMatrix(TM);
        freedMatrix(MZ);
        freecVector(fa);
        freecVector(fb);
        freecVector(tm);
        freedVector(cfl2);
        freec4Darray(IU);
        freec4Darray(IC);
        freecMatrix(IFa);
        freecMatrix(IFb);
        freecMatrix(ITM);
        freecVector(Ifa);
        freecVector(Ifb);
        freecVector(Itm);
        freec3Darray(ICT);
        freecMatrix(Uxbt);
        freecMatrix(Uzbt);
        freecMatrix(Uxb);
        freecMatrix(Uzb);
        freec4Darray(AU);
        freec4Darray(IAU);
        freec4Darray(AC);
        freec4Darray(IAC);
        freec5Darray(MC);
        freec5Darray(MIC);
        freec4Darray(MU);
        freec4Darray(MIU);
        freec4Darray(LU);
        freec4Darray(LIU);
        freecMatrix(AUxb);
        freecMatrix(AUzb);
        freecMatrix(IUzb);
        freecMatrix(IUxb);
        freecMatrix(IAUxb);
        return (NO_MEM);
    }

    if ((grad = cMatrix(Nz, Nx / 2)) == NULL) {
        printf("ERROR: No memory for grad  array.\n");
        freec4Darray(U);
        freec4Darray(C);
        freec3Darray(CT);
        freecMatrix(Fa);
        freecMatrix(Fb);
        freed3Darray(M);
        freecMatrix(TM);
        freedMatrix(MZ);
        freecVector(fa);
        freecVector(fb);
        freecVector(tm);
        freedVector(cfl2);
        freec4Darray(IU);
        freec4Darray(IC);
        freecMatrix(IFa);
        freecMatrix(IFb);
        freecMatrix(ITM);
        freecVector(Ifa);
        freecVector(Ifb);
        freecVector(Itm);
        freec3Darray(ICT);
        freecMatrix(Uxbt);
        freecMatrix(Uzbt);
        freecMatrix(Uxb);
        freecMatrix(Uzb);
        freec4Darray(AU);
        freec4Darray(IAU);
        freec4Darray(AC);
        freec4Darray(IAC);
        freec5Darray(MC);
        freec5Darray(MIC);
        freec4Darray(MU);
        freec4Darray(MIU);
        freec4Darray(LU);
        freec4Darray(LIU);
        freecMatrix(AUxb);
        freecMatrix(AUzb);
        freecMatrix(IUzb);
        freecMatrix(IUxb);
        freecMatrix(IAUxb);
        freecMatrix(IAUzb);
        return (NO_MEM);
    }


    if ((GUxb = cMatrix(Nz, Nx / 2)) == NULL) {
        printf("ERROR: No memory for Guxb  array.\n");
        freec4Darray(U);
        freec4Darray(C);
        freec3Darray(CT);
        freecMatrix(Fa);
        freecMatrix(Fb);
        freed3Darray(M);
        freecMatrix(TM);
        freedMatrix(MZ);
        freecVector(fa);
        freecVector(fb);
        freecVector(tm);
        freedVector(cfl2);
        freec4Darray(IU);
        freec4Darray(IC);
        freecMatrix(IFa);
        freecMatrix(IFb);
        freecMatrix(ITM);
        freecVector(Ifa);
        freecVector(Ifb);
        freecVector(Itm);
        freec3Darray(ICT);
        freecMatrix(Uxbt);
        freecMatrix(Uzbt);
        freecMatrix(Uxb);
        freecMatrix(Uzb);
        freec4Darray(AU);
        freec4Darray(IAU);
        freec4Darray(AC);
        freec4Darray(IAC);
        freec5Darray(MC);
        freec5Darray(MIC);
        freec4Darray(MU);
        freec4Darray(MIU);
        freec4Darray(LU);
        freec4Darray(LIU);
        freecMatrix(AUxb);
        freecMatrix(AUzb);
        freecMatrix(IUzb);
        freecMatrix(IUxb);
        freecMatrix(IAUxb);
        freecMatrix(IAUzb);
        freecMatrix(grad);
        return (NO_MEM);
    }

    if ((GUzb = cMatrix(Nz, Nx / 2)) == NULL) {
        printf("ERROR: No memory for Guzb  array.\n");
        freec4Darray(U);
        freec4Darray(C);
        freec3Darray(CT);
        freecMatrix(Fa);
        freecMatrix(Fb);
        freed3Darray(M);
        freecMatrix(TM);
        freedMatrix(MZ);
        freecVector(fa);
        freecVector(fb);
        freecVector(tm);
        freedVector(cfl2);
        freec4Darray(IU);
        freec4Darray(IC);
        freecMatrix(IFa);
        freecMatrix(IFb);
        freecMatrix(ITM);
        freecVector(Ifa);
        freecVector(Ifb);
        freecVector(Itm);
        freec3Darray(ICT);
        freecMatrix(Uxbt);
        freecMatrix(Uzbt);
        freecMatrix(Uxb);
        freecMatrix(Uzb);
        freec4Darray(AU);
        freec4Darray(IAU);
        freec4Darray(AC);
        freec4Darray(IAC);
        freec5Darray(MC);
        freec5Darray(MIC);
        freec4Darray(MU);
        freec4Darray(MIU);
        freec4Darray(LU);
        freec4Darray(LIU);
        freecMatrix(AUxb);
        freecMatrix(AUzb);
        freecMatrix(IUzb);
        freecMatrix(IUxb);
        freecMatrix(IAUxb);
        freecMatrix(IAUzb);
        freecMatrix(grad);
        freecMatrix(GUxb);
        return (NO_MEM);
    }

    if ((GIUxb = cMatrix(Nz, Nx / 2)) == NULL) {
        printf("ERROR: No memory for GIuxb  array.\n");
        freec4Darray(U);
        freec4Darray(C);
        freec3Darray(CT);
        freecMatrix(Fa);
        freecMatrix(Fb);
        freed3Darray(M);
        freecMatrix(TM);
        freedMatrix(MZ);
        freecVector(fa);
        freecVector(fb);
        freecVector(tm);
        freedVector(cfl2);
        freec4Darray(IU);
        freec4Darray(IC);
        freecMatrix(IFa);
        freecMatrix(IFb);
        freecMatrix(ITM);
        freecVector(Ifa);
        freecVector(Ifb);
        freecVector(Itm);
        freec3Darray(ICT);
        freecMatrix(Uxbt);
        freecMatrix(Uzbt);
        freecMatrix(Uxb);
        freecMatrix(Uzb);
        freec4Darray(AU);
        freec4Darray(IAU);
        freec4Darray(AC);
        freec4Darray(IAC);
        freec5Darray(MC);
        freec5Darray(MIC);
        freec4Darray(MU);
        freec4Darray(MIU);
        freec4Darray(LU);
        freec4Darray(LIU);
        freecMatrix(AUxb);
        freecMatrix(AUzb);
        freecMatrix(IUzb);
        freecMatrix(IUxb);
        freecMatrix(IAUxb);
        freecMatrix(IAUzb);
        freecMatrix(grad);
        freecMatrix(GUxb);
        freecMatrix(GUzb);
        return (NO_MEM);
    }

    if ((GIUzb = cMatrix(Nz, Nx / 2)) == NULL) {
        printf("ERROR: No memory for GIUzb  array.\n");
        freec4Darray(U);
        freec4Darray(C);
        freec3Darray(CT);
        freecMatrix(Fa);
        freecMatrix(Fb);
        freed3Darray(M);
        freecMatrix(TM);
        freedMatrix(MZ);
        freecVector(fa);
        freecVector(fb);
        freecVector(tm);
        freedVector(cfl2);
        freec4Darray(IU);
        freec4Darray(IC);
        freecMatrix(IFa);
        freecMatrix(IFb);
        freecMatrix(ITM);
        freecVector(Ifa);
        freecVector(Ifb);
        freecVector(Itm);
        freec3Darray(ICT);
        freecMatrix(Uxbt);
        freecMatrix(Uzbt);
        freecMatrix(Uxb);
        freecMatrix(Uzb);
        freec4Darray(AU);
        freec4Darray(IAU);
        freec4Darray(AC);
        freec4Darray(IAC);
        freec5Darray(MC);
        freec5Darray(MIC);
        freec4Darray(MU);
        freec4Darray(MIU);
        freec4Darray(LU);
        freec4Darray(LIU);
        freecMatrix(AUxb);
        freecMatrix(AUzb);
        freecMatrix(IUzb);
        freecMatrix(IUxb);
        freecMatrix(IAUxb);
        freecMatrix(IAUzb);
        freecMatrix(grad);
        freecMatrix(GUxb);
        freecMatrix(GUzb);
        freecMatrix(GIUxb);
        return (NO_MEM);
    }

    if ((HUxb = cMatrix(Nz, Nx / 2)) == NULL) {
        printf("ERROR: No memory for HUxb  array.\n");
        freec4Darray(U);
        freec4Darray(C);
        freec3Darray(CT);
        freecMatrix(Fa);
        freecMatrix(Fb);
        freed3Darray(M);
        freecMatrix(TM);
        freedMatrix(MZ);
        freecVector(fa);
        freecVector(fb);
        freecVector(tm);
        freedVector(cfl2);
        freec4Darray(IU);
        freec4Darray(IC);
        freecMatrix(IFa);
        freecMatrix(IFb);
        freecMatrix(ITM);
        freecVector(Ifa);
        freecVector(Ifb);
        freecVector(Itm);
        freec3Darray(ICT);
        freecMatrix(Uxbt);
        freecMatrix(Uzbt);
        freecMatrix(Uxb);
        freecMatrix(Uzb);
        freec4Darray(AU);
        freec4Darray(IAU);
        freec4Darray(AC);
        freec4Darray(IAC);
        freec5Darray(MC);
        freec5Darray(MIC);
        freec4Darray(MU);
        freec4Darray(MIU);
        freec4Darray(LU);
        freec4Darray(LIU);
        freecMatrix(AUxb);
        freecMatrix(AUzb);
        freecMatrix(IUzb);
        freecMatrix(IUxb);
        freecMatrix(IAUxb);
        freecMatrix(IAUzb);
        freecMatrix(grad);
        freecMatrix(GUxb);
        freecMatrix(GUzb);
        freecMatrix(GIUxb);
        freecMatrix(GIUzb);
        return (NO_MEM);
    }

    if ((HUzb = cMatrix(Nz, Nx / 2)) == NULL) {
        printf("ERROR: No memory for HUzb  array.\n");
        freec4Darray(U);
        freec4Darray(C);
        freec3Darray(CT);
        freecMatrix(Fa);
        freecMatrix(Fb);
        freed3Darray(M);
        freecMatrix(TM);
        freedMatrix(MZ);
        freecVector(fa);
        freecVector(fb);
        freecVector(tm);
        freedVector(cfl2);
        freec4Darray(IU);
        freec4Darray(IC);
        freecMatrix(IFa);
        freecMatrix(IFb);
        freecMatrix(ITM);
        freecVector(Ifa);
        freecVector(Ifb);
        freecVector(Itm);
        freec3Darray(ICT);
        freecMatrix(Uxbt);
        freecMatrix(Uzbt);
        freecMatrix(Uxb);
        freecMatrix(Uzb);
        freec4Darray(AU);
        freec4Darray(IAU);
        freec4Darray(AC);
        freec4Darray(IAC);
        freec5Darray(MC);
        freec5Darray(MIC);
        freec4Darray(MU);
        freec4Darray(MIU);
        freec4Darray(LU);
        freec4Darray(LIU);
        freecMatrix(AUxb);
        freecMatrix(AUzb);
        freecMatrix(IUzb);
        freecMatrix(IUxb);
        freecMatrix(IAUxb);
        freecMatrix(IAUzb);
        freecMatrix(grad);
        freecMatrix(GUxb);
        freecMatrix(GUzb);
        freecMatrix(GIUxb);
        freecMatrix(GIUzb);
        freecMatrix(HUxb);
        return (NO_MEM);
    }
    if ((HAUxb = cMatrix(Nz, Nx / 2)) == NULL) {
        printf("ERROR: No memory for HIUxb  array.\n");
        freec4Darray(U);
        freec4Darray(C);
        freec3Darray(CT);
        freecMatrix(Fa);
        freecMatrix(Fb);
        freed3Darray(M);
        freecMatrix(TM);
        freedMatrix(MZ);
        freecVector(fa);
        freecVector(fb);
        freecVector(tm);
        freedVector(cfl2);
        freec4Darray(IU);
        freec4Darray(IC);
        freecMatrix(IFa);
        freecMatrix(IFb);
        freecMatrix(ITM);
        freecVector(Ifa);
        freecVector(Ifb);
        freecVector(Itm);
        freec3Darray(ICT);
        freecMatrix(Uxbt);
        freecMatrix(Uzbt);
        freecMatrix(Uxb);
        freecMatrix(Uzb);
        freec4Darray(AU);
        freec4Darray(IAU);
        freec4Darray(AC);
        freec4Darray(IAC);
        freec5Darray(MC);
        freec5Darray(MIC);
        freec4Darray(MU);
        freec4Darray(MIU);
        freec4Darray(LU);
        freec4Darray(LIU);
        freecMatrix(AUxb);
        freecMatrix(AUzb);
        freecMatrix(IUzb);
        freecMatrix(IUxb);
        freecMatrix(IAUxb);
        freecMatrix(IAUzb);
        freecMatrix(grad);
        freecMatrix(GUxb);
        freecMatrix(GUzb);
        freecMatrix(GIUxb);
        freecMatrix(GIUzb);
        freecMatrix(HUxb);
        freecMatrix(HUzb);
        return (NO_MEM);
    }

    if ((HAUzb = cMatrix(Nz, Nx / 2)) == NULL) {
        printf("ERROR: No memory for HAUzb  array.\n");
        freec4Darray(U);
        freec4Darray(C);
        freec3Darray(CT);
        freecMatrix(Fa);
        freecMatrix(Fb);
        freed3Darray(M);
        freecMatrix(TM);
        freedMatrix(MZ);
        freecVector(fa);
        freecVector(fb);
        freecVector(tm);
        freedVector(cfl2);
        freec4Darray(IU);
        freec4Darray(IC);
        freecMatrix(IFa);
        freecMatrix(IFb);
        freecMatrix(ITM);
        freecVector(Ifa);
        freecVector(Ifb);
        freecVector(Itm);
        freec3Darray(ICT);
        freecMatrix(Uxbt);
        freecMatrix(Uzbt);
        freecMatrix(Uxb);
        freecMatrix(Uzb);
        freec4Darray(AU);
        freec4Darray(IAU);
        freec4Darray(AC);
        freec4Darray(IAC);
        freec5Darray(MC);
        freec5Darray(MIC);
        freec4Darray(MU);
        freec4Darray(MIU);
        freec4Darray(LU);
        freec4Darray(LIU);
        freecMatrix(AUxb);
        freecMatrix(AUzb);
        freecMatrix(IUzb);
        freecMatrix(IUxb);
        freecMatrix(IAUxb);
        freecMatrix(IAUzb);
        freecMatrix(grad);
        freecMatrix(GUxb);
        freecMatrix(GUzb);
        freecMatrix(GIUxb);
        freecMatrix(GIUzb);
        freecMatrix(HUxb);
        freecMatrix(HUzb);
        freecMatrix(HAUxb);
        return (NO_MEM);
    }
    if ((hess = cMatrix(Nz, Nx / 2)) == NULL) {
        printf("ERROR: No memory for hess array.\n");
        freec4Darray(U);
        freec4Darray(C);
        freec3Darray(CT);
        freecMatrix(Fa);
        freecMatrix(Fb);
        freed3Darray(M);
        freecMatrix(TM);
        freedMatrix(MZ);
        freecVector(fa);
        freecVector(fb);
        freecVector(tm);
        freedVector(cfl2);
        freec4Darray(IU);
        freec4Darray(IC);
        freecMatrix(IFa);
        freecMatrix(IFb);
        freecMatrix(ITM);
        freecVector(Ifa);
        freecVector(Ifb);
        freecVector(Itm);
        freec3Darray(ICT);
        freecMatrix(Uxbt);
        freecMatrix(Uzbt);
        freecMatrix(Uxb);
        freecMatrix(Uzb);
        freec4Darray(AU);
        freec4Darray(IAU);
        freec4Darray(AC);
        freec4Darray(IAC);
        freec5Darray(MC);
        freec5Darray(MIC);
        freec4Darray(MU);
        freec4Darray(MIU);
        freec4Darray(LU);
        freec4Darray(LIU);
        freecMatrix(AUxb);
        freecMatrix(AUzb);
        freecMatrix(IUzb);
        freecMatrix(IUxb);
        freecMatrix(IAUxb);
        freecMatrix(IAUzb);
        freecMatrix(grad);
        freecMatrix(GUxb);
        freecMatrix(GUzb);
        freecMatrix(GIUxb);
        freecMatrix(GIUzb);
        freecMatrix(HUxb);
        freecMatrix(HUzb);
        freecMatrix(HAUxb);
        freecMatrix(HAUzb);
        return (NO_MEM);
    }

    return (NO_ERR);
}

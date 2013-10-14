/***********************************************************
New piece of code to solve the incremental adjoint  equation for
 (Kx, Kz)=(0,0). The incremental state equation will be discretized exactly
as the state equation.
***********************************************************************/

#include <math.h>
#include <stdio.h>
#include <string.h>
#include "minChnl.h"
#include "mvOps.h"

/* project when (Kx,Kz) = (0,0) */
void increadjproject0(int k, int n, int count,func_force_t force )
{
    /* External Variables */
  extern int qpts, dimR, Nx, Nz;
    extern double dt, re, mpg;
    extern mcomplex *Ifa, *Ifb, *Itm, ishear;
    extern double **R, **Rp, **Rw,**Rpw,  **Rs, **Rps, *Rp0, *W;
    extern double **MZ;
    extern mcomplex **Uxb, **Uzb;
    extern mcomplex **Uxbt, **Uzbt;
    extern mcomplex **IAUxb, **IAUzb;
    extern double *Qy;
    extern mcomplex ****U, ****C;
    extern mcomplex ****IU, ****IC;
    extern mcomplex ****AU, ****AC;
    extern mcomplex ****IAU, ****IAC;
    extern mcomplex *****MC, *****MIC; 
    extern mcomplex ****LU, ****LIU;
    //static double e[3] = { 0, 1./3, 7./15};
    static double e[3] = { 1./3., 1./2, 1.};
    /* Local variables */
    int i, j;
    mcomplex tmp[dimR], tmp2[dimR], add[dimR];
    double flux_t;
    static double a[3] = {1./3., -1./2, 1./3.};
    static double b[3] = {2./3, 1./6, 0.};
    static double c[3] = { 1., 1./3., 1./2.,};
    static double d[3] = { 0.,  -2./3, -1./6};


    if(force !=NULL)
      {
	force(n, k,1, tmp, tmp2);
      }

    /* Create matrices for solving linear system. 
       Right hand side of system:  If this is the first step in the Runge-Kutta
       scheme, compute 
          [M0 + (1/RE)a[k]dt*D0]*C[0][ALPHA] + dt*c[k]Fu. 
          [M0 + (1/RE)a[k]dt*D0]*C[0][BETA] + dt*c[k]Fw. 
       Otherwise, compute C[0][ALPHA] + dt*c[k]Fu, 
                          C[0][BETA] + dt*c[k]Fw, 
       Lhs matrix:  M0 - (1/RE)b[k]dt*D0 				*/

    if (k == 0)		/* first step */
    {
        for (i = 0; i < dimR; ++i)
        {					/* MZ = M0 + (1/RE)a[k]dt*D0 */
            for (j = 0; j < T_RSDIAG; ++j)
	    {
                MZ[i][j] = Rs[i][j] - re*a[0]*dt*Rps[i][j];
            }
        }

        /* for alpha */
        smMult0(MZ, IAC[0][ALPHA], Itm, RSDIAG-1, RSDIAG-1, dimR);
        for (i = 0; i < dimR; ++i)
	{
	    Re(IAC[0][ALPHA][i][0]) = Re(Itm[i]);
	    Im(IAC[0][ALPHA][i][0]) = Im(Itm[i]);
	}

	smMult0(Rps, MIC[count][0][ALPHA], Itm, RSDIAG-1, RSDIAG-1, dimR);
	for (i = 0; i < dimR; ++i)
	{
	  Re(IAC[0][ALPHA][i][0]) -= Re(Itm[i])*dt*re*(a[0]+b[0]);
	  Im(IAC[0][ALPHA][i][0]) -= Im(Itm[i])*dt*re*(a[0]+b[0]);
	}

	/*beta */
        smMult0(MZ, IAC[0][BETA], Itm, RSDIAG-1, RSDIAG-1, dimR);
        for (i = 0; i < dimR; ++i)
	{
	    Re(IAC[0][BETA][i][0]) = Re(Itm[i]);
	    Im(IAC[0][BETA][i][0]) = Im(Itm[i]);
	}
	smMult0(Rps, MIC[count][0][BETA], Itm, RSDIAG-1, RSDIAG-1, dimR);
        for (i = 0; i < dimR; ++i)
	{
	  Re(IAC[0][BETA][i][0]) -= Re(Itm[i])*dt*re*(a[0]+b[0]);
	  Im(IAC[0][BETA][i][0]) -= Im(Itm[i])*dt*re*(a[0]+b[0]);
	}
    }
 
 /* MZ = M0 - (1/RE)b[k]dt*D0 */
    for (i = 0; i < dimR; ++i)	
    {
        for (j = 0; j < T_RSDIAG; ++j)
        {
            MZ[i][j] = Rs[i][j] + re*b[k]*dt*Rps[i][j];
        }
    }

    memset(add, 0, dimR*sizeof(mcomplex));
    for (i = 0; i < dimR; ++i)
    {
      for (j = 0; j < qpts; ++j)
	{
	  Re(add[i]) += Rw[i][j];
	}
    }
    /* Compute add  */
    bsolve_0(MZ, add, RSDIAG-1, RSDIAG-1, dimR);


    /* MZ = M0 - (1/RE)b[k]dt*D0 */
    for (i = 0; i < dimR; ++i)	
    {
        for (j = 0; j < T_RSDIAG; ++j)
        {
            MZ[i][j] = Rs[i][j] + re*b[k]*dt*Rps[i][j];
        }
    }

    /* Finish computing the right hand size */
    /* array fa */
    memset(Ifa, 0, dimR*sizeof(mcomplex));
    for (i = 0; i < dimR; ++i)
    {
      for (j = 0; j < qpts; ++j)
	{
	  Re(Ifa[i]) += Rw[i][j]*(Re(IU[0][HXEL][j][0]))+Rpw[i][j]*Re(IU[0][DZEL][j][0]);
	}
      if(force !=NULL)
	  {
	    Re(Ifa[i])+=Re(tmp[i]) ;
	    Im(Ifa[i])+=Im(tmp[i]);
	  }

      	Re(IAC[0][ALPHA][i][0]) += dt*c[k]*Re(Ifa[i]);
    }

   for (i = 0; i < dimR; ++i)
    {
      for (j = 0; j < qpts; ++j)
	{
	  Re(IAC[0][ALPHA][i][0]) += Rw[i][j]*(1-Qy[j])*0.5*(Re(Uxbt[0][0])-Re(Uxb[0][0]));
	  Im(IAC[0][ALPHA][i][0]) += Rw[i][j]*(1-Qy[j])*0.5*(Im(Uxbt[0][0])-Im(Uxb[0][0]));
	}
    }

    /* Compute a's */
    bsolve0(MZ, IAC[0][ALPHA], RSDIAG-1, RSDIAG-1, dimR);


    /* MZ = M0 - (1/RE)b[k]dt*D0 */
    for (i = 0; i < dimR; ++i)
    {
      for (j = 0; j < T_RSDIAG; ++j)
        {
            MZ[i][j] = Rs[i][j] + re*b[k]*dt*Rps[i][j];
        }
    }

    /* array fb */
    memset(Ifb, 0, dimR*sizeof(mcomplex));
    for (i = 0; i < dimR; ++i)
    {
      for (j = 0; j < qpts; ++j)
	{
	  Re(Ifb[i]) += Rw[i][j] * Re(IU[0][HZEL][j][0])-Rpw[i][j]*Re(IU[0][DXEL][j][0]);
	}
      if(force !=NULL)
	{
	  Re(Ifb[i])+=Re(tmp2[i]);
	  Im(Ifb[i])+=Im(tmp2[i]);
	  }

      Re(IAC[0][BETA][i][0]) += dt*c[k]*Re(Ifb[i]);
    }

    for (i = 0; i < dimR; ++i)
      {
	for (j = 0; j < qpts; ++j)
	  {
	    Re(IAC[0][BETA][i][0]) += Rw[i][j]*(1-Qy[j])*0.5*(Re(Uzbt[0][0])-Re(Uzb[0][0]));
	    Im(IAC[0][BETA][i][0]) += Rw[i][j]*(1-Qy[j])*0.5*(Im(Uzbt[0][0])-Im(Uzb[0][0]));
	}
    }

    /* Compute b's */
    bsolve0(MZ, IAC[0][BETA], RSDIAG-1, RSDIAG-1, dimR);


    /* NOW COMPUTE U HATS */
    /* ux_hat = R*a */
    for (i = 0; i < qpts; ++i)
    {
	Re(IAU[0][XEL][i][0]) = 0.0;
	Im(IAU[0][XEL][i][0]) = 0.0;
        for (j = 0; j < dimR; ++j)
	{
	    Re(IAU[0][XEL][i][0]) += R[i][j]*Re(IAC[0][ALPHA][j][0]);
	    Im(IAU[0][XEL][i][0]) += R[i][j]*Im(IAC[0][ALPHA][j][0]);
	}
	Re(IAU[0][XEL][i][0]) +=Re(Uxb[0][0])*(1-Qy[i])*0.5;
 	Im(IAU[0][XEL][i][0]) +=Im(Uxb[0][0])*(1-Qy[i])*0.5;
    }

    flux_t=0;
    for (i = 0; i < qpts; ++i)
      {
	flux_t += W[i]*Re(IAU[0][XEL][i][0]);
      }

    flux_t=-flux_t;

     for (j = 0; j < dimR; ++j)
      {
	Re(IAC[0][ALPHA][j][0])= Re(IAC[0][ALPHA][j][0])+flux_t*Re(add[j]);
	Im(IAC[0][ALPHA][j][0])= Im(IAC[0][ALPHA][j][0])+flux_t*Im(add[j]);
	}

    for (i = 0; i < qpts; ++i)
    {
	Re(IAU[0][XEL][i][0]) = 0.0;
	Im(IAU[0][XEL][i][0]) = 0.0;
        for (j = 0; j < dimR; ++j)
	{
	    Re(IAU[0][XEL][i][0]) += R[i][j]*Re(IAC[0][ALPHA][j][0]);
	    Im(IAU[0][XEL][i][0]) += R[i][j]*Im(IAC[0][ALPHA][j][0]);
	}
	Re(IAU[0][XEL][i][0]) +=Re(Uxb[0][0])*(1-Qy[i])*0.5;
 	Im(IAU[0][XEL][i][0]) +=Im(Uxb[0][0])*(1-Qy[i])*0.5;
    }



    /* dux_hat = Rp*a */
    for (i = 0; i < qpts; ++i)
    {
	Re(IAU[0][DXEL][i][0]) = 0.0;
	Im(IAU[0][DXEL][i][0]) = 0.0;
        for (j = 0; j < dimR; ++j)
	{
	    Re(IAU[0][DXEL][i][0]) += Rp[i][j]*Re(IAC[0][ALPHA][j][0]);
	    Im(IAU[0][DXEL][i][0]) += Rp[i][j]*Im(IAC[0][ALPHA][j][0]);
	}
	Re(IAU[0][DXEL][i][0]) += -Re(Uxb[0][0])/2.;
	Im(IAU[0][DXEL][i][0]) += -Im(Uxb[0][0])/2.;

    }
    /* uy_hat = 0 */
    for (i = 0; i < qpts; ++i)
    {
	Re(IAU[0][YEL][i][0]) = 0.0;
	Im(IAU[0][YEL][i][0]) = 0.0;
    }

    /* uz_hat = R*b */
    for (i = 0; i < qpts; ++i)
    {
	Re(IAU[0][ZEL][i][0]) = 0.0;
	Im(IAU[0][ZEL][i][0]) = 0.0;
        for (j = 0; j < dimR; ++j)
	{
	    Re(IAU[0][ZEL][i][0]) += R[i][j]*Re(IAC[0][BETA][j][0]);
	    Im(IAU[0][ZEL][i][0]) += R[i][j]*Im(IAC[0][BETA][j][0]);
	}
	Re(IAU[0][ZEL][i][0]) +=Re(Uzb[0][0])*(1-Qy[i])*0.5;
	Im(IAU[0][ZEL][i][0]) +=Im(Uzb[0][0])*(1-Qy[i])*0.5;
    }

    /* duz_hat = Rp*b */
    for (i = 0; i < qpts; ++i)
    {
	Re(IAU[0][DZEL][i][0]) = 0.0;
	Im(IAU[0][DZEL][i][0]) = 0.0;
        for (j = 0; j < dimR; ++j)
	{
	    Re(IAU[0][DZEL][i][0]) += Rp[i][j]*Re(IAC[0][BETA][j][0]);
	    Im(IAU[0][DZEL][i][0]) += Rp[i][j]*Im(IAC[0][BETA][j][0]);
	}
	Re(IAU[0][DZEL][i][0]) += -Re(Uzb[0][0])/2.;
 	Im(IAU[0][DZEL][i][0]) += -Im(Uzb[0][0])/2.;
    }

     memset(IAUxb[0], 0, Nz*(Nx/2)*sizeof(fftw_complex));
    memset(IAUzb[0], 0, Nz*(Nx/2)*sizeof(fftw_complex));
    for (j=0; j< dimR; ++j)
      {
	Re(IAUxb[0][0]) += Rp0[j]*Re(IAC[0][ALPHA][j][0]);
	Im(IAUxb[0][0]) += Rp0[j]*Im(IAC[0][ALPHA][j][0]);

	Re(IAUzb[0][0]) += Rp0[j]*Re(IAC[0][BETA][j][0]);
	Im(IAUzb[0][0]) += Rp0[j]*Im(IAC[0][BETA][j][0]);
      }

    Re(IAUxb[0][0])  -=Re(Uxb[0][0])/2.;
    Im(IAUxb[0][0])  -=Im(Uxb[0][0])/2.;
    Re(IAUzb[0][0])  -=Re(Uzb[0][0])/2.;
    Im(IAUzb[0][0])  -=Im(Uzb[0][0])/2.;

    Re(ishear)=0;
    Im(ishear)=0;
  for (j=0; j< dimR; ++j)
      {
	Re(ishear) += Rp0[j]*Re(IAC[0][ALPHA][j][0]);
	Im(ishear) += Rp0[j]*Im(IAC[0][ALPHA][j][0]);
      }
	Re(ishear) += -Re(Uxb[0][0])/2.;
	Im(ishear) += -Im(Uxb[0][0])/2.;


    /* UPDATE RHS FOR NEXT TIME */
    if (k != 2)		/* not last step */
    {
      
    if(force !=NULL)
      {
	force(n, k+1,1, tmp, tmp2);
      }
        /* MZ = M0 + (1/RE)a[k+1]dt*D0 */
        for (i = 0; i < dimR; ++i)
        {
            for (j = 0; j < T_RSDIAG; ++j)
	    {
                MZ[i][j] = Rs[i][j] - re*a[k+1]*dt*Rps[i][j];
            }
        }

        /* compute [M0 + (1/RE)a[k+1]dt*D0]*C[0][ALPHA].  Then update
	C[0][ALPHA].  Do the same for the b's. */
        smMult0(MZ, IAC[0][ALPHA], Itm, RSDIAG-1, RSDIAG-1, dimR);

	memset(Ifa, 0, dimR*sizeof(mcomplex));
	for (i = 0; i < dimR; ++i)
	  {
	    for (j = 0; j < qpts; ++j)
	      {
		Re(Ifa[i]) += Rw[i][j]*(Re(LIU[0][HXEL][j][0]))+Rpw[i][j]*Re(LIU[0][DZEL][j][0]);
	      }
	    if(force !=NULL)
	      {
		Re(Ifa[i])+=Re(tmp[i]) ;
		Im(Ifa[i])+=Im(tmp[i]) ;
	      }

	  }
	for (i = 0; i < dimR; ++i)
	{
	    Re(IAC[0][ALPHA][i][0]) = Re(Itm[i]) + dt*d[k+1]*Re(Ifa[i]);
	    Im(IAC[0][ALPHA][i][0]) = Im(Itm[i]) + dt*d[k+1]*Im(Ifa[i]);
	}
	smMult0(Rps, MIC[count-1][0][ALPHA], Itm, RSDIAG-1, RSDIAG-1, dimR);
	for (i = 0; i < dimR; ++i)
	{
	  Re(IAC[0][ALPHA][i][0]) -= Re(Itm[i])*dt*re*(a[k+1]+b[k+1]);
	  Im(IAC[0][ALPHA][i][0]) -= Im(Itm[i])*dt*re*(a[k+1]+b[k+1]);
	}

	/* now b's */
        smMult0(MZ, IAC[0][BETA], Itm, RSDIAG-1, RSDIAG-1, dimR);
	memset(Ifb, 0, dimR*sizeof(mcomplex));
	for (i = 0; i < dimR; ++i)
	  {
	    for (j = 0; j < qpts; ++j)
	      {
		Re(Ifb[i]) += Rw[i][j] * Re(LIU[0][HZEL][j][0])-Rpw[i][j]*Re(LIU[0][DXEL][j][0]);
	      }
	    if(force !=NULL)
	      {
		Re(Ifb[i])+=Re(tmp2[i]) ;
		Im(Ifb[i])+=Im(tmp2[i]) ;
	      }

	  }


	for (i = 0; i < dimR; ++i)
	{
	    Re(IAC[0][BETA][i][0]) = Re(Itm[i]) + dt*d[k+1]*Re(Ifb[i]);
	    Im(IAC[0][BETA][i][0]) = Im(Itm[i]) + dt*d[k+1]*Im(Ifb[i]);
	}
	smMult0(Rps, MIC[count-1][0][BETA], Itm, RSDIAG-1, RSDIAG-1, dimR);
        for (i = 0; i < dimR; ++i)
	{
	  Re(IAC[0][BETA][i][0]) -= Re(Itm[i])*dt*re*(a[k+1]+b[k+1]);
	  Im(IAC[0][BETA][i][0]) -= Im(Itm[i])*dt*re*(a[k+1]+b[k+1]);
	}
    }
}    /* end project zero (Kx,Kz) */

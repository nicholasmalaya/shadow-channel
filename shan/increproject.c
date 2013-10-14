/***************************************************************
New piece of code used to solve the linear system for (kx, kx)!=(0,0)

***********************************************************************/



#include <stdio.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include "minChnl.h"
#include "mvOps.h"

/* project when (Kx,Kz) != (0,0) */
void increproject(int count, int k, int z, int x0, int n,  func_force_t force)
{

  /* External Variables */
  extern int qpts, dimR, dimQ, Nx, Nz;
    extern double dt, re;
    extern double *Kx, *Kz, **K2;
    extern mcomplex **IFa, **IFb, **ITM;
    extern double **Q, **Qp, **Qpp, **R, **Rp, **Qw, **Qpw, **Rw, **Qs, **Qps,
      **Qpps, **Rs, **Rps;
    extern double ***M;
    extern mcomplex ****IU, ****IC;
    extern mcomplex **Uxb, **Uzb, **Uxbt, **Uzbt;
    extern double *Uadd, *Vadd, *Vpadd, *Qy;
    extern mcomplex *****MIC;
    /* Local variables */
    int i, j, x;
    double s, t[2];

    /* static double a[3] = {29./96., -3./40.,   1./6.};
     static double b[3] = {37./160., 5./24.,   1./6.};
    static double c[3] = { 8./15.,  5./12.,   3./4.};
    static double d[3] = { 0.,    -17./60.,  -5./12.};*/

     static double a[3] = {1./3., -1./2, 1./3.};
    static double b[3] = {1./6, 2./3, 0.};
    static double c[3] = { 1./2, 1./3., 1.,};
    static double d[3] = { 0.,  -1./6, -2./3};
    
    static double e[4] = { 4./3., -4./15., -4./15., 4./35.}; 
    static double r[8] = { -16./15., 16./35., 32./105., -32./105., -16./315., 16./231.,0., 0.};
    static double xx[8] = { 8./35., -8./315., -8./105., 8./385., 8./385., -8./1001., -8./3003., 8./6435.};
    //   static  double h[3] = { 0., 1./2., 2./3.};
    //   static  double h[3] = { 0., 8./15., 2./3.};

    mcomplex tmp[Nx/2][dimR], tmp2[Nx/2][dimQ];
     if(force!=NULL)
      {
	force(n, k, z, tmp[0], tmp2[0]);
      }

  /* FIRST COMPUTE ALPHAS */
    /* Create matrices for solving linear system. 
       Right hand side of system:  If this is the first step in the Runge-Kutta
       scheme, compute 
          [Mv + (1/RE)a[k]dt*Dv]*C[z][ALPHA] + dt*c[k]Fv. 
       Otherwise, compute C[z][ALPHA] + dt*c[k]Fv. 
       Lhs matrix:  Mv - (1/RE)b[k]dt*Dv           */

    if (k == 0)		/* first step */
    {
        memset(M[0][0], 0, dimR*9*(Nx/2)*sizeof(double));
        for (i = 0; i < dimQ; ++i) 		/* M = Mv + (1/RE)a[k]dt*Dv */
        {
            for (j = 0; j < T_QSDIAG; ++j)
	    {
                for (x = x0; x < Nx/2; ++x)
	        {
	            s = K2[z][x]*K2[z][x];
                    M[i][j][x] = -(K2[z][x]*Qs[i][j] + Qps[i][j]) +
	              re*a[0]*dt*(s*Qs[i][j] + 2.*K2[z][x]*Qps[i][j] +
		      Qpps[i][j]); 
	        }
            }
        }
        /* compute [Mv + (1/RE)a[k]dt*Dv]*IC[z][ALPHA] and store the result
	   in TM.  Then transfer the result back to IC[z][ALPHA]. */
        smMult(M, IC[z][ALPHA], ITM, QSDIAG-1, QSDIAG-1, dimQ, Nx/2, x0);
	for (i = 0; i < dimQ; ++i)
	{
	    memcpy(&IC[z][ALPHA][i][x0], &ITM[i][x0], (Nx/2-x0)*sizeof(mcomplex));
	}
    }

 /* Left hand side M = Mv - (1/RE)b[k]dt*Dv */
    memset(M[0][0], 0, dimR*9*(Nx/2)*sizeof(double));
    for (i = 0; i < dimQ; ++i)
    {
        for (j = 0; j < T_QSDIAG; ++j)
        {
            for (x = x0; x < Nx/2; ++x)
            {
                s = K2[z][x]*K2[z][x];
                M[i][j][x] = -(K2[z][x]*Qs[i][j] + Qps[i][j]) -
	          re*b[k]*dt*(s*Qs[i][j] + 2.*K2[z][x]*Qps[i][j] + Qpps[i][j]); 
	    }
        }
    }
  
    /* Finish computing the right hand size */
    /* array IFa */
    memset(IFa[0], 0, dimR*(Nx/2)*sizeof(mcomplex));
    for (i = 0; i < dimQ; ++i)
    {
        for (j = 0; j < qpts; ++j)
	{
	    for (x = x0; x < Nx/2; ++x)
	    {
	        Re(IFa[i][x]) += ( Qpw[i][j] *
	          (-Kx[x]*Im(IU[z][HXEL][j][x]) - Kz[z]*Im(IU[z][HZEL][j][x])) -
	          K2[z][x]*(Qw[i][j]*Re(IU[z][HYEL][j][x])) );
	        Im(IFa[i][x]) += ( Qpw[i][j] *
		  (Kx[x]*Re(IU[z][HXEL][j][x]) + Kz[z]*Re(IU[z][HZEL][j][x])) -
	           K2[z][x]*(Qw[i][j]*Im(IU[z][HYEL][j][x])) );
	    }
	}
    }

    if(force!=NULL)
      {
	for (x = x0; x < Nx/2; ++x)
	    {
	      for (i = 0; i < dimQ; ++i)
		{
		  Im(IFa[i][x]) += Im(tmp2[x][i]);
		  Re(IFa[i][x]) += Re(tmp2[x][i]);
		}
	    }
      }

    for (i = 0; i < dimQ; ++i)
    {
	for (x = x0; x < Nx/2; ++x)
	{
            Re(IC[z][ALPHA][i][x]) += dt*c[k]*Re(IFa[i][x]);
            Im(IC[z][ALPHA][i][x]) += dt*c[k]*Im(IFa[i][x]);
	}
    }

   
   for (x = x0; x < Nx/2; ++x)
     {
       for (i=0; i<8 &&i <dimQ; ++i)
	 {
	   Re(IC[z][ALPHA][i][x]) +=(r[i]/2.-K2[z][x]*xx[i])*(Re(Uzbt[z][x])-Re(Uzb[z][x]))+
	     dt*re*(a[k]*Re(Uzbt[z][x])+b[k]*Re(Uzb[z][x]))*(-K2[z][x]*r[i]+K2[z][x]*K2[z][x]*xx[i]);
	   Im(IC[z][ALPHA][i][x]) +=(r[i]/2.-K2[z][x]*xx[i])*(Im(Uzbt[z][x])-Im(Uzb[z][x]))+
	     dt*re*(a[k]*Im(Uzbt[z][x])+b[k]*Im(Uzb[z][x]))*(-K2[z][x]*r[i]+K2[z][x]*K2[z][x]*xx[i]);
	 }

     }

    /* Compute alphas */
    bsolve(M, IC[z][ALPHA], QSDIAG-1, QSDIAG-1, dimQ, Nx/2, x0);
 
  /* NOW COMPUTE BETAS */
    /* Create matrices for solving linear system. 
       Right hand side of system:  If this is the first step in the Runge-Kutta
       scheme, compute 
          [Mg + (1/RE)a[k]dt*Dg]*C[z][BETA] + dt*c[k]Fg. 
       Otherwise, compute C[z][BETA] + dt*c[k]Fg. 
       Lhs matrix:  Mg - (1/RE)b[k]dt*Dg 				*/

    if (k == 0)		/* first step */
    {
        /* M = Mg + (1/RE)a[k]dt*Dg */
        memset(M[0][0], 0, dimR*9*(Nx/2)*sizeof(double));
        for (i = 0; i < dimR; ++i)	
        {
            for (j = 0; j < T_RSDIAG; ++j)
            {
		for (x = x0; x < Nx/2; ++x)
		{
                    M[i][j][x] = Rs[i][j] - 
		      re*a[0]*dt*(Rps[i][j] + K2[z][x]*Rs[i][j]);
		}
            }
	}

        /* compute [Mv + (1/RE)a[k]dt*Dv]*IC[z][BETA] and store the result
	   in TM.  Then transfer the result back to IC[z][ALPHA]. */
        smMult(M, IC[z][BETA], ITM, RSDIAG-1, RSDIAG-1, dimR, Nx/2, x0);
	for (i = 0; i < dimR; ++i)
	{
	    memcpy(&IC[z][BETA][i][x0], &ITM[i][x0], (Nx/2-x0)*sizeof(mcomplex));
	}
    }

    /* left hand side M = Mg - (1/RE)b[k]dt*Dg */
    memset(M[0][0], 0, dimR*9*(Nx/2)*sizeof(double));
    for (i = 0; i < dimR; ++i)		/* M = Mg - (1/RE)b[k]dt*Dg */
    {
        for (j = 0; j < T_RSDIAG; ++j)
        {
	    for (x = x0; x < Nx/2; ++x)
	    {
                    M[i][j][x] = Rs[i][j] + 
		      re*b[k]*dt*(Rps[i][j] + K2[z][x]*Rs[i][j]);
            }
	}
    }

   /* Finish computing the right hand size */
    /* array Fb */
    memset(IFb[0], 0, dimR*(Nx/2)*sizeof(mcomplex));
    for (i = 0; i < dimR; ++i)
    {
        for (j = 0; j < qpts; ++j)
	{
	    for (x = x0; x < Nx/2; ++x)
	    {
	        Re(IFb[i][x]) += Rw[i][j] *
		  (-Kz[z]*Im(IU[z][HXEL][j][x]) + Kx[x]*Im(IU[z][HZEL][j][x]));
	        Im(IFb[i][x]) += Rw[i][j] *
		  (Kz[z]*Re(IU[z][HXEL][j][x]) - Kx[x]*Re(IU[z][HZEL][j][x]));
	    }
	}
    }
 
    if(force!=NULL)
      {
	for (x = x0; x < Nx/2; ++x)
	  {
	    for (i = 0; i < dimR; ++i)
	      {
		Im(IFb[i][x]) += Im(tmp[x][i]);
		Re(IFb[i][x]) += Re(tmp[x][i]);
	      }
	  }
      }
    for (i = 0; i < dimR; ++i)
    {
	for (x = x0; x < Nx/2; ++x)
	{
	    Re(IC[z][BETA][i][x]) += dt*c[k]*Re(IFb[i][x]);
	    Im(IC[z][BETA][i][x]) += dt*c[k]*Im(IFb[i][x]);
	}
    }

    for (x = x0; x < Nx/2; ++x)
	  {
	    for (i=0; i<4 && i<dimR; ++i)
	      {
		Re(IC[z][BETA][i][x]) +=e[i]/2.*(Re(Uxbt[z][x])-Re(Uxb[z][x]))+
		  dt*re*(a[k]*Re(Uxbt[z][x])+b[k]*Re(Uxb[z][x]))*(-K2[z][x]*e[i]/2.);
		Im(IC[z][BETA][i][x]) +=e[i]/2.*(Im(Uxbt[z][x])-Im(Uxb[z][x]))+
		  dt*re*(a[k]*Im(Uxbt[z][x])+b[k]*Im(Uxb[z][x]))*(-K2[z][x]*e[i]/2.);
	      }
	  }

    /* Compute betas */
    bsolve(M, IC[z][BETA], RSDIAG-1, RSDIAG-1, dimR, Nx/2, x0);
 
    /* NOW COMPUTE IU HATS */
    /* v = uy_hat+c2/4*(1-y)^2*(1+y). */
    for (i = 0; i < qpts; ++i)
    {
        memset(&IU[z][YEL][i][x0], 0, (Nx/2-x0)*sizeof(mcomplex));
    }
    for (i = 0; i < qpts; ++i)
    {
      for (x = x0; x < Nx/2; ++x)
	{
	  for (j = 0; j < dimQ; ++j)
	    {
	        Re(IU[z][YEL][i][x]) += Q[i][j]*Re(IC[z][ALPHA][j][x]);
	        Im(IU[z][YEL][i][x]) += Q[i][j]*Im(IC[z][ALPHA][j][x]);
	    }
	  Re(IU[z][YEL][i][x]) += Re(Uzb[z][x])*Vadd[i]/4.;
	  Im(IU[z][YEL][i][x]) += Im(Uzb[z][x])*Vadd[i]/4.;
	}
    }


   /* f = -dv/dy and store temporarily in XEL position of array U. */
    /*f=-dv/dy-c2/4*(3y^2-2y-1) */
    for (i = 0; i < qpts; ++i)
    {
        memset(&IU[z][XEL][i][x0], 0, (Nx/2-x0)*sizeof(mcomplex));
    }
    for (i = 0; i < qpts; ++i)
    {
      for (x = x0; x < Nx/2; ++x)
	{
	  for (j = 0; j < dimQ; ++j)
	    {
	        Re(IU[z][XEL][i][x]) -= Qp[i][j]*Re(IC[z][ALPHA][j][x]);
	        Im(IU[z][XEL][i][x]) -= Qp[i][j]*Im(IC[z][ALPHA][j][x]);
	    }
	  Re(IU[z][XEL][i][x]) -= Re(Uzb[z][x])*Vpadd[i]/4.;
	  Im(IU[z][XEL][i][x]) -= Im(Uzb[z][x])*Vpadd[i]/4.;
	}
    }

    /* sum(Q''alpha) and store in DXEL position. */
    /*df/dy=d^2 v/dy^2+c2/4*(6y-2)=d^2 v/dy^2+c2/2*(-3+3y+2)*/
    for (i = 0; i < qpts; ++i)
    {
        memset(&IU[z][DXEL][i][x0], 0, (Nx/2-x0)*sizeof(mcomplex));
    }
    for (i = 0; i < qpts; ++i)
    {
      for (x = x0; x < Nx/2; ++x)
	{
	  for (j = 0; j < dimQ; ++j)
	    {
	        Re(IU[z][DXEL][i][x]) += Qpp[i][j]*Re(IC[z][ALPHA][j][x]);
	        Im(IU[z][DXEL][i][x]) += Qpp[i][j]*Im(IC[z][ALPHA][j][x]);
	    }
	  Re(IU[z][DXEL][i][x]) += Re(Uzb[z][x])*(-3*Uadd[i]+2)/2.;
	  Im(IU[z][DXEL][i][x]) += Im(Uzb[z][x])*(-3*Uadd[i]+2)/2.;
	}
    }

   /* Compute g = sum(beta*R) and store temporarily in ZEL position of 
       array U. */
    for (i = 0; i < qpts; ++i)
    {
        memset(&IU[z][ZEL][i][x0], 0, (Nx/2-x0)*sizeof(mcomplex));
    }
    for (i = 0; i < qpts; ++i)
    {
      for (x = x0; x < Nx/2; ++x)
	{
	  for (j = 0; j < dimR; ++j)
	    {
	        Re(IU[z][ZEL][i][x]) += R[i][j]*Re(IC[z][BETA][j][x]);
	        Im(IU[z][ZEL][i][x]) += R[i][j]*Im(IC[z][BETA][j][x]);
	    }

       	  Re(IU[z][ZEL][i][x]) += Re(Uxb[z][x])*Uadd[i]/2.;
	  Im(IU[z][ZEL][i][x]) += Im(Uxb[z][x])*Uadd[i]/2.;
	  }
    }


    /* Compute sum(beta*R') and store temporarily in DZEL position of
       array IU. */ 
    for (i = 0; i < qpts; ++i)
    {
        memset(&IU[z][DZEL][i][x0], 0, (Nx/2-x0)*sizeof(mcomplex));
    }
    for (i = 0; i < qpts; ++i)
    {
      for (x = x0; x < Nx/2; ++x)
	{
	  for (j = 0; j < dimR; ++j)
	    {
	        Re(IU[z][DZEL][i][x]) += Rp[i][j]*Re(IC[z][BETA][j][x]);
	        Im(IU[z][DZEL][i][x]) += Rp[i][j]*Im(IC[z][BETA][j][x]);
	    }
	  Re(IU[z][DZEL][i][x]) += -Re(Uxb[z][x])/2.;
	  Im(IU[z][DZEL][i][x]) += -Im(Uxb[z][x])/2.;
	}
    }

  /* now compute Iux hat, Iuz hat */
    for (i = 0; i < qpts; ++i)
    {
	for (x = x0; x < Nx/2; ++x)
	{
	   t[0] = Re(IU[z][XEL][i][x]);	/* real part of f */
	   t[1] = Re(IU[z][ZEL][i][x]);	/* real part of g */

	   Re(IU[z][XEL][i][x]) = (Kx[x]*Im(IU[z][XEL][i][x]) +
	      Kz[z]*Im(IU[z][ZEL][i][x])) / K2[z][x];
	   Re(IU[z][ZEL][i][x]) = (-Kx[x]*Im(IU[z][ZEL][i][x]) +
	      Kz[z]*Im(IU[z][XEL][i][x])) / K2[z][x];
	   Im(IU[z][XEL][i][x]) = -(Kx[x]*t[0] + Kz[z]*t[1]) / K2[z][x];
	   Im(IU[z][ZEL][i][x]) = (Kx[x]*t[1] - Kz[z]*t[0]) / K2[z][x];
	}
    }

    /* dux hat, duz hat */
    for (i = 0; i < qpts; ++i)
    {
	for (x = x0; x < Nx/2; ++x)
	{
	   t[0] = Re(IU[z][DXEL][i][x]);	/* real part of Q''alpha */
	   t[1] = Re(IU[z][DZEL][i][x]);	/* real part of R'beta */

	   Re(IU[z][DXEL][i][x]) = (-Kx[x]*Im(IU[z][DXEL][i][x]) +
	      Kz[z]*Im(IU[z][DZEL][i][x])) / K2[z][x];
	   Re(IU[z][DZEL][i][x]) = -(Kx[x]*Im(IU[z][DZEL][i][x]) +
	      Kz[z]*Im(IU[z][DXEL][i][x])) / K2[z][x];
	   Im(IU[z][DXEL][i][x]) = (Kx[x]*t[0] - Kz[z]*t[1]) / K2[z][x];
	   Im(IU[z][DZEL][i][x]) = (Kx[x]*t[1] + Kz[z]*t[0]) / K2[z][x];
	}
    }
 

    for (i=0; i< dimQ; i++)
      {
	for (x=x0; x<Nx/2; x++)
	  {
	    Re(MIC[count][z][ALPHA][i][x])=Re(IC[z][ALPHA][i][x]);
	    Im(MIC[count][z][ALPHA][i][x])=Im(IC[z][ALPHA][i][x]);
	  }
      }
    for (i=0; i< dimR; i++)
      {
	for (x=x0; x<Nx/2; x++)
	  {
	    Re(MIC[count][z][BETA][i][x])=Re(IC[z][BETA][i][x]);
	    Im(MIC[count][z][BETA][i][x])=Im(IC[z][BETA][i][x]);
	  }
      }
   /* UPDATE RHS FOR NEXT TIME */
    if (k != 2)		/* not last step */
    {
	/* first alphas */
        memset(M[0][0], 0, dimR*9*(Nx/2)*sizeof(double));
        for (i = 0; i < dimQ; ++i)		/* M = Mv + (1/RE)a[k+1]dt*Dv */
        {
            for (j = 0; j < T_QSDIAG; ++j)
	    {
                for (x = x0; x < Nx/2; ++x)
	        {
	            s = K2[z][x]*K2[z][x];
                    M[i][j][x] = -(K2[z][x]*Qs[i][j] + Qps[i][j]) +
	              re*a[k+1]*dt*(s*Qs[i][j] + 2.*K2[z][x]*Qps[i][j] +
		      Qpps[i][j]); 
		}
            }
        }

  /* compute [Mv + (1/RE)a[k]dt*Dv]*C[z][ALPHA].  Then update
	C[z][ALPHA] */
        smMult(M, IC[z][ALPHA], ITM, QSDIAG-1, QSDIAG-1, dimQ, Nx/2, x0);
	for (i = 0; i < dimQ; ++i)
	{
	    for (x = x0; x < Nx/2; ++x)
	    {
		Re(IC[z][ALPHA][i][x]) = Re(ITM[i][x]) + dt*d[k+1]*Re(IFa[i][x]);
		Im(IC[z][ALPHA][i][x]) = Im(ITM[i][x]) + dt*d[k+1]*Im(IFa[i][x]);
	    }
	}

	/* now betas */
        memset(M[0][0], 0, dimR*9*(Nx/2)*sizeof(double));
        for (i = 0; i < dimR; ++i)		/* M = Mg + (1/RE)a[k+1]dt*Dg */
        {
            for (j = 0; j < T_RSDIAG; ++j)
	    {
                for (x = x0; x < Nx/2; ++x)
	        {
                    M[i][j][x] = Rs[i][j] - 
		      re*a[k+1]*dt*(Rps[i][j] + K2[z][x]*Rs[i][j]);
	        }
            }
        }

        /* compute [Mg + (1/RE)a[k+1]dt*Dg]*C[z][BETA] and then
	   update C[z][BETA] */
        smMult(M, IC[z][BETA], ITM, RSDIAG-1, RSDIAG-1, dimR, Nx/2, x0);
	for (i = 0; i < dimR; ++i)
	{
	    for (x = x0; x < Nx/2; ++x)
	    {
		Re(IC[z][BETA][i][x]) = Re(ITM[i][x]) + dt*d[k+1]*Re(IFb[i][x]);
		Im(IC[z][BETA][i][x]) = Im(ITM[i][x]) + dt*d[k+1]*Im(IFb[i][x]);
	    }
	}
    }

}

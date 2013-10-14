
#include <stdio.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include "minChnl.h"
#include "mvOps.h"

/* project when (Kx,Kz) != (0,0) */
void increadjproject(int n, int k, int z, int x0, int count,func_force_t force )
{
    /* External Variables */
  extern int qpts, dimR, dimQ, Nx, Nz;
    extern double dt, re;
    extern double *Kx, *Kz, **K2;
    extern mcomplex **IFa, **IFb, **ITM;
    extern double **Q, **Qp, **Qpp, **R, **Rp, **Qw, **Qpw, **Rw, **Qs, **Qps,
      **Qpps, **Rs, **Rps, *Rp0, **Qppw, **Rpw, *Qy;
    extern double ***M;
    extern mcomplex **Uxb, **Uzb;
    extern mcomplex **IAUxb, **IAUzb;
    extern mcomplex **Uxbt, **Uzbt;
    extern mcomplex **AUxb, **AUzb;
    extern mcomplex ****U, ****C;
    extern mcomplex ****AU, ****AC;
    extern mcomplex ****IU, ****IC;
    extern mcomplex ****IAU, ****IAC;
    extern mcomplex *****MC, ****LIU;
    extern mcomplex *****MIC, *****MC;;
    /* Local variables */
    int i, j, x;
    double s, t[2];
    int sign;

    static double a[3] = {1./3., -1./2, 1./3.};
    static double b[3] = {2./3, 1./6, 0.};
    static double c[3] = { 1., 1./3., 1./2.,};
    static double d[3] = { 0.,  -2./3, -1./6};
    static double h[3]={ 0., 1., 2./3.};
   static double e[3] = { 1./3., 1./2, 1.};
    static double g[4] = { 4./3., -4./15., -4./15., 4./35.}; 
    static double r[8] = { -16./15., 16./35., 32./105., -32./105., -16./315., 16./231.,0., 0.};
    static double xx[8] = { 8./35., -8./315., -8./105., 8./385., 8./385., -8./1001., -8./3003., 8./6435.};
   static double ss[3] = {1.,2./3., 1.};
    mcomplex tmp[Nx/2][dimR], tmp2[Nx/2][dimQ];
    double k11=2*M_PI;

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
       Lhs matrix:  Mv - (1/RE)b[k]dt*Dv 				*/
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
        /* compute [Mv + (1/RE)a[k]dt*Dv]*C[z][ALPHA] and store the result
	   in TM.  Then transfer the result back to C[z][ALPHA]. */
        smMult(M, IAC[z][ALPHA], ITM, QSDIAG-1, QSDIAG-1, dimQ, Nx/2, x0);
	for (i = 0; i < dimQ; ++i)
	{
	    memcpy(&IAC[z][ALPHA][i][x0], &ITM[i][x0], (Nx/2-x0)*sizeof(mcomplex));
	}

	memset(M[0][0], 0, dimR*9*(Nx/2)*sizeof(double));
	for (i = 0; i < dimQ; ++i)
        {
            for (j = 0; j < T_QSDIAG; ++j)
	    {
                for (x = x0; x < Nx/2; ++x)
	        {
	            s = K2[z][x]*K2[z][x];
                    M[i][j][x] = (s*Qs[i][j] + 2.*K2[z][x]*Qps[i][j] +
		      Qpps[i][j]); 
	        }
            }
        }
	smMult(M, MIC[count][z][ALPHA], ITM, QSDIAG-1, QSDIAG-1, dimQ, Nx/2, x0);
	for (i = 0; i < dimQ; ++i)
	{
	  for (x=x0; x< Nx/2; x++)
	    {
	      Re(IAC[z][ALPHA][i][x]) += Re(ITM[i][x])*re*dt*(a[0]+b[0]);
	      Im(IAC[z][ALPHA][i][x]) += Im(ITM[i][x])*re*dt*(a[0]+b[0]); 
	    }
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
    /* array Fa */
    memset(IFa[0], 0, dimR*(Nx/2)*sizeof(mcomplex));
    for (i = 0; i < dimQ; ++i)
    {
        for (j = 0; j < qpts; ++j)
	{
	    for (x = x0; x < Nx/2; ++x)
	    {
	        Re(IFa[i][x]) += ( Qpw[i][j] *
	          (-Kx[x]*Im(IU[z][HXEL][j][x]) - Kz[z]*Im(IU[z][HZEL][j][x])) -
	          K2[z][x]*(Qw[i][j]*Re(IU[z][HYEL][j][x])) )+ Qppw[i][j]*
		  (Kz[z]*Im(IU[z][DXEL][j][x])-Kx[x]*Im(IU[z][DZEL][j][x]));

	        Im(IFa[i][x]) += ( Qpw[i][j] *
		  (Kx[x]*Re(IU[z][HXEL][j][x]) + Kz[z]*Re(IU[z][HZEL][j][x])) -
	           K2[z][x]*(Qw[i][j]*Im(IU[z][HYEL][j][x])) )+ Qppw[i][j]*
		  (-Kz[z]*Re(IU[z][DXEL][j][x])+Kx[x]*Re(IU[z][DZEL][j][x]));
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
            Re(IAC[z][ALPHA][i][x]) += dt*c[k]*Re(IFa[i][x]);
            Im(IAC[z][ALPHA][i][x]) += dt*c[k]*Im(IFa[i][x]);
	}
    }
   for (x = x0; x < Nx/2; ++x)
     {
       for (i=0; i<8 &&i <dimQ; ++i)
	 {
	   Re(IAC[z][ALPHA][i][x]) +=(r[i]/2.-K2[z][x]*xx[i])*(Re(Uzbt[z][x])-Re(Uzb[z][x]))+
	     dt*re*(a[k]*Re(Uzbt[z][x])+b[k]*Re(Uzb[z][x]))*(-K2[z][x]*r[i]+K2[z][x]*K2[z][x]*xx[i])
	     +dt*re*(a[k]+b[k])*Re(AUzb[z][x])*(-K2[z][x]*r[i]+K2[z][x]*K2[z][x]*xx[i]);
	   Im(IAC[z][ALPHA][i][x]) +=(r[i]/2.-K2[z][x]*xx[i])*(Im(Uzbt[z][x])-Im(Uzb[z][x]))+
	     dt*re*(a[k]*Im(Uzbt[z][x])+b[k]*Im(Uzb[z][x]))*(-K2[z][x]*r[i]+K2[z][x]*K2[z][x]*xx[i])
	     +  dt*re*(a[k]+b[k])*Im(AUzb[z][x])*(-K2[z][x]*r[i]+K2[z][x]*K2[z][x]*xx[i]);
	 }

     }

    /* Compute alphas */
    bsolve(M, IAC[z][ALPHA], QSDIAG-1, QSDIAG-1, dimQ, Nx/2, x0);

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

        /* compute [Mv + (1/RE)a[k]dt*Dv]*C[z][BETA] and store the result
	   in TM.  Then transfer the result back to C[z][ALPHA]. */
        smMult(M, IAC[z][BETA], ITM, RSDIAG-1, RSDIAG-1, dimR, Nx/2, x0);
	for (i = 0; i < dimR; ++i)
	{
	    memcpy(&IAC[z][BETA][i][x0], &ITM[i][x0], (Nx/2-x0)*sizeof(mcomplex));
	}

	memset(M[0][0], 0, dimR*9*(Nx/2)*sizeof(double));
	for (i = 0; i < dimR; ++i) 	
        {
            for (j = 0; j < T_RSDIAG; ++j)
	    {
                for (x = x0; x < Nx/2; ++x)
	        {
                    M[i][j][x] = -(Rps[i][j] + K2[z][x]*Rs[i][j]); 
	        }
            }
        }

	smMult(M, MIC[count][z][BETA], ITM, RSDIAG-1, RSDIAG-1, dimR, Nx/2, x0);
	for (i = 0; i < dimR; ++i)
	{
	  for (x=x0; x< Nx/2; x++)
	    {
	      Re(IAC[z][BETA][i][x]) += Re(ITM[i][x])*re*dt*(a[0]+b[0]);
	      Im(IAC[z][BETA][i][x]) += Im(ITM[i][x])*re*dt*(a[0]+b[0]); 
	    }
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
		  (-Kz[z]*Im(IU[z][HXEL][j][x]) + Kx[x]*Im(IU[z][HZEL][j][x]))+
		  Rpw[i][j]*(-Kz[z]*Im(IU[z][DZEL][j][x])-Kx[x]*Im(IU[z][DXEL][j][x]));
	        Im(IFb[i][x]) += Rw[i][j] *
		  (Kz[z]*Re(IU[z][HXEL][j][x]) - Kx[x]*Re(IU[z][HZEL][j][x]))+
		  Rpw[i][j]*(Kz[z]*Re(IU[z][DZEL][j][x])+Kx[x]*Re(IU[z][DXEL][j][x]));
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
	    Re(IAC[z][BETA][i][x]) += dt*c[k]*Re(IFb[i][x]);
	    Im(IAC[z][BETA][i][x]) += dt*c[k]*Im(IFb[i][x]);
	}
    }
 
    for (x = x0; x < Nx/2; ++x)
	  {
	    for (i=0; i<4 && i<dimR; ++i)
	      {
		Re(IAC[z][BETA][i][x]) +=g[i]/2.*(Re(Uxbt[z][x])-Re(Uxb[z][x]))+
		  dt*re*(a[k]*Re(Uxbt[z][x])+b[k]*Re(Uxb[z][x]))*(-K2[z][x]*g[i]/2.)
		  + dt*re*(a[k]+b[k])*Re(AUxb[z][x])*(-K2[z][x]*g[i]/2.) ;
		Im(IAC[z][BETA][i][x]) +=g[i]/2.*(Im(Uxbt[z][x])-Im(Uxb[z][x]))+
		  dt*re*(a[k]*Im(Uxbt[z][x])+b[k]*Im(Uxb[z][x]))*(-K2[z][x]*g[i]/2.)
		  +  dt*re*(a[k]+b[k])*Im(AUxb[z][x])*(-K2[z][x]*g[i]/2.);
	      }
	  }


    /* Compute betas */
    bsolve(M, IAC[z][BETA], RSDIAG-1, RSDIAG-1, dimR, Nx/2, x0);
 
    /* NOW COMPUTE U HATS */
    /* v = uy_hat. */
    for (i = 0; i < qpts; ++i)
    {
        memset(&IAU[z][YEL][i][x0], 0, (Nx/2-x0)*sizeof(mcomplex));
    }

    for (i = 0; i < qpts; ++i)
      {
	for (x = x0; x < Nx/2; ++x)
	  {
	    for (j = 0; j < dimQ; ++j)
	      {
		Re(IAU[z][YEL][i][x]) += Q[i][j]*Re(IAC[z][ALPHA][j][x]);
	        Im(IAU[z][YEL][i][x]) += Q[i][j]*Im(IAC[z][ALPHA][j][x]);
	      }
	    Re(IAU[z][YEL][i][x]) += Re(Uzb[z][x])*(1-Qy[i])*(1-Qy[i])*(1+Qy[i])*0.25;
	    Im(IAU[z][YEL][i][x]) += Im(Uzb[z][x])*(1-Qy[i])*(1-Qy[i])*(1+Qy[i])*0.25;
	}
    }

    /* f = -dv/dy and store temporarily in XEL position of array U. */
    for (i = 0; i < qpts; ++i)
    {
        memset(&IAU[z][XEL][i][x0], 0, (Nx/2-x0)*sizeof(mcomplex));
    }

    for (i = 0; i < qpts; ++i)
      {
	for (x = x0; x < Nx/2; ++x)
	  {
	    for (j = 0; j < dimQ; ++j)
	      {
	        Re(IAU[z][XEL][i][x]) -= Qp[i][j]*Re(IAC[z][ALPHA][j][x]);
	        Im(IAU[z][XEL][i][x]) -= Qp[i][j]*Im(IAC[z][ALPHA][j][x]);
	      }
	    Re(IAU[z][XEL][i][x]) -= Re(Uzb[z][x])*(3*Qy[i]*Qy[i]-2*Qy[i]-1)*0.25;
	    Im(IAU[z][XEL][i][x]) -= Im(Uzb[z][x])*(3*Qy[i]*Qy[i]-2*Qy[i]-1)*0.25;
	}
    }
 
    /* sum(Q''alpha) and store in DXEL position. */
    for (i = 0; i < qpts; ++i)
    {
        memset(&IAU[z][DXEL][i][x0], 0, (Nx/2-x0)*sizeof(mcomplex));
    }

   for (i = 0; i < qpts; ++i)
      {
	for (x = x0; x < Nx/2; ++x)
	  {
	    for (j = 0; j < dimQ; ++j)
	      {
	        Re(IAU[z][DXEL][i][x]) += Qpp[i][j]*Re(IAC[z][ALPHA][j][x]);
	        Im(IAU[z][DXEL][i][x]) += Qpp[i][j]*Im(IAC[z][ALPHA][j][x]);
	    }
	    Re(IAU[z][DXEL][i][x]) += Re(Uzb[z][x])*(3*Qy[i]-1)*0.5;
	    Im(IAU[z][DXEL][i][x]) += Im(Uzb[z][x])*(3*Qy[i]-1)*0.5;
	}
    }

    /* Compute g = sum(beta*R) and store temporarily in ZEL position of 
       array U. */
    for (i = 0; i < qpts; ++i)
    {
        memset(&IAU[z][ZEL][i][x0], 0, (Nx/2-x0)*sizeof(mcomplex));
    }
    for (i = 0; i < qpts; ++i)
    {
        for (x = x0; x < Nx/2; ++x)
	  {
	    for (j = 0; j < dimQ; ++j)
	      {
	        Re(IAU[z][ZEL][i][x]) += R[i][j]*Re(IAC[z][BETA][j][x]);
	        Im(IAU[z][ZEL][i][x]) += R[i][j]*Im(IAC[z][BETA][j][x]);
	      }
	    Re(IAU[z][ZEL][i][x]) += Re(Uxb[z][x])*(1-Qy[i])*0.5;
	    Im(IAU[z][ZEL][i][x]) += Im(Uxb[z][x])*(1-Qy[i])*0.5;
	}
    }
  
    /* Compute sum(beta*R') and store temporarily in DZEL position of
       array U. */ 
    for (i = 0; i < qpts; ++i)
    {
        memset(&IAU[z][DZEL][i][x0], 0, (Nx/2-x0)*sizeof(mcomplex));
    }
    for (i = 0; i < qpts; ++i)
    {
        for (x = x0; x < Nx/2; ++x)
	  {
	    for (j = 0; j < dimQ; ++j)
	      {
	        Re(IAU[z][DZEL][i][x]) += Rp[i][j]*Re(IAC[z][BETA][j][x]);
	        Im(IAU[z][DZEL][i][x]) += Rp[i][j]*Im(IAC[z][BETA][j][x]);
	      }
	  Re(IAU[z][DZEL][i][x]) += -Re(Uxb[z][x])/2.;
	  Im(IAU[z][DZEL][i][x]) += -Im(Uxb[z][x])/2.;
	}
    }

    /* now compute ux hat, uz hat */
    for (i = 0; i < qpts; ++i)
    {
	for (x = x0; x < Nx/2; ++x)
	{
	   t[0] = Re(IAU[z][XEL][i][x]);	/* real part of f */
	   t[1] = Re(IAU[z][ZEL][i][x]);	/* real part of g */

	   Re(IAU[z][XEL][i][x]) = (Kx[x]*Im(IAU[z][XEL][i][x]) +
	      Kz[z]*Im(IAU[z][ZEL][i][x])) / K2[z][x];
	   Re(IAU[z][ZEL][i][x]) = (-Kx[x]*Im(IAU[z][ZEL][i][x]) +
	      Kz[z]*Im(IAU[z][XEL][i][x])) / K2[z][x];
	   Im(IAU[z][XEL][i][x]) = -(Kx[x]*t[0] + Kz[z]*t[1]) / K2[z][x];
	   Im(IAU[z][ZEL][i][x]) = (Kx[x]*t[1] - Kz[z]*t[0]) / K2[z][x];
	}
    }

    /* dux hat, duz hat */
    for (i = 0; i < qpts; ++i)
    {
	for (x = x0; x < Nx/2; ++x)
	{
	   t[0] = Re(IAU[z][DXEL][i][x]);	/* real part of Q''alpha */
	   t[1] = Re(IAU[z][DZEL][i][x]);	/* real part of R'beta */

	   Re(IAU[z][DXEL][i][x]) = (-Kx[x]*Im(IAU[z][DXEL][i][x]) +
	      Kz[z]*Im(IAU[z][DZEL][i][x])) / K2[z][x];
	   Re(IAU[z][DZEL][i][x]) = -(Kx[x]*Im(IAU[z][DZEL][i][x]) +
	      Kz[z]*Im(IAU[z][DXEL][i][x])) / K2[z][x];
	   Im(IAU[z][DXEL][i][x]) = (Kx[x]*t[0] - Kz[z]*t[1]) / K2[z][x];
	   Im(IAU[z][DZEL][i][x]) = (Kx[x]*t[1] + Kz[z]*t[0]) / K2[z][x];
	}
    }

   for (x=x0; x< Nx/2; ++x)
      {
	for (j=0; j < dimR; ++j)
	  {
	    Re(IAUxb[z][x]) += (Rp0[j]*Im(IAC[z][BETA][j][x]))*Kz[z]/K2[z][x];
	    Im(IAUxb[z][x]) +=(-Rp0[j]*Re(IAC[z][BETA][j][x]))*Kz[z]/K2[z][x];

	    Re(IAUzb[z][x]) +=-(Rp0[j]*Im(IAC[z][BETA][j][x]))*Kx[x]/K2[z][x];
	    Im(IAUzb[z][x]) += (Rp0[j]*Re(IAC[z][BETA][j][x]))*Kx[x]/K2[z][x];

	  }
	Re(IAUxb[z][x]) +=-Im(Uxb[z][x])/2.*Kz[z]/K2[z][x];
	Im(IAUxb[z][x]) +=Re(Uxb[z][x])/2.*Kz[z]/K2[z][x];
	Re(IAUzb[z][x]) +=Im(Uxb[z][x])/2.*Kx[x]/K2[z][x];
	Im(IAUzb[z][x]) +=-Re(Uxb[z][x])/2.*Kx[x]/K2[z][x];
      }

 for (x=x0; x< Nx/2; ++x)
      {
	sign=-1;
	for (j=0; j < dimQ; ++j)
	  {
	    sign=sign*(-1);
	    Re(IAUxb[z][x]) -= (sign*8.*Im(IAC[z][ALPHA][j][x]))*Kx[x]/K2[z][x];
	    Im(IAUxb[z][x]) -=(-sign*8.*Re(IAC[z][ALPHA][j][x]))*2*Kx[x]/K2[z][x];

	    Re(IAUzb[z][x]) -= (sign*8.*Im(IAC[z][ALPHA][j][x]))*Kz[z]/K2[z][x];
	    Im(IAUzb[z][x]) -=(-sign*8.*Re(IAC[z][ALPHA][j][x]))*2*Kz[z]/K2[z][x];
	  }
	Re(IAUxb[z][x]) -=-Im(Uzb[z][x])*2.*Kx[x]/K2[z][x];
	Im(IAUxb[z][x]) -= Re(Uzb[z][x])*2.*Kx[x]/K2[z][x];
	Re(IAUzb[z][x]) -= -Im(Uzb[z][x])*2.*Kz[z]/K2[z][x];
	Im(IAUzb[z][x]) -= Re(Uzb[z][x])*2*Kz[z]/K2[z][x];
      }

    /* UPDATE RHS FOR NEXT TIME */
    if (k != 2)		/* not last step */
    {
      if(force!=NULL)
      {
	force(n, k+1, z, tmp[0], tmp2[0]);
      }
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
        smMult(M, IAC[z][ALPHA], ITM, QSDIAG-1, QSDIAG-1, dimQ, Nx/2, x0);
	memset(IFa[0], 0, dimR*(Nx/2)*sizeof(mcomplex));
	for (i = 0; i < dimQ; ++i)
	  {
	    for (j = 0; j < qpts; ++j)
	      {
		for (x = x0; x < Nx/2; ++x)
		  {
		    Re(IFa[i][x]) += ( Qpw[i][j] *
				      (-Kx[x]*Im(LIU[z][HXEL][j][x]) - Kz[z]*Im(LIU[z][HZEL][j][x])) -
				      K2[z][x]*(Qw[i][j]*Re(LIU[z][HYEL][j][x])) )+ Qppw[i][j]*
		      (Kz[z]*Im(LIU[z][DXEL][j][x])-Kx[x]*Im(LIU[z][DZEL][j][x]));

		    Im(IFa[i][x]) += ( Qpw[i][j] *
				      (Kx[x]*Re(LIU[z][HXEL][j][x]) + Kz[z]*Re(LIU[z][HZEL][j][x])) -
				      K2[z][x]*(Qw[i][j]*Im(LIU[z][HYEL][j][x])) )+ Qppw[i][j]*
		      (-Kz[z]*Re(LIU[z][DXEL][j][x])+Kx[x]*Re(LIU[z][DZEL][j][x]));
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
		Re(IAC[z][ALPHA][i][x]) = Re(ITM[i][x]) + dt*d[k+1]*Re(IFa[i][x]);
		Im(IAC[z][ALPHA][i][x]) = Im(ITM[i][x]) + dt*d[k+1]*Im(IFa[i][x]);
	    }
	}
	memset(M[0][0], 0, dimR*9*(Nx/2)*sizeof(double));
	for (i = 0; i < dimQ; ++i) 	
        {
            for (j = 0; j < T_QSDIAG; ++j)
	    {
                for (x = x0; x < Nx/2; ++x)
	        {
	            s = K2[z][x]*K2[z][x];
                    M[i][j][x] = (s*Qs[i][j] + 2.*K2[z][x]*Qps[i][j] +
		      Qpps[i][j]); 
	        }
            }
        }
	smMult(M, MIC[count-1][z][ALPHA], ITM, QSDIAG-1, QSDIAG-1, dimQ, Nx/2, x0);

	for (i = 0; i < dimQ; ++i)
	{
	  for (x=x0; x< Nx/2; x++)
	    {
	      Re(IAC[z][ALPHA][i][x]) += Re(ITM[i][x])*re*dt*(a[k+1]+b[k+1]);
	      Im(IAC[z][ALPHA][i][x]) += Im(ITM[i][x])*re*dt*(a[k+1]+b[k+1]); 
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
        smMult(M, IAC[z][BETA], ITM, RSDIAG-1, RSDIAG-1, dimR, Nx/2, x0);
	memset(IFb[0], 0, dimR*(Nx/2)*sizeof(mcomplex));
	for (i = 0; i < dimR; ++i)
	  {
	    for (j = 0; j < qpts; ++j)
	      {
		for (x = x0; x < Nx/2; ++x)
		  {
		    Re(IFb[i][x]) += Rw[i][j] *
		      (-Kz[z]*Im(LIU[z][HXEL][j][x]) + Kx[x]*Im(LIU[z][HZEL][j][x]))+
		      Rpw[i][j]*(-Kz[z]*Im(LIU[z][DZEL][j][x])-Kx[x]*Im(LIU[z][DXEL][j][x]));
		    Im(IFb[i][x]) += Rw[i][j] *
		      (Kz[z]*Re(LIU[z][HXEL][j][x]) - Kx[x]*Re(LIU[z][HZEL][j][x]))+
		      Rpw[i][j]*(Kz[z]*Re(LIU[z][DZEL][j][x])+Kx[x]*Re(LIU[z][DXEL][j][x]));
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
		Re(IAC[z][BETA][i][x]) = Re(ITM[i][x]) + dt*d[k+1]*Re(IFb[i][x]);
		Im(IAC[z][BETA][i][x]) = Im(ITM[i][x]) + dt*d[k+1]*Im(IFb[i][x]);
	      }
	  }

	memset(M[0][0], 0, dimR*9*(Nx/2)*sizeof(double));
	for (i = 0; i < dimR; ++i)
        {
            for (j = 0; j < T_RSDIAG; ++j)
	    {
                for (x = x0; x < Nx/2; ++x)
	        {
                    M[i][j][x] = -(Rps[i][j] + K2[z][x]*Rs[i][j]);
	        }
            }
        }
	smMult(M, MIC[count-1][z][BETA], ITM, RSDIAG-1, RSDIAG-1, dimR, Nx/2, x0);
	for (i = 0; i < dimR; ++i)
	{
	  for (x=x0; x< Nx/2; x++)
	    {
	      Re(IAC[z][BETA][i][x]) += Re(ITM[i][x])*re*dt*(a[k+1]+b[k+1]);
	      Im(IAC[z][BETA][i][x]) += Im(ITM[i][x])*re*dt*(a[k+1]+b[k+1]); 
	    }
	}

    }
}    /* end project for nonzero (Kx,Kz) */

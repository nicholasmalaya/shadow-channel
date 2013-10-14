/***********************************************************
New piece of code to compute the boundary condition of incremental
state equation:
             IUx=dux*(W.n)
             IUy=0          for y=-1,
             IUz=duz*(W.n)
After we solve the state equation, dux hat and duz hat evaluated at
y=-1 is alread stored in Uxb and Uzb. Now we use FFT to compute IUx and IUz

***********************************************************************/

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include "fftw.h"
#include "rfftw.h"
#include "minChnl.h"

int increBoundary_s(void)
{
 /* External Variables */
    extern int Nx, Nz;
    extern fftw_complex ***CT;  /* 6-by-(3Nz/2)-by-(3*Nx/4+1) */
    extern mcomplex  **Uxb, **Uzb;
    extern fftw_plan pf1, pf2;
    extern rfftwnd_plan pr1, pr2;
    extern double *Kx, *Kz;

    int  x,i,z, idx;
    double norm, tmp1, tmp2,tmp3;
    fftw_real *RT;     /* real to complex transform */
    fftw_complex *fout = NULL;
    fftw_real *rout = NULL;

    idx = (3*Nz/2)*(3*Nx/2+2);
    RT = (fftw_real *)CT[0][0];
    norm = 1.0 / ((3.*Nx/2.)*(3.*Nz/2.));
 
    memset(CT[0][0], 0, MAXT*(3*Nz/2)*(3*Nx/4+1)*sizeof(fftw_complex));

      /* store Uxb hat and Uzb hat and w hat on CT for inverse FFT */
      for (z=0; z<Nz/2; ++z)
      {
	/* CT[0] store the data of Uxb, CT[1] storedata for Uzb */
	memcpy(CT[0][z], Uxb[z], (Nx/2)*sizeof(fftw_complex));
	memcpy(CT[1][z], Uzb[z], (Nx/2)*sizeof(fftw_complex));
      }

    for (z=Nz/2+1; z<Nz; ++z)
      {
	memcpy(CT[0][z+Nz/2], Uxb[z], (Nx/2)*sizeof(fftw_complex));
	memcpy(CT[1][z+Nz/2], Uzb[z], (Nx/2)*sizeof(fftw_complex));

      }

    // Re(CT[2][1][0])=1.;
    //  Re(CT[2][3*Nz/2-1][0])=1.;
    Re(CT[2][0][0])=1.;

    /* inverse Fourier transform */
	for (i = 0; i < 3; ++i)
	{
	   /* Each column of CT[i] */
	   fftw(pf1, Nx/2, CT[i][0], 3*Nx/4+1, 1, fout, -1, -1);

	   /* Each row of CT[i] */
	   rfftwnd_complex_to_real(pr1, 3*Nz/2, CT[i][0], 1, 3*Nx/4+1, rout, -1, -1);
	}

      	/* compute (dux)*(w.n) and (duz)*(w.n) */
	for (z = 0; z < (3*Nz/2); ++z)
	{
	    for (x = 0; x < 3*Nx/2; ++x)
	    {
	      RT[(z*(3*Nx/2+2)+x)] = RT[(z*(3*Nx/2+2)+x)]*RT[2*idx+(z*(3*Nx/2+2)+x)];
	      RT[idx+(z*(3*Nx/2+2)+x)] = RT[idx+(z*(3*Nx/2+2)+x)]*RT[2*idx+(z*(3*Nx/2+2)+x)];
	    }
	}


	/* Fourier transform to get Uxb hats and Uzb hats. */
	for (i = 0; i < 3; ++i)
	{

	    /* Each row of RT[i] */
	    rfftwnd_real_to_complex(pr2, 3*Nz/2, RT+(i*idx), 1, 3*Nx/2+2, 
		fout, -1, -1);

	    /* Each column of CT[i] */
	    fftw(pf2, Nx/2, CT[i][0], 3*Nx/4+1, 1, fout, -1, -1);

	    /* constant of FFT */
	    for (z = 0; z < Nz/2; ++z)
	    {
	        for (x = 0; x < Nx/2; ++x)
	        {
		    Re(CT[i][z][x]) = norm*Re(CT[i][z][x]);
		    Im(CT[i][z][x]) = norm*Im(CT[i][z][x]);
	        }
	    }

	    for (z = Nz+1; z < 3*Nz/2; ++z)
	    {
	        for (x = 0; x < Nx/2; ++x)
	        {
		    Re(CT[i][z][x]) = norm*Re(CT[i][z][x]);
		    Im(CT[i][z][x]) = norm*Im(CT[i][z][x]);
	        }
	    }
	}

	/*put date back in array Uxb and Uzb*/
	memset(Uxb[0], 0,  Nz*(Nx/2)*sizeof(mcomplex));
	memset(Uzb[0], 0, Nz*(Nx/2)*sizeof(mcomplex));
	for (z=0; z<Nz/2; ++z)
	  {
	    memcpy(Uxb[z], CT[0][z], Nx/2*sizeof(fftw_complex));
	    memcpy(Uzb[z], CT[1][z], Nx/2*sizeof(fftw_complex));
	  }
	for (z=Nz+1; z<3*Nz/2; ++z)
	  {
	    memcpy(Uxb[z-Nz/2], CT[0][z], Nx/2*sizeof(fftw_complex));
	    memcpy(Uzb[z-Nz/2], CT[1][z], Nx/2*sizeof(fftw_complex));
	  }

	/* further computation to get c1, c2, c3, c4 as in the note and results are rewritten in
	   Uxb and Uzb:
	   c1=g hat=iKz*Uxb-iKx*Uzb;-------rewrittin in Uxb
	   c2=du_y=-iKx*Uxb-iKz*Uzb;-------rewrittin in Uzb
	   c3=U=Uxb(0,0);           -------rewritten in Uxb[0][0]
	   c4=Uzb(0,0)              -------rewritten in Uzb[0][0]*/

	for (z = 0; z < Nz; ++z)
	    {
	        for (x = 0; x < Nx/2; ++x)
	        {
		  if(z*z+x*x>0)
		    {
		        tmp1=-Kz[z]*Im(Uxb[z][x])+Kx[x]*Im(Uzb[z][x]);
			tmp2= Kz[z]*Re(Uxb[z][x])-Kx[x]*Re(Uzb[z][x]);
			tmp3= Kx[x]*Im(Uxb[z][x])+Kz[z]*Im(Uzb[z][x]);

			Im(Uzb[z][x])=-Kx[x]*Re(Uxb[z][x])-Kz[z]*Re(Uzb[z][x]);
			Re(Uzb[z][x])=tmp3;
			Re(Uxb[z][x])=tmp1;
			Im(Uxb[z][x])=tmp2;
		    }
		}
	    }


	return(NO_ERR);
}

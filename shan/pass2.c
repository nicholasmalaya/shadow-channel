/*************************************
function used to compute the nonlinear term for adjoint and incremental adjoint systems. 
Now at the end of the code, H is stored in two parts in HXEL, HYEL and HZEl. and another part 
which need to be take care of in assembling the Fs, ( y derivatives) and they are stored in DXEL
and DZEL.

*****************************************************************************************/



/* Compute H =\nabla x (u x lambda ) + lambda x (\nabla x u)  for adjoint
           IH= \nabla x (u x Ilambda) +Ilambda x (\nabla x u)
	    + \nabla x (Iu x lambda)+  lambda x (\nabla x Iu) for incremental adjoint */

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include "fftw.h"
#include "rfftw.h"
#include "minChnl.h"

int pass2(int dctr, int n)
{
    /* External Variables */
    extern int qpts, Nx, Nz;
    extern double dt, cfl1, cfl3, *cfl2;
    extern double *Kx, *Kz, *Qy;
    extern fftw_complex ***CT, ***ICT;  /* 9-by-(3Nz/2)-by-(3*Nx/4+1) */
    extern mcomplex ****U, ****IU;
    extern mcomplex ****AU, ****IAU;
    extern fftw_plan pf1, pf2;
    extern rfftwnd_plan pr1, pr2;

    extern fftw_plan Ipf1, Ipf2;

    /* Local variables */
    int x, y, z, i, j, idx;
    double uy, norm, t, CFL, Iuy, wy, Iwy;
    fftw_real *RT, *IRT;     /* real to complex transform */
    fftw_complex *fout = NULL;
    fftw_real *rout = NULL;

    idx = (3*Nz/2)*(3*Nx/2+2);
     RT = (fftw_real *)CT[0][0];
    IRT = (fftw_real *)ICT[0][0];
    norm = 1.0 / ((3.*Nx/2.)*(3.*Nz/2.));
    static double h[3]={ 0., 1., 2./3.};
    for (y = 0; y < qpts; ++y)
    {

      /* first compute omega= \nabla x u */

	/* Clear the array CT */

      memset(CT[0][0], 0, MAXTT*(3*Nz/2)*(3*Nx/4+1)*sizeof(fftw_complex));

	/* Load "first half" of array */
	for (i = 0; i < MAXU; ++i)
	{
	    for (z = 0; z < Nz/2; ++z)
	    {
		memcpy(CT[i][z], U[z][i][y], (Nx/2)*sizeof(fftw_complex));
	    }
	}
	for (z = 0; z < Nz/2; ++z)
	{
	    for (x = 0; x < (Nx/2); ++x)
	    {
	        /* omega_x hat */
		Re(CT[WXEL][z][x]) = Re(CT[DZEL][z][x]) +
		   Kz[z]*Im(CT[YEL][z][x]);
		Im(CT[WXEL][z][x]) = Im(CT[DZEL][z][x]) -
		   Kz[z]*Re(CT[YEL][z][x]);
	        /* omega_z hat */
		Re(CT[WZEL][z][x]) = -(Re(CT[DXEL][z][x]) +
		   Kx[x]*Im(CT[YEL][z][x]));
		Im(CT[WZEL][z][x]) = -Im(CT[DXEL][z][x]) +
		   Kx[x]*Re(CT[YEL][z][x]);
	        /* omega_y hat */
		Re(CT[WYEL][z][x]) = Kx[x]*Im(CT[ZEL][z][x]) -
		   Kz[z]*Im(CT[XEL][z][x]);
		Im(CT[WYEL][z][x]) = Kz[z]*Re(CT[XEL][z][x]) -
		   Kx[x]*Re(CT[ZEL][z][x]);

	    }
	}

	/* Load "second half" of array */
	for (i = 0; i < MAXU; ++i)
	{
	    for (z = Nz/2+1; z < Nz; ++z)
	    {
		memcpy(CT[i][z+Nz/2], U[z][i][y], (Nx/2)*sizeof(fftw_complex));
	    }
	}
	for (z = Nz+1; z < (3*Nz/2); ++z)
	{
	    for (x = 0; x < (Nx/2); ++x)
	    {
	        /* omega_x hat */
		Re(CT[WXEL][z][x]) = Re(CT[DZEL][z][x]) + 
		   Kz[z-Nz/2]*Im(CT[YEL][z][x]);
		Im(CT[WXEL][z][x]) = Im(CT[DZEL][z][x]) - 
		   Kz[z-Nz/2]*Re(CT[YEL][z][x]);
	        /* omega_z hat */
		Re(CT[WZEL][z][x]) = -(Re(CT[DXEL][z][x]) + 
		   Kx[x]*Im(CT[YEL][z][x]));
		Im(CT[WZEL][z][x]) = -Im(CT[DXEL][z][x]) +
		   Kx[x]*Re(CT[YEL][z][x]);
	        /* omega_y hat */
		Re(CT[WYEL][z][x]) = Kx[x]*Im(CT[ZEL][z][x]) -   
		   Kz[z-Nz/2]*Im(CT[XEL][z][x]);
		Im(CT[WYEL][z][x]) = Kz[z-Nz/2]*Re(CT[XEL][z][x]) -
		   Kx[x]*Re(CT[ZEL][z][x]);
	    }
	}

	/*store AU_x, AU_y, AU_z to the last three components of CT */
	for (i = 0; i < 3; ++i)
	{
	    for (z = 0; z < Nz/2; ++z)
	    {
		memcpy(CT[i+MAXT][z], AU[z][i][y], (Nx/2)*sizeof(fftw_complex));
	    }
	}
	for (i = 0; i < 3; ++i)
	{
	    for (z = Nz/2+1; z < Nz; ++z)
	    {
		memcpy(CT[i+MAXT][z+Nz/2], AU[z][i][y], (Nx/2)*sizeof(fftw_complex));
	    }
	}

	/* inverse Fourier transform */
	for (i = 0; i < MAXTT; ++i)
	{
	   /* Each column of CT[i] */
	   fftw(pf1, Nx/2, CT[i][0], 3*Nx/4+1, 1, fout, -1, -1);

	   /* Each row of CT[i] */
	   rfftwnd_complex_to_real(pr1, 3*Nz/2, CT[i][0], 1, 3*Nx/4+1, 
	       rout, -1, -1);
	}
	
        /* now compute Iomega= \nabla x Iu */

	/* Clear the array ICT */
	memset(ICT[0][0], 0, MAXTT*(3*Nz/2)*(3*Nx/4+1)*sizeof(fftw_complex));

	/* Load "first half" of array */
	for (i = 0; i < MAXU; ++i)
	{
	    for (z = 0; z < Nz/2; ++z)
	    {
		memcpy(ICT[i][z], IU[z][i][y], (Nx/2)*sizeof(fftw_complex));
	    }
	}
	for (z = 0; z < Nz/2; ++z)
	{
	    for (x = 0; x < (Nx/2); ++x)
	    {
	        /* Iomega_x hat */
		Re(ICT[WXEL][z][x]) = Re(ICT[DZEL][z][x]) + 
		   Kz[z]*Im(ICT[YEL][z][x]);
		Im(ICT[WXEL][z][x]) = Im(ICT[DZEL][z][x]) - 
		   Kz[z]*Re(ICT[YEL][z][x]);
	        /* Iomega_z hat */
		Re(ICT[WZEL][z][x]) = -(Re(ICT[DXEL][z][x]) + 
		   Kx[x]*Im(ICT[YEL][z][x]));
		Im(ICT[WZEL][z][x]) = -Im(ICT[DXEL][z][x]) +
		   Kx[x]*Re(ICT[YEL][z][x]);
	        /* Iomega_y hat */
		Re(ICT[WYEL][z][x]) = Kx[x]*Im(ICT[ZEL][z][x]) -   
		   Kz[z]*Im(ICT[XEL][z][x]);
		Im(ICT[WYEL][z][x]) = Kz[z]*Re(ICT[XEL][z][x]) -
		   Kx[x]*Re(ICT[ZEL][z][x]);
	    }
	}

	/* Load "second half" of array */
	for (i = 0; i < MAXU; ++i)
	{
	    for (z = Nz/2+1; z < Nz; ++z)
	    {
		memcpy(ICT[i][z+Nz/2], IU[z][i][y], (Nx/2)*sizeof(fftw_complex));
	    }
	}
	for (z = Nz+1; z < (3*Nz/2); ++z)
	{
	    for (x = 0; x < (Nx/2); ++x)
	    {
	        /* omega_x hat */
		Re(ICT[WXEL][z][x]) = Re(ICT[DZEL][z][x]) + 
		   Kz[z-Nz/2]*Im(ICT[YEL][z][x]);
		Im(ICT[WXEL][z][x]) = Im(ICT[DZEL][z][x]) - 
		   Kz[z-Nz/2]*Re(ICT[YEL][z][x]);
	        /* omega_z hat */
		Re(ICT[WZEL][z][x]) = -(Re(ICT[DXEL][z][x]) + 
		   Kx[x]*Im(ICT[YEL][z][x]));
		Im(ICT[WZEL][z][x]) = -Im(ICT[DXEL][z][x]) +
		   Kx[x]*Re(ICT[YEL][z][x]);
	        /* omega_y hat */
		Re(ICT[WYEL][z][x]) = Kx[x]*Im(ICT[ZEL][z][x]) -   
		   Kz[z-Nz/2]*Im(ICT[XEL][z][x]);
		Im(ICT[WYEL][z][x]) = Kz[z-Nz/2]*Re(ICT[XEL][z][x]) -
		   Kx[x]*Re(ICT[ZEL][z][x]);
	    }
	}


	/*store IAU_x, IAU_y, IAU_z to the last three components of ICT */
	for (i = 0; i < 3; ++i)
	{
	    for (z = 0; z < Nz/2; ++z)
	    {
		memcpy(ICT[i+MAXT][z], IAU[z][i][y], (Nx/2)*sizeof(fftw_complex));
	    }
	}
	for (i = 0; i < 3; ++i)
	{
	    for (z = Nz/2+1; z < Nz; ++z)
	      {
		memcpy(ICT[i+MAXT][z+Nz/2], IAU[z][i][y], (Nx/2)*sizeof(fftw_complex));
	      }
	}

	/* inverse Fourier transform */
	for (i = 0; i < MAXTT; ++i)
	{
	   /* Each column of ICT[i] */
	   fftw(Ipf1, Nx/2, ICT[i][0], 3*Nx/4+1, 1, fout, -1, -1);

	   /* Each row of ICT[i] */
	   rfftwnd_complex_to_real(pr1, 3*Nz/2, ICT[i][0], 1, 3*Nx/4+1, 
	       rout, -1, -1);
	}


	for (z = 0; z < (3*Nz/2); ++z)
	{
	    for (x = 0; x < 3*Nx/2; ++x)
	    {
                Iuy = IRT[YEL*idx+(z*(3*Nx/2+2)+x)];
         	uy = RT[YEL*idx+(z*(3*Nx/2+2)+x)];

		wy=RT[WYEL*idx+(z*(3*Nx/2+2)+x)];
		Iwy=IRT[WYEL*idx+(z*(3*Nx/2+2)+x)];

		//compute  IS= (u x Ilambda +Iu x lambda ) and store result in SXEL, SYEL, SZEL,
		IRT[SXEL*idx+(z*(3*Nx/2+2)+x)] =
		  RT[YEL*idx+(z*(3*Nx/2+2)+x)]*IRT[AZEL*idx+(z*(3*Nx/2+2)+x)] -
		  RT[ZEL*idx+(z*(3*Nx/2+2)+x)]*IRT[AYEL*idx+(z*(3*Nx/2+2)+x)]+
		  IRT[YEL*idx+(z*(3*Nx/2+2)+x)]*RT[AZEL*idx+(z*(3*Nx/2+2)+x)] -
		  IRT[ZEL*idx+(z*(3*Nx/2+2)+x)]*RT[AYEL*idx+(z*(3*Nx/2+2)+x)];

		IRT[SYEL*idx+(z*(3*Nx/2+2)+x)] =
		   RT[ZEL*idx+(z*(3*Nx/2+2)+x)]*IRT[AXEL*idx+(z*(3*Nx/2+2)+x)] -
		   RT[XEL*idx+(z*(3*Nx/2+2)+x)]*IRT[AZEL*idx+(z*(3*Nx/2+2)+x)]+
		   IRT[ZEL*idx+(z*(3*Nx/2+2)+x)]*RT[AXEL*idx+(z*(3*Nx/2+2)+x)] -
		   IRT[XEL*idx+(z*(3*Nx/2+2)+x)]*RT[AZEL*idx+(z*(3*Nx/2+2)+x)];

		IRT[SZEL*idx+(z*(3*Nx/2+2)+x)] =
		   RT[XEL*idx+(z*(3*Nx/2+2)+x)]*IRT[AYEL*idx+(z*(3*Nx/2+2)+x)] -
		   uy*IRT[AXEL*idx+(z*(3*Nx/2+2)+x)]+
		  IRT[XEL*idx+(z*(3*Nx/2+2)+x)]*RT[AYEL*idx+(z*(3*Nx/2+2)+x)] -
		   Iuy*RT[AXEL*idx+(z*(3*Nx/2+2)+x)];

;
		// compute IT= lambda x Iomega+ Ilambda x omega and store results in TXEL, TYEL, TZEL
		IRT[TXEL*idx+(z*(3*Nx/2+2)+x)] =
	           RT[AYEL*idx+(z*(3*Nx/2+2)+x)]*IRT[WZEL*idx+(z*(3*Nx/2+2)+x)] -
		   RT[AZEL*idx+(z*(3*Nx/2+2)+x)]*IRT[WYEL*idx+(z*(3*Nx/2+2)+x)] +
		   IRT[AYEL*idx+(z*(3*Nx/2+2)+x)]*RT[WZEL*idx+(z*(3*Nx/2+2)+x)] -
		   IRT[AZEL*idx+(z*(3*Nx/2+2)+x)]*RT[WYEL*idx+(z*(3*Nx/2+2)+x)];

		IRT[TYEL*idx+(z*(3*Nx/2+2)+x)] =
		   RT[AZEL*idx+(z*(3*Nx/2+2)+x)]*IRT[WXEL*idx+(z*(3*Nx/2+2)+x)] -
		   RT[AXEL*idx+(z*(3*Nx/2+2)+x)]*IRT[WZEL*idx+(z*(3*Nx/2+2)+x)]+
		  IRT[AZEL*idx+(z*(3*Nx/2+2)+x)]*RT[WXEL*idx+(z*(3*Nx/2+2)+x)] -
		  IRT[AXEL*idx+(z*(3*Nx/2+2)+x)]*RT[WZEL*idx+(z*(3*Nx/2+2)+x)];

		IRT[TZEL*idx+(z*(3*Nx/2+2)+x)] =
		   RT[AXEL*idx+(z*(3*Nx/2+2)+x)]*Iwy-
		  RT[AYEL*idx+(z*(3*Nx/2+2)+x)]*IRT[WXEL*idx+(z*(3*Nx/2+2)+x)]+
		  IRT[AXEL*idx+(z*(3*Nx/2+2)+x)]*wy-
		  IRT[AYEL*idx+(z*(3*Nx/2+2)+x)]*RT[WXEL*idx+(z*(3*Nx/2+2)+x)];

		// compute S= (u x lambda ) and store result in SXEL, SYEL, SZEL,
		RT[SXEL*idx+(z*(3*Nx/2+2)+x)] =
	           RT[YEL*idx+(z*(3*Nx/2+2)+x)]*RT[AZEL*idx+(z*(3*Nx/2+2)+x)] -
		   RT[ZEL*idx+(z*(3*Nx/2+2)+x)]*RT[AYEL*idx+(z*(3*Nx/2+2)+x)];
		RT[SYEL*idx+(z*(3*Nx/2+2)+x)] =
		   RT[ZEL*idx+(z*(3*Nx/2+2)+x)]*RT[AXEL*idx+(z*(3*Nx/2+2)+x)] -
		   RT[XEL*idx+(z*(3*Nx/2+2)+x)]*RT[AZEL*idx+(z*(3*Nx/2+2)+x)];
		RT[SZEL*idx+(z*(3*Nx/2+2)+x)] =
		   RT[XEL*idx+(z*(3*Nx/2+2)+x)]*RT[AYEL*idx+(z*(3*Nx/2+2)+x)] -
		   uy*RT[AXEL*idx+(z*(3*Nx/2+2)+x)];

		// compute T= lambda x omega and store results in TXEL, TYEL, TZEL
		RT[TXEL*idx+(z*(3*Nx/2+2)+x)] =
	           RT[AYEL*idx+(z*(3*Nx/2+2)+x)]*RT[WZEL*idx+(z*(3*Nx/2+2)+x)] -
		   RT[AZEL*idx+(z*(3*Nx/2+2)+x)]*RT[WYEL*idx+(z*(3*Nx/2+2)+x)];
		RT[TYEL*idx+(z*(3*Nx/2+2)+x)] =
		   RT[AZEL*idx+(z*(3*Nx/2+2)+x)]*RT[WXEL*idx+(z*(3*Nx/2+2)+x)] -
		   RT[AXEL*idx+(z*(3*Nx/2+2)+x)]*RT[WZEL*idx+(z*(3*Nx/2+2)+x)];
		RT[TZEL*idx+(z*(3*Nx/2+2)+x)] =
		   RT[AXEL*idx+(z*(3*Nx/2+2)+x)]*wy-
		  RT[AYEL*idx+(z*(3*Nx/2+2)+x)]*RT[WXEL*idx+(z*(3*Nx/2+2)+x)];

	    }
	}

	/* Fourier transform to get H hats and IH hats. */
	for (i = 0; i < 6; ++i)
	{
	  /* start of H hats */
	    /* Each row of RT[i] */
	    rfftwnd_real_to_complex(pr2, 3*Nz/2, RT+(i*idx), 1, 3*Nx/2+2, 
		fout, -1, -1);

	    /* Each column of CT[i] */
	    fftw(pf2, Nx/2, CT[i][0], 3*Nx/4+1, 1, fout, -1, -1);

	    /* Put data back in array U */
	    for (z = 0; z < Nz/2; ++z)
	    {
	        for (x = 0; x < Nx/2; ++x)
	        {
		    Re(CT[i][z][x]) = norm*Re(CT[i][z][x]);
		    Im(CT[i][z][x]) = norm*Im(CT[i][z][x]);
	        }

		//memcpy(U[z][i][y], CT[i][z], Nx/2*sizeof(fftw_complex));
	    }
	    for (z = Nz+1; z < 3*Nz/2; ++z)
	    {
	        for (x = 0; x < Nx/2; ++x)
	        {
		    Re(CT[i][z][x]) = norm*Re(CT[i][z][x]);
		    Im(CT[i][z][x]) = norm*Im(CT[i][z][x]);
	        }

		//	memcpy(U[z - Nz/2][i][y], CT[i][z], Nx/2*sizeof(fftw_complex));
	    }
            /* end of H hats */

	 /* start of IH hats */

	    /* Each row of RT[i] */
	    rfftwnd_real_to_complex(pr2, 3*Nz/2, IRT+(i*idx), 1, 3*Nx/2+2, 
		fout, -1, -1);

	    /* Each column of CT[i] */
	    fftw(Ipf2, Nx/2, ICT[i][0], 3*Nx/4+1, 1, fout, -1, -1);

	    /* Put data back in array U */
	    for (z = 0; z < Nz/2; ++z)
	    {
	        for (x = 0; x < Nx/2; ++x)
	        {
		    Re(ICT[i][z][x]) = norm*Re(ICT[i][z][x]);
		    Im(ICT[i][z][x]) = norm*Im(ICT[i][z][x]);
	        }

		//	memcpy(IU[z][i][y], ICT[i][z], Nx/2*sizeof(fftw_complex));
	    }
	    for (z = Nz+1; z < 3*Nz/2; ++z)
	    {
	        for (x = 0; x < Nx/2; ++x)
	        {
		    Re(ICT[i][z][x]) = norm*Re(ICT[i][z][x]);
		    Im(ICT[i][z][x]) = norm*Im(ICT[i][z][x]);
	        }

		//memcpy(IU[z - Nz/2][i][y], ICT[i][z], Nx/2*sizeof(fftw_complex));
	    }
	  /*End of IH hats */
	}

	 for (z = 0; z < Nz/2; ++z)
	    {
	      memcpy(U[z][DXEL][y], CT[SXEL][z], Nx/2*sizeof(fftw_complex));
	      memcpy(U[z][DZEL][y], CT[SZEL][z], Nx/2*sizeof(fftw_complex));
	      memcpy(IU[z][DXEL][y], ICT[SXEL][z], Nx/2*sizeof(fftw_complex));
	      memcpy(IU[z][DZEL][y], ICT[SZEL][z], Nx/2*sizeof(fftw_complex));
	      for (x = 0; x < Nx/2; ++x)
		{
		  Re(U[z][HXEL][y][x])=-Re(CT[TXEL][z][x])- Kz[z]*Im(CT[SYEL][z][x]);
		  Im(U[z][HXEL][y][x])=-Im(CT[TXEL][z][x])+Kz[z]*Re(CT[SYEL][z][x]);
		  Re(U[z][HYEL][y][x])=-Re(CT[TYEL][z][x])+Kz[z]*Im(CT[SXEL][z][x])- Kx[x]*Im(CT[SZEL][z][x]);
		  Im(U[z][HYEL][y][x])=-Im(CT[TYEL][z][x])-Kz[z]*Re(CT[SXEL][z][x])+Kx[x]*Re(CT[SZEL][z][x]);
		  Re(U[z][HZEL][y][x])=-Re(CT[TZEL][z][x])+ Kx[x]*Im(CT[SYEL][z][x]);
		  Im(U[z][HZEL][y][x])=-Im(CT[TZEL][z][x])- Kx[x]*Re(CT[SYEL][z][x]);

		  Re(IU[z][HXEL][y][x])=-Re(ICT[TXEL][z][x])- Kz[z]*Im(ICT[SYEL][z][x]);
		  Im(IU[z][HXEL][y][x])=-Im(ICT[TXEL][z][x])+Kz[z]*Re(ICT[SYEL][z][x]);
		  Re(IU[z][HYEL][y][x])=-Re(ICT[TYEL][z][x])+Kz[z]*Im(ICT[SXEL][z][x])- Kx[x]*Im(ICT[SZEL][z][x]);
		  Im(IU[z][HYEL][y][x])=-Im(ICT[TYEL][z][x])-Kz[z]*Re(ICT[SXEL][z][x])+Kx[x]*Re(ICT[SZEL][z][x]);
		  Re(IU[z][HZEL][y][x])=-Re(ICT[TZEL][z][x])+ Kx[x]*Im(ICT[SYEL][z][x]);
		  Im(IU[z][HZEL][y][x])=-Im(ICT[TZEL][z][x])- Kx[x]*Re(ICT[SYEL][z][x]);
		}
	    }

	 for (z = Nz+1; z <3* Nz/2; ++z)
	    {
	      memcpy(U[z-Nz/2][DXEL][y], CT[SXEL][z], Nx/2*sizeof(fftw_complex));
	      memcpy(U[z-Nz/2][DZEL][y], CT[SZEL][z], Nx/2*sizeof(fftw_complex));
	      memcpy(IU[z-Nz/2][DXEL][y], ICT[SXEL][z], Nx/2*sizeof(fftw_complex));
	      memcpy(IU[z-Nz/2][DZEL][y], ICT[SZEL][z], Nx/2*sizeof(fftw_complex));
	      for (x = 0; x < Nx/2; ++x)
		{
		  Re(U[z-Nz/2][HXEL][y][x])=-Re(CT[TXEL][z][x])- Kz[z-Nz/2]*Im(CT[SYEL][z][x]);
		  Im(U[z-Nz/2][HXEL][y][x])=-Im(CT[TXEL][z][x])+Kz[z-Nz/2]*Re(CT[SYEL][z][x]);
		  Re(U[z-Nz/2][HYEL][y][x])=-Re(CT[TYEL][z][x])+Kz[z-Nz/2]*Im(CT[SXEL][z][x])- Kx[x]*Im(CT[SZEL][z][x]);
		  Im(U[z-Nz/2][HYEL][y][x])=-Im(CT[TYEL][z][x])-Kz[z-Nz/2]*Re(CT[SXEL][z][x])+Kx[x]*Re(CT[SZEL][z][x]);
		  Re(U[z-Nz/2][HZEL][y][x])=-Re(CT[TZEL][z][x])+ Kx[x]*Im(CT[SYEL][z][x]);
		  Im(U[z-Nz/2][HZEL][y][x])=-Im(CT[TZEL][z][x])- Kx[x]*Re(CT[SYEL][z][x]);

		  Re(IU[z-Nz/2][HXEL][y][x])=-Re(ICT[TXEL][z][x])- Kz[z-Nz/2]*Im(ICT[SYEL][z][x]);
		  Im(IU[z-Nz/2][HXEL][y][x])=-Im(ICT[TXEL][z][x])+Kz[z-Nz/2]*Re(ICT[SYEL][z][x]);
		  Re(IU[z-Nz/2][HYEL][y][x])=-Re(ICT[TYEL][z][x])+Kz[z-Nz/2]*Im(ICT[SXEL][z][x])- Kx[x]*Im(ICT[SZEL][z][x]);
		  Im(IU[z-Nz/2][HYEL][y][x])=-Im(ICT[TYEL][z][x])-Kz[z-Nz/2]*Re(ICT[SXEL][z][x])+Kx[x]*Re(ICT[SZEL][z][x]);
		  Re(IU[z-Nz/2][HZEL][y][x])=-Re(ICT[TZEL][z][x])+ Kx[x]*Im(ICT[SYEL][z][x]);
		  Im(IU[z-Nz/2][HZEL][y][x])=-Im(ICT[TZEL][z][x])- Kx[x]*Re(ICT[SYEL][z][x]);
		}
	    }
    }   /* end for y... */

    return(NO_ERR);
}    /* end pass2 */

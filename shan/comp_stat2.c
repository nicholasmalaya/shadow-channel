#include <math.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include "fftw.h"
#include "rfftw.h"
#include "minChnl.h"

int comp_stat2(int n)
{
 /* External Variables */
  extern int Nx, Nz, qpts, nums;
  extern mcomplex ****AU,****IAU, **GIUxb, ishear, **GIUzb, **GUxb, **GUzb;
  extern double **uu, **us, *W;
    extern FILE *afp, *afp2, *afp3, *afp4, *afp5, *afp6, *afp7, *afp8, *afp9, *afp10;
    extern fftw_complex ***CT;
    extern fftw_plan pf1, pf2;
    extern fftw_plan Ipf1, Ipf2;
    extern rfftwnd_plan pr1, pr2;
    extern double *Kx, *Kz, **K2;
    extern double re, dt, atau, iatau;
    int  x,i,j,y, z, idx;
    fftw_real *RT;     /* real to complex transform */
    fftw_complex *fout = NULL;
    fftw_real *rout = NULL;
    extern mcomplex ****U,****IU;
    double norm; 
   

    idx = (3*Nz/2)*(3*Nx/2+2);
    RT = (fftw_real *)CT[0][0];
    norm = 1.0 / ((3.*Nx/2.)*(3.*Nz/2.));
    double u1, u2, u3, u1u2, u1u3, u2u3;
    double iu1, iu2, iu3, iu1u2, iu1u3, iu2u3, u1u1, u2u2, u3u3, iu1u1, iu2u2, iu3u3;
    double au1, au2, au3, au1u2, au1u3, au2u3, au1u1, au2u2, au3u3;
    double iau1, iau2, iau3, iau1u2, iau1u3, iau2u3, iau1u1, iau2u2, iau3u3;

    atau+=Re(GIUxb[0][0]);
    iatau+=Re(ishear);

    if((n-1) % 100==0)
      {
    /* record wall shear stress for adjoint and incremental adjoint */
	fprintf(afp7, "%f\n", atau/100.);
	atau=0;
	fprintf(afp8, "%f\n", iatau/100.);
	iatau=0;

	au1=0; au2=0; au3=0; au1u2=0; au2u3=0; au1u3=0; au1u1=0; au2u2=0; au3u3=0;
	iau1=0; iau2=0; iau3=0; iau1u2=0; iau2u3=0; iau1u3=0; iau1u1=0; iau2u2=0; iau3u3=0;

	nums=nums+1;
		/*fprintf(fp, "the time averge %d, %d\n", (n+1)-50, (n+1));
	for(j=0; j< qpts; j++)
	  {
	    fprintf(fp, "%d, %f %f %f  \n", j, uu[0][j]/50., uu[1][j]/50., uu[2][j]/50.);
	  }

	fprintf(fp2, "the time averge %d, %d\n", (n+1)-50, (n+1));
	for(j=0; j< qpts; j++)
	  {
	    fprintf(fp2, "%d, %f %f %f  \n", j, uu[3][j]/50., uu[4][j]/50., uu[5][j]/50.);
	    }

	    memset(uu[0], 0, 6*(qpts)*sizeof(double));*/


	memset(CT[0][0], 0,  MAXTT*(3*Nz/2)*(3*Nx/4+1)*sizeof(fftw_complex));
	for (z = 0; z < Nz/2; ++z)
	  {
	    memcpy(CT[0][z], GIUxb[z], (Nx/2)*sizeof(fftw_complex));
	    memcpy(CT[1][z], GIUzb[z], (Nx/2)*sizeof(fftw_complex));
	    memcpy(CT[2][z], GUxb[z], (Nx/2)*sizeof(fftw_complex));
	    memcpy(CT[3][z], GUzb[z], (Nx/2)*sizeof(fftw_complex));
	  }
	for (z = Nz/2+1; z < Nz; ++z)
	  {
	    memcpy(CT[0][z+Nz/2], GIUxb[z], (Nx/2)*sizeof(fftw_complex));
	    memcpy(CT[1][z+Nz/2], GIUzb[z], (Nx/2)*sizeof(fftw_complex));
	    memcpy(CT[2][z+Nz/2], GUxb[z], (Nx/2)*sizeof(fftw_complex));
	    memcpy(CT[3][z+Nz/2], GUzb[z], (Nx/2)*sizeof(fftw_complex));
	  }

	for (i = 0; i < 4; ++i)
	      {
		/* Each column of CT[i] */
		fftw(pf1, Nx/2, CT[i][0], 3*Nx/4+1, 1, fout, -1, -1);

		/* Each row of CT[i] */
		rfftwnd_complex_to_real(pr1, 3*Nz/2, CT[i][0], 1, 3*Nx/4+1, 
					rout, -1, -1);
	      }
	z=5;
	for(x=0; x< Nx/2; x++)
	  {
	    fprintf(afp, " %f, %f,  %f \n",  RT[(z*(3*Nx/2+2)+x)], RT[idx+(z*(3*Nx/2+2)+x)],
		    RT[(z*(3*Nx/2+2)+x)]*RT[idx*2+(z*(3*Nx/2+2)+x)]
		    + RT[idx+(z*(3*Nx/2+2)+x)]*RT[idx*3+(z*(3*Nx/2+2)+x)]);
	  }

	z=20;
	for(x=0; x< Nx/2; x++)
	  {
	    fprintf(afp2, " %f, %f,  %f \n",RT[(z*(3*Nx/2+2)+x)], RT[idx+(z*(3*Nx/2+2)+x)],
		    RT[(z*(3*Nx/2+2)+x)]*RT[idx*2+(z*(3*Nx/2+2)+x)]
		    + RT[idx+(z*(3*Nx/2+2)+x)]*RT[idx*3+(z*(3*Nx/2+2)+x)]);
	  }

	for (y = 0; y < qpts; ++y)
	  {
	    memset(CT[0][0], 0, MAXTT*(3*Nz/2)*(3*Nx/4+1)*sizeof(fftw_complex));

	    for (i = 0; i < 3; ++i)
	      {
		for (z = 0; z < Nz/2; ++z)
		  {
		    memcpy(CT[i][z], AU[z][i][y], (Nx/2)*sizeof(fftw_complex));
		    memcpy(CT[i+3][z], IAU[z][i][y], (Nx/2)*sizeof(fftw_complex));
		  }
		for (z = Nz/2+1; z < Nz; ++z)
		  {
		    memcpy(CT[i][z+Nz/2], AU[z][i][y], (Nx/2)*sizeof(fftw_complex));
		    memcpy(CT[i+3][z+Nz/2], IAU[z][i][y], (Nx/2)*sizeof(fftw_complex));
		  }
	      }

	    for (i = 0; i < 6; ++i)
	      {
		/* Each column of CT[i] */
		fftw(pf1, Nx/2, CT[i][0], 3*Nx/4+1, 1, fout, -1, -1);

		/* Each row of CT[i] */
		rfftwnd_complex_to_real(pr1, 3*Nz/2, CT[i][0], 1, 3*Nx/4+1, 
					rout, -1, -1);
	      }

	    u1=0;
	    u2=0;
	    u3=0;
	    u1u2=0;
	    u1u3=0;
	    u2u3=0;
	    u1u1=0;
	    u2u2=0;
	    u3u3=0;
	    for (z = 0; z < (3*Nz/2); ++z)
	      {
		for (x = 0; x < 3*Nx/2; ++x)
		  {
		    u1 += RT[XEL*idx+(z*(3*Nx/2+2)+x)];
		    u2 += RT[YEL*idx+(z*(3*Nx/2+2)+x)];
		    u3 += RT[ZEL*idx+(z*(3*Nx/2+2)+x)];
		    u1u2 += RT[XEL*idx+(z*(3*Nx/2+2)+x)]* RT[YEL*idx+(z*(3*Nx/2+2)+x)];
		    u1u3 += RT[XEL*idx+(z*(3*Nx/2+2)+x)]* RT[ZEL*idx+(z*(3*Nx/2+2)+x)];
		    u2u3 += RT[YEL*idx+(z*(3*Nx/2+2)+x)]* RT[ZEL*idx+(z*(3*Nx/2+2)+x)];
		    u1u1 += RT[XEL*idx+(z*(3*Nx/2+2)+x)]* RT[XEL*idx+(z*(3*Nx/2+2)+x)];
		    u2u2 += RT[YEL*idx+(z*(3*Nx/2+2)+x)]* RT[YEL*idx+(z*(3*Nx/2+2)+x)];
		    u3u3 += RT[ZEL*idx+(z*(3*Nx/2+2)+x)]* RT[ZEL*idx+(z*(3*Nx/2+2)+x)];
		  }
	      }

	    u1=u1/(3*Nz/2)/(3*Nx/2);
	    u2=u2/(3*Nz/2)/(3*Nx/2);
	    u3=u3/(3*Nz/2)/(3*Nx/2);
	    u1u2=u1u2/(3*Nz/2)/(3*Nx/2);
	    u1u3=u1u3/(3*Nz/2)/(3*Nx/2);
	    u2u3=u2u3/(3*Nz/2)/(3*Nx/2);
	    u1u1=u1u1/(3*Nz/2)/(3*Nx/2);
	    u2u2=u2u2/(3*Nz/2)/(3*Nx/2);
	    u3u3=u3u3/(3*Nz/2)/(3*Nx/2);

	    us[0][y]+=u1;
	    us[1][y]+=u2;
	    us[2][y]+=u3;
	    us[3][y]+=u1u2;
	    us[4][y]+=u1u3;
	    us[5][y]+=u2u3;
	    us[12][y]+=u1u1;
	    us[13][y]+=u2u2;
	    us[14][y]+=u3u3;
	    // fprintf(fp, "%d, %f %f %f  %f %f  %f %f %f %f\n", y, u1, u2, u3, u1u2, u1u3, u2u3, u1u1, u2u2, u3u3);

	     if(y==5)
	      {
		fprintf(afp3, " %f %f %f  %f %f  %f  %f  %f  %f %f\n", u1, u2, u3, u1u2, u1u3, u2u3, u1u1, u2u2, u3u3, 0.5*(u1u1+u2u2+u3u3));
	      }

	    if(y==40)
	      {
		fprintf(afp4, " %f %f %f  %f %f  %f  %f  %f  %f %f \n", u1, u2, u3, u1u2, u1u3, u2u3, u1u1, u2u2, u3u3, 0.5*(u1u1+u2u2+u3u3));
		}

	    iu1=0;
	    iu2=0;
	    iu3=0;
	    iu1u2=0;
	    iu1u3=0;
	    iu2u3=0;
	    iu1u1=0;
	    iu2u2=0;
	    iu3u3=0;
	    for (z = 0; z < (3*Nz/2); ++z)
	      {
		for (x = 0; x < 3*Nx/2; ++x)
		  {
		    iu1 += RT[(XEL+3)*idx+(z*(3*Nx/2+2)+x)];
		    iu2 += RT[(YEL+3)*idx+(z*(3*Nx/2+2)+x)];
		    iu3 += RT[(ZEL+3)*idx+(z*(3*Nx/2+2)+x)];
		    iu1u2 += RT[(XEL+3)*idx+(z*(3*Nx/2+2)+x)]* RT[(YEL+3)*idx+(z*(3*Nx/2+2)+x)];
		    iu1u3 += RT[(XEL+3)*idx+(z*(3*Nx/2+2)+x)]* RT[(ZEL+3)*idx+(z*(3*Nx/2+2)+x)];
		    iu2u3 += RT[(YEL+3)*idx+(z*(3*Nx/2+2)+x)]* RT[(ZEL+3)*idx+(z*(3*Nx/2+2)+x)];
		    iu1u1 += RT[(XEL+3)*idx+(z*(3*Nx/2+2)+x)]* RT[(XEL+3)*idx+(z*(3*Nx/2+2)+x)];
		    iu2u2 += RT[(YEL+3)*idx+(z*(3*Nx/2+2)+x)]* RT[(YEL+3)*idx+(z*(3*Nx/2+2)+x)];
		    iu3u3 += RT[(ZEL+3)*idx+(z*(3*Nx/2+2)+x)]* RT[(ZEL+3)*idx+(z*(3*Nx/2+2)+x)];
		  }
	      }

	    iu1=iu1/(3*Nz/2)/(3*Nx/2);
	    iu2=iu2/(3*Nz/2)/(3*Nx/2);
	    iu3=iu3/(3*Nz/2)/(3*Nx/2);
	    iu1u2=iu1u2/(3*Nz/2)/(3*Nx/2);
	    iu1u3=iu1u3/(3*Nz/2)/(3*Nx/2);
	    iu2u3=iu2u3/(3*Nz/2)/(3*Nx/2);
	    iu1u1=iu1u1/(3*Nz/2)/(3*Nx/2);
	    iu2u2=iu2u2/(3*Nz/2)/(3*Nx/2);
	    iu3u3=iu3u3/(3*Nz/2)/(3*Nx/2);

	    us[6][y]+=iu1;
	    us[7][y]+=iu2;
	    us[8][y]+=iu3;
	    us[9][y]+=iu1u2;
	    us[10][y]+=iu1u3;
	    us[11][y]+=iu2u3;
	    us[15][y]+=iu1u1;
	    us[16][y]+=iu2u2;
	    us[17][y]+=iu3u3;
	    //fprintf(fp2, "%d, %f %f %f  %f %f  %f %f  %f  %f \n", y, iu1, iu2, iu3, iu1u2, iu1u3, iu2u3, iu1u1, iu2u2, iu3u3);

	     if(y==5)
	      {
		fprintf(afp5, " %f %f %f  %f %f  %f %f  %f  %f  %f \n", iu1, iu2, iu3, iu1u2, iu1u3, iu2u3, iu1u1, iu2u2, iu3u3, 0.5*(iu1u1+iu2u2+iu3u3));
	      }
	    if(y==40)
	      {
		fprintf(afp6, " %f %f %f  %f %f  %f %f  %f  %f  %f\n", iu1, iu2, iu3, iu1u2, iu1u3, iu2u3, iu1u1, iu2u2, iu3u3,0.5*(iu1u1+iu2u2+iu3u3) );
		}

	    au1+=u1*W[y];     au2+=u2*W[y]; 	    au3+=u3*W[y];
	    au1u2+=u1u2*W[y];  au2u3+=u2u3*W[y];  au1u3+=u1u3*W[y];
	    au1u1+=u1u1*W[y];  au2u2+=u2u2*W[y];  au3u3+=u3u3*W[y];

	    iau1+=iu1*W[y];     iau2+=iu2*W[y]; 	    iau3+=iu3*W[y];
	    iau1u2+=iu1u2*W[y];  iau2u3+=iu2u3*W[y];  iau1u3+=iu1u3*W[y];
	    iau1u1+=iu1u1*W[y];  iau2u2+=iu2u2*W[y];  iau3u3+=iu3u3*W[y];
	  }
	fprintf(afp9, " %f %f %f  %f %f  %f  %f  %f  %f \n", au1, au2, au3, au1u2, au1u3, au2u3, au1u1, au2u2, au3u3);
	fprintf(afp10, " %f %f %f  %f %f  %f %f  %f  %f \n",  iau1, iau2, iau3, iau1u2, iau1u3, iau2u3, iau1u1, iau2u2, iau3u3);
      }

 return(NO_ERR);
}

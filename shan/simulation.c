/* manufacture solution test 2 */

#include "main.h"
#include "hdf5.h"

int main(int argc, char **argv)
{
    /* External Variables.  All external variables are defined in main.h */
  extern int    qpts, dimR, dimQ, Nx, Nz, nums;
  extern double dt, re, mpg, flux,maxtstep;
    extern double *Kx, *Kz, **K2, *cfl2;
    extern double **Q, **Qp, **Qpp, **R, **Rp, **Qw, **Qpw, **Rw, **Qs, **Qps,
      **Qpps, **Rs, **Rps, *Rp0, *Rpp0;
    extern double **MZ;
    extern double ***M;
    extern double *Uadd, *Vadd, *Vpadd;
    extern double *Qy, **uu, **us;
    extern double *W;

    extern mcomplex ****U, ****C;         /* state variables */
    extern mcomplex **Fa, **Fb,	**TM, ishear;
    extern mcomplex *fa, *fb, *tm;
    extern double **MZ;
    extern double ***M;
    extern double tau, itau, atau, iatau;

    extern mcomplex ****IU, ****IC;       /* incremental state variables */
    extern mcomplex **IFa, **IFb, **ITM;
    extern mcomplex *Ifa, *Ifb, *Itm;

    extern mcomplex ****AU, ****AC;        /* adjoint variables and will use the same other variables
					      used in state equations */

    extern mcomplex ****IAU, ****IAC;      /* incremental adjoint variables */

    extern mcomplex **Uxbt, **Uzb;        /* variables used to store dux duz evaluated at y=-1 used for
					     computing boundary conditions for incremental state equations*/
    extern mcomplex **Uxb, **Uzb;         /* variables used to store dux duz evaluated at y=-1 from previous 
					     state used for boundary conditions for incremental state equations*/ 
     extern mcomplex **AUxb, **AUzb;         /* variables used to store dux duz evaluated at y=-1 from previous 
					     state used for boundary conditions for incremental state equations*/
     extern mcomplex **IUxb, **IUzb;
    extern mcomplex **IAUxb, **IAUzb;
    extern mcomplex **GUxb, **GUzb;
    extern mcomplex **GIUxb, **GIUzb;
    extern mcomplex **grad, **hess;
    extern  mcomplex **HUxb, **HUzb;
    extern  mcomplex **HAUxb, **HAUzb;
    extern fftw_complex ***CT, ***ICT;      /* variables used in fft */
    extern fftw_plan pf1, pf2;
    extern fftw_plan Ipf1, Ipf2;
    extern rfftwnd_plan pr1, pr2;

    extern mcomplex *****MC, *****MIC;     /* variables used to store state and incremental state
					      solutions between two check points. */

    extern mcomplex ****MU, ****MIU;      /* variables used to store manufacture solutions*/


    extern mcomplex ****LU, ****LIU;

    extern FILE *fp,*fp2, *fp3, *fp4, *fp5, *fp6, *fp7, *fp8, *fp9, *fp10, *fp11;
    extern FILE *afp,*afp2, *afp3, *afp4, *afp5, *afp6, *afp7, *afp8, *afp9, *afp10;

    static double e[3] = { 1./2., 2./3, 1.};
    static double h[3]={ 0., 1., 2./3.};
    static double ee[3] = { 1./3., 1./2., 1.};
    static double ss[3] = {1.,2./3., 1.};
      /* Local Variables */
    int n, z, dctr, tsteps, j, i, x, y;
    double k11=2*M_PI*1.;
    double Lx, Lz, Ny, r1,r2;
    fftw_complex *fout = NULL;
    double err,err2, norm;
   int  sizeRealTransform;
   int restart_flag;
   int count;
   int checknum, checkstep;

   maxtstep=1;

   tau=0;
   itau=0;
    if (argc != 11)
    {
	printf("Required arguments are Nx,Ny,Nz,Lx,Lz,dt,tsteps,mpg,Re.\n");
	return(EXIT_FAILURE);
    }

    Nx = atoi(argv[1]);          /* number of terms in the truncated expansion
                                    in x */
    Ny = atoi(argv[2]);          /* number of terms in the truncated expansion
                                    in y */
    Nz = atoi(argv[3]);          /* number of terms in the truncated expansion
                                    in z */
    Lx = atof(argv[4]);          /* length of the interval in x */
    Lz = atof(argv[5]);          /* length of the interval in z */
    dt = atof(argv[6]);		 /* value of time step */
    tsteps = atoi(argv[7]);      /* number of time steps to take */
    mpg = atof(argv[8]);         /* mean pressure gradient */
    re = atof(argv[9]);          /* Reynolds number */
    restart_flag=atoi(argv[10]);      /* whether to restart, 0, starting from initial conditions, 
				    nonzero, reading in stored results on hdf5 and restart from this time step */

    re = 1./re;                  /* time step routines assume I pass 1/Re */
    qpts = 3*Ny/2;               /* number of quadrature points 
                                    (see page 9 of Moser's notes) */
    dimR = Ny-2;                 /* dimR and dimQ denote the number of terms */
    dimQ = Ny-4;                 /* in the truncated expansions for the */
                                 /* functions v_hat, g_hat, U, W */
				 /* (see page 5 of Moser's notes). */
    sizeRealTransform = 3*Nx/2;  /* for the FFTs */
    // printf("Re=%f\n", re);


    //assert((Nx%4==0)&&(Nz%2==0));

    if(Nx%4!=0)
    {
	printf("Required arguments Ny/4==0\n");
	return(EXIT_FAILURE);
    }
  if(Nz%2!=0)
    {
	printf("Required arguments Nz/2==0\n");
	return(EXIT_FAILURE);
    }
  if(Ny-4<0)
    {
	printf("Required arguments Nzy>4\n");
	return(EXIT_FAILURE);
    }
    /* Create matrices using Legendre polynomials */
    if ( LegendreSetup() != NO_ERR )
    {
	return(EXIT_FAILURE);
    }

     /* Compute wave numbers */
    if ( waveNums(Nx/2, Nz, Lx, Lz) != NO_ERR )
    {
        freedMatrix(Q); freedMatrix(Qp); freedMatrix(Qpp); freedMatrix(R); 
        freedMatrix(Rp); freedMatrix(Qw); freedMatrix(Qpw); freedMatrix(Rw);
        freedMatrix(Qs); freedMatrix(Qps); freedMatrix(Qpps); freedMatrix(Rs);
        freedMatrix(Rps);freedVector(Rp0); freedVector(Vadd);freedVector(Vpadd);
	freedVector(Uadd);freedMatrix(Rpw);freedVector(Rpp0);
	return(EXIT_FAILURE);
    }

    /* get memory for 4D arrays and other matrices */
    if ( getMem() != NO_ERR )
    {
        freedMatrix(Q); freedMatrix(Qp); freedMatrix(Qpp); freedMatrix(R); 
        freedMatrix(Rp); freedMatrix(Qw); freedMatrix(Qpw); freedMatrix(Rw);
        freedMatrix(Qs); freedMatrix(Qps); freedMatrix(Qpps); freedMatrix(Rs);
        freedMatrix(Rps); freedVector(Kx); freedVector(Kz); freedMatrix(K2); 
        freedVector(Rp0); freedVector(Vadd);freedVector(Vpadd);freedMatrix(Rpw);
	freedVector(Uadd); freedVector(Qy);freedVector(Rpp0);
	return(EXIT_FAILURE);
    }

    /* Create plans for FFTs */
    pf1 = fftw_create_plan_specific(3*Nz/2, FFTW_BACKWARD, 
	FFTW_MEASURE|FFTW_IN_PLACE, CT[0][0], 3*Nx/4+1, fout, -1);
    pf2 = fftw_create_plan_specific(3*Nz/2, FFTW_FORWARD, 
	FFTW_MEASURE|FFTW_IN_PLACE, CT[0][0], 3*Nx/4+1, fout, -1);
    pr1 = rfftwnd_create_plan(1, &sizeRealTransform, FFTW_COMPLEX_TO_REAL, 
	FFTW_MEASURE|FFTW_IN_PLACE);
    pr2 = rfftwnd_create_plan(1, &sizeRealTransform, FFTW_REAL_TO_COMPLEX,
 	FFTW_MEASURE|FFTW_IN_PLACE);


    /* Create plans for FFTs */
    Ipf1 = fftw_create_plan_specific(3*Nz/2, FFTW_BACKWARD, 
	FFTW_MEASURE|FFTW_IN_PLACE, ICT[0][0], 3*Nx/4+1, fout, -1);
    Ipf2 = fftw_create_plan_specific(3*Nz/2, FFTW_FORWARD, 
	FFTW_MEASURE|FFTW_IN_PLACE, ICT[0][0], 3*Nx/4+1, fout, -1);


    /* set variables for checking CFL condition */
    if (cflVars(Lx, Lz) != 0)
    {
        printf("Error creating CF variables\n");
        fftw_destroy_plan(pf1);   fftw_destroy_plan(pf2);
        rfftwnd_destroy_plan(pr1);  rfftwnd_destroy_plan(pr2);
        freedMatrix(Q); freedMatrix(Qp); freedMatrix(Qpp); freedMatrix(R); 
        freedMatrix(Rp); freedMatrix(Qw); freedMatrix(Qpw); freedMatrix(Rw);
        freedMatrix(Qs); freedMatrix(Qps); freedMatrix(Qpps); freedMatrix(Rs);
        freedMatrix(Rps); freedVector(Kx); freedVector(Kz); freedMatrix(K2); 
        freec4Darray(U); freec4Darray(C); freec3Darray(CT); freec3Darray(ICT);
        freecMatrix(Fa); freecMatrix(Fb); freed3Darray(M); freecMatrix(TM); 
        freedMatrix(MZ); freecVector(fa); freecVector(fb); freecVector(tm);
        freedVector(cfl2);   freedVector(Rp0);freedVector(Vadd);freedVector(Vpadd);
	freedVector(Uadd);freecMatrix(ITM); freecMatrix(IFa); freecMatrix(IFb);  
         freecMatrix(Uxb); freecMatrix(Uzb);freecMatrix(Uzbt); freecMatrix(Uxbt); 
        freecVector(Ifa); freecVector(Ifb); freecVector(Itm);	freedVector(Qy);
	freec4Darray(IU); freec4Darray(IC);freec4Darray(AU); freec4Darray(AC);
	freec4Darray(IAU); freec4Darray(IAC); freec5Darray(MC); freec5Darray(MIC);
	freedMatrix(Rpw);freec4Darray(LU);freec4Darray(LIU);
	freecMatrix(AUxb);freecMatrix(AUzb); freecMatrix(IUzb);  freecMatrix(IUxb);
	freecMatrix(IAUxb);	freecMatrix(IAUzb);	freecMatrix(grad);
	freecMatrix(GUxb);freecMatrix(GUzb);freecMatrix(GIUxb);freecMatrix(GIUzb);
	freecMatrix(HUxb);freecMatrix(HUzb);	freecMatrix(HAUxb);	freecMatrix(HAUzb);
	freecMatrix(hess);freedVector(Rpp0);freedMatrix(uu);freedMatrix(us);
        return(EXIT_FAILURE);
    }

     /* initalize part */
    memset(C[0][0][0], 0, (Nz)*2*(Ny-2)*(Nx/2)*sizeof(mcomplex));
    memset(IC[0][0][0], 0, (Nz)*2*(Ny-2)*(Nx/2)*sizeof(mcomplex));
    memset(AC[0][0][0], 0, (Nz)*2*(Ny-2)*(Nx/2)*sizeof(mcomplex));
    memset(IAC[0][0][0], 0, (Nz)*2*(Ny-2)*(Nx/2)*sizeof(mcomplex));

    memset(MC[0][0][0][0], 0, (MAXNUM*3+1)*(Nz)*2*(Ny-2)*(Nx/2)*sizeof(mcomplex));
    memset(MIC[0][0][0][0], 0, (MAXNUM*3+1)*(Nz)*2*(Ny-2)*(Nx/2)*sizeof(mcomplex));

    memset(U[0][0][0], 0, (Nz)*5*qpts*(Nx/2)*sizeof(mcomplex));
    memset(AU[0][0][0], 0, (Nz)*5*qpts*(Nx/2)*sizeof(mcomplex));
    memset(IU[0][0][0], 0, (Nz)*5*qpts*(Nx/2)*sizeof(mcomplex));
    memset(IAU[0][0][0], 0, (Nz)*5*qpts*(Nx/2)*sizeof(mcomplex));
    memset(LU[0][0][0], 0, (Nz)*5*qpts*(Nx/2)*sizeof(mcomplex));
    memset(LIU[0][0][0], 0, (Nz)*5*qpts*(Nx/2)*sizeof(mcomplex));
    memset(Uxb[0], 0, Nz*(Nx/2)*sizeof(mcomplex));
    memset(Uzb[0], 0, Nz*(Nx/2)*sizeof(mcomplex));
    memset(AUxb[0], 0, Nz*(Nx/2)*sizeof(mcomplex));
    memset(AUzb[0], 0, Nz*(Nx/2)*sizeof(mcomplex));
    memset(grad[0], 0, Nz*(Nx/2)*sizeof(mcomplex));
    memset(hess[0], 0, Nz*(Nx/2)*sizeof(mcomplex));
    memset(HUxb[0], 0, Nz*(Nx/2)*sizeof(mcomplex));
    memset(HUzb[0], 0, Nz*(Nx/2)*sizeof(mcomplex));
    memset(uu[0], 0, 6*(qpts)*sizeof(double));
     memset(us[0], 0, 18*(qpts)*sizeof(double));
    count=0;
    nums=0;

    Re(C[0][ALPHA][0][0])=1.;

     for (x=0; x< Nx/8; x++)
      {
	for(z=0; z< Nz/4; z++)
	  {
	    for(y=0; y< 1; y++)
	      {
		if(x*x+z*z>0)
		  {
		      	r1 = 5*(   (double)rand() / ((double)(RAND_MAX)+(double)(1)) );
			r2=5*(   (double)rand() / ((double)(RAND_MAX)+(double)(1)) );
			//	printf("random perturbation for alpha (Kx, kz)=(%d, %d) is %f, %f\n", x, z, r1,r2); 
			Re(C[z][ALPHA][y][x])=r1;
			Im(C[z][ALPHA][y][x])=r2;
			r1 = 5*(   (double)rand() / ((double)(RAND_MAX)+(double)(1)) );
			Re(C[z][BETA][y][x])=r1;
			r2=5*(   (double)rand() / ((double)(RAND_MAX)+(double)(1)) );
			Im(C[z][BETA][y][x])=r2;
			//printf("random perturbation for beta (Kx, kz)=(%d, %d) is %f, %f\n", x, z, r1,r2); 
		  }
	      }
	  }
      }

    initAlphaBeta2();
    if (increBoundary()!= NO_ERR)
      {
	printf("increBoundary failure\n");
      }
    incre_initAlphaBeta2();
 if( restart_flag >=0)
      {
	restart2(restart_flag);
	checknum=restart_flag/MAXSTEP;
	fp=fopen("all.txt", "a+");
	fp2=fopen("home_1.txt", "a+");
	fp3=fopen("stat_uy1.txt", "a+");
	fp4=fopen("stat_uy2.txt", "a+");
	fp5=fopen("stat_iuy1.txt", "a+");
	fp6=fopen("stat_iuy2.txt", "a+");
	fp7=fopen("wall_shear.txt", "a+");
	fp8=fopen("i_wall_shear.txt", "a+");
	fp9=fopen("aver_u.txt", "a+");
	fp10=fopen("i_aver_u.txt", "a+");
	fp11=fopen("home_2.txt", "a+");
      }
    else
      {
	 write_data2(-1);
	 checknum=0;
	 fp=fopen("all.txt", "w+");
	 fp2=fopen("home_1", "w+");
	  fp3=fopen("stat_uy1.txt", "w+");
	 fp4=fopen("stat_uy2.txt", "w+");
	 fp5=fopen("stat_iuy1.txt", "w+");
	 fp6=fopen("stat_iuy2.txt", "w+");
	 fp7=fopen("wall_shear.txt", "w+");
	 fp8=fopen("i_wall_shear.txt", "w+");
	 fp9=fopen("aver_u.txt", "w+");
	 fp10=fopen("i_aver_u.txt", "w+");
	 fp11=fopen("home_2", "w+");
      }

    afp=fopen("adj_u_wall1.txt", "w+");
    afp2=fopen("adj_u_wall2.txt", "w+");
    afp3=fopen("adj_uy1.txt", "w+");
    afp4=fopen("adj_uy2.txt", "w+");
    afp5=fopen("adj_iuy1.txt", "w+");
    afp6=fopen("adj_iuy2.txt", "w+");
    afp7=fopen("adj_wall_shear.txt", "w+");
    afp8=fopen("i_adj_wall_shear.txt", "w+");
    afp9=fopen("aver_adj.txt", "w+");
    afp10=fopen("i_aver_adj.txt", "w+");

  memcpy(MC[count][0][0][0],C[0][0][0], (Nz)*2*(Ny-2)*(Nx/2)*sizeof(mcomplex));
  memcpy(MIC[count][0][0][0],IC[0][0][0], (Nz)*2*(Ny-2)*(Nx/2)*sizeof(mcomplex));

  flux=0.;

   for (i = 0; i < dimR; ++i)
    {
      for (j = 0; j < qpts; ++j)
	{
	  flux += Rw[i][j]*Re(C[0][ALPHA][i][0]); 
	}
    }
   //   printf("flux=%f\n", flux);

   /* time step for forward problem */
    for (n = restart_flag; n < tsteps; ++n)	     /* loop for each timestep */
    {
	for (dctr = 0; dctr < 3; ++dctr)     /* RK steps */
	{

	  count=(count+1)% 300;
	  // count=count+1;

	  /* copy the result to Uxbt, Uzbt. Uxb and Uzb will be used later for boundary condition
	     of current time stage */
	  memcpy(Uxbt[0], Uxb[0], (Nz)*(Nx/2)*sizeof(fftw_complex));
	  memcpy(Uzbt[0], Uzb[0], (Nz)*(Nx/2)*sizeof(fftw_complex));
	  memset(Uxb[0], 0, Nz*(Nx/2)*sizeof(fftw_complex));
	  memset(Uzb[0], 0, Nz*(Nx/2)*sizeof(fftw_complex));


	    /* do FFTs to get H_hats.  After this we have for each (Kx,y,Kz)
                  Hx_hat   -->  U[z][HXEL][y][x]
                  Hy_hat   -->  U[z][HYEL][y][x]
                  Hz_hat   -->  U[z][HZEL][y][x]
	    */
	  if (pass1(dctr, n) != NO_ERR)
	    {
		printf("Pass1 failure\n");
		n = tsteps;
		break;
	    }

	  project0(count, dctr, n, NULL);           /* (kx, kz)=(0, 0), solve for a, b */
	  project(count, n, dctr, 0, 1, NULL);      /* (kx, kz)!=(0, 0), solve for alpha, beta */
	    for (z = 1; z < Nz; ++z)
	      {
		if (z == Nz/2) 
		{
		 memset(U[z][0][0], 0, 5*qpts*(Nx/2)*sizeof(mcomplex)); 
		 continue; 
		}

		project(count, n, dctr, z, 0, NULL); 
	      }
	    memcpy(GUxb[0], Uxb[0], (Nz)*(Nx/2)*sizeof(fftw_complex));
	    memcpy(GUzb[0], Uzb[0], (Nz)*(Nx/2)*sizeof(fftw_complex));
	    /*now update the boundary condition using current time step solution of state equation*/
	    if (increBoundary()!= NO_ERR)
	      {
		printf("increBoundary failure\n");
		n = tsteps;
		break;
		}

	    increproject0(count, dctr, n, 1, NULL);
	    increproject(count, dctr, 0, 1, n, NULL);
	    for (z = 1; z < Nz; ++z)
	      {
		if (z == Nz/2) 
		{
		    // SET U[z][XEL,YEL,ZEL,DXEL,DZEL] TO ZEROS 
		  memset(IU[z][0][0], 0, 5*qpts*(Nx/2)*sizeof(mcomplex)); 
		    continue; 
		}

		increproject(count, dctr, z, 0,n, NULL); 
	      }

	}   /* end for dctr... */

	  comp_stat(n);

	/* now writing the results to HDF file at selected time steps*/
	//if (((n+1) % MAXSTEP==0) && (n+1< tsteps))           /* when we need to store the current time step results*/
	  if (((n+1) % MAXSTEP==0))
	  {
	    write_data2(n);
	    checknum=checknum+1;
	    count=0;
	    memcpy(MC[count][0][0][0],C[0][0][0], (Nz)*2*(Ny-2)*(Nx/2)*sizeof(mcomplex));
	    memcpy(MIC[count][0][0][0],IC[0][0][0], (Nz)*2*(Ny-2)*(Nx/2)*sizeof(mcomplex));
	  }
    }   /* end for n... */
    fprintf(fp, " averge from step:%d to %d\n",  restart_flag, tsteps);
    for(j=0; j< qpts; j++)
      {
	fprintf(fp, "%f %f %f %f %f %f %f %f %f %f \n",  us[0][j]/nums, us[1][j]/nums, us[2][j]/nums, us[3][j]/nums, us[4][j]/nums, us[5][j]/nums, us[12][j]/nums, us[13][j]/nums, us[14][j]/nums,us[18][j]/nums);
      }
    for(j=0; j< qpts; j++)
      {
	fprintf(fp,  " %f %f %f %f %f %f %f %f %f  \n",us[6][j]/nums, us[7][j]/nums, us[8][j]/nums, us[9][j]/nums, us[10][j]/nums, us[11][j]/nums, us[15][j]/nums, us[16][j]/nums, us[17][j]/nums);
      }

  fclose(fp3);    fclose(fp4);     fclose(fp5);     fclose(fp6);
    fclose(fp7);  fclose(fp8);  fclose(fp9);  fclose(fp10);
	/* start restart and solving for adjoint system */
   memset(us[0], 0, 18*(qpts)*sizeof(double));
  for(checkstep=checknum-1; checkstep>=0; checkstep--)
    {

      restart_flag=checkstep*MAXSTEP;
      printf("\n");
      printf(" restarting from time step: %d\n", restart_flag);
      restart2(restart_flag);
      count=0;

      memcpy(MC[count][0][0][0],C[0][0][0], (Nz)*2*(Ny-2)*(Nx/2)*sizeof(mcomplex));
      memcpy(MIC[count][0][0][0],IC[0][0][0], (Nz)*2*(Ny-2)*(Nx/2)*sizeof(mcomplex));
      for (n = restart_flag; n <restart_flag+MAXSTEP; ++n)	     /* loop for each timestep */
	{
	  for (dctr = 0; dctr < 3; ++dctr)     /* RK steps */
	    {

	      count=count+1;

	      /* copy the result to Uxbt, Uzbt. Uxb and Uzb will be used later for boundary condition
		 of current time stage */
	      memcpy(Uxbt[0], Uxb[0], (Nz)*(Nx/2)*sizeof(fftw_complex));
	      memcpy(Uzbt[0], Uzb[0], (Nz)*(Nx/2)*sizeof(fftw_complex));
	      memset(Uxb[0], 0, Nz*(Nx/2)*sizeof(fftw_complex));
	      memset(Uzb[0], 0, Nz*(Nx/2)*sizeof(fftw_complex));

	      if (pass1(dctr, n) != NO_ERR)
		{
		  printf("Pass1 failure\n");
		  checkstep=-1;
		return(EXIT_FAILURE);
		}

	      project0(count, dctr, n, NULL);           /* (kx, kz)=(0, 0), solve for a, b */
	      project(count, n, dctr, 0, 1, NULL);      /* (kx, kz)!=(0, 0), solve for alpha, beta */
	      for (z = 1; z < Nz; ++z)
		{
		  if (z == Nz/2) 
		    {
		      memset(U[z][0][0], 0, 5*qpts*(Nx/2)*sizeof(mcomplex)); 
		      continue; 
		    }

		  project(count, n, dctr, z, 0, NULL); 
		}

	      /*now update the boundary condition using current time step solution of state equation*/
	      if (increBoundary()!= NO_ERR)
		{
		  printf("increBoundary failure\n");
		  checkstep=-1;
		  return(EXIT_FAILURE);
		}

	      increproject0(count, dctr, n, 1, NULL);
	      increproject(count, dctr, 0, 1, n, NULL);
	      for (z = 1; z < Nz; ++z)
		{
		  if (z == Nz/2) 
		    {
		    // SET U[z][XEL,YEL,ZEL,DXEL,DZEL] TO ZEROS 
		      memset(IU[z][0][0], 0, 5*qpts*(Nx/2)*sizeof(mcomplex)); 
		      continue; 
		    }

		  increproject(count, dctr, z, 0,n, NULL); 
		}

	    }   /* end for dctr... */
	}

	  memcpy(Uxb[0], AUxb[0], (Nz)*(Nx/2)*sizeof(fftw_complex));
	  memcpy(Uzb[0], AUzb[0], (Nz)*(Nx/2)*sizeof(fftw_complex));
	  printf(" start adjoint solve at time step %d\n", restart_flag+MAXSTEP);

    /* time step for backward equations */
	  for (n =restart_flag+MAXSTEP; n >restart_flag; --n)	     /* loop for each timestep */
	    {
	      for (dctr = 0; dctr < 3; ++dctr)     /* RK steps */
		{

	  /* copy the result to Uxbt, Uzbt. Uxb and Uzb will be used later for boundary condition
	     of current time stage */
		  memcpy(Uxbt[0], Uxb[0], (Nz)*(Nx/2)*sizeof(fftw_complex));
		  memcpy(Uzbt[0], Uzb[0], (Nz)*(Nx/2)*sizeof(fftw_complex));
		  memset(Uxb[0], 0, Nz*(Nx/2)*sizeof(fftw_complex));
		  memset(Uzb[0], 0, Nz*(Nx/2)*sizeof(fftw_complex));
		  memset(AUxb[0], 0, Nz*(Nx/2)*sizeof(fftw_complex));
		  memset(AUzb[0], 0, Nz*(Nx/2)*sizeof(fftw_complex));

		  count=count-1;

		  /*read data from memery */
		  if(dctr<2)
		    {
		      memcpy(C[0][0][0],MC[count-1][0][0][0], (Nz)*2*(Ny-2)*(Nx/2)*sizeof(mcomplex));
		      memcpy(IC[0][0][0],MIC[count-1][0][0][0], (Nz)*2*(Ny-2)*(Nx/2)*sizeof(mcomplex));


		      /*reconstruct the state and incremental state solution u, iu from alpha and beta */
		      initAlphaBeta2();
		      if (increBoundary()!= NO_ERR)
			{
			  printf("increBoundary failure\n");
			}
		      incre_initAlphaBeta2();

		      if (pass2(dctr, n) != NO_ERR)
			{
			  printf("Pass2 failure\n");
			  n =restart_flag-1;
			  return(EXIT_FAILURE);
			}
		      memcpy(LU[0][0][0], U[0][0][0], (Nz)*5*qpts*(Nx/2)*sizeof(mcomplex));
		      memcpy(LIU[0][0][0], IU[0][0][0], (Nz)*5*qpts*(Nx/2)*sizeof(mcomplex));

		    }
		  /*read data from memery */
		  memcpy(C[0][0][0],MC[count][0][0][0], (Nz)*2*(Ny-2)*(Nx/2)*sizeof(mcomplex));
		  memcpy(IC[0][0][0],MIC[count][0][0][0], (Nz)*2*(Ny-2)*(Nx/2)*sizeof(mcomplex));


		  /*reconstruct the state and incremental state solution u, iu from alpha and beta */
		  initAlphaBeta2();
		  memcpy(GUxb[0], Uxb[0], (Nz)*(Nx/2)*sizeof(fftw_complex));
		  memcpy(GUzb[0], Uzb[0], (Nz)*(Nx/2)*sizeof(fftw_complex));

		  if (increBoundary()!= NO_ERR)
		    {
		      printf("increBoundary failure\n");
		    }
		  incre_initAlphaBeta2();

		  memcpy(AUxb[0], Uxb[0], (Nz)*(Nx/2)*sizeof(fftw_complex));
		  memcpy(AUzb[0], Uzb[0], (Nz)*(Nx/2)*sizeof(fftw_complex));

		  //  memset(IU[0][0][0], 0, (Nz)*5*qpts*(Nx/2)*sizeof(mcomplex));
		  if (pass2(dctr, n) != NO_ERR)
		    {
		      printf("Pass2 failure\n");
		      n=restart_flag-1;
		      return(EXIT_FAILURE);
		    }

		  memset(Uxb[0],  0, Nz*(Nx/2)*sizeof(mcomplex));
		  memset(Uzb[0], 0, Nz*(Nx/2)*sizeof(mcomplex));

		  adjproject0(dctr, n, count, NULL);
		  adjproject(n, dctr, 0, 1, count, NULL);
		  for (z = 1; z < Nz; ++z)
		    {
		      if (z == Nz/2) 
			{
			  memset(AU[z][0][0], 0, 5*qpts*(Nx/2)*sizeof(mcomplex)); 
			  continue; 
			}

		      adjproject(n, dctr, z, 0, count, NULL); 
		    }

		  memcpy(GIUxb[0], Uxb[0], (Nz)*(Nx/2)*sizeof(fftw_complex));
		  memcpy(GIUzb[0], Uzb[0], (Nz)*(Nx/2)*sizeof(fftw_complex));
		  /*now update the boundary condition using current time step solution of state equation*/
		  if (increBoundary()!= NO_ERR)
		    {
		      printf("increBoundary failure\n");
		      n=restart_flag-1;
		      return(EXIT_FAILURE);
		    }

		  increadjproject0(dctr, n, count,NULL);
		  increadjproject(n, dctr, 0, 1, count, NULL);
		  for (z = 1; z < Nz; ++z)
		    {
		      if (z == Nz/2) 
			{
			  // SET U[z][XEL,YEL,ZEL,DXEL,DZEL] TO ZEROS 
			  memset(IAU[z][0][0], 0, 5*qpts*(Nx/2)*sizeof(mcomplex)); 
			  continue; 
			}

		      increadjproject(n, dctr, z, 0, count, NULL); 
		    }


	      if(dctr==2)
	      {
		comp_gradient(n, dctr);
		comp_hess(n, dctr);
		}

		  if( n==restart_flag+1 && dctr==2)
		    {
		      memcpy(AUxb[0], Uxb[0], (Nz)*(Nx/2)*sizeof(fftw_complex));
		      memcpy(AUzb[0], Uzb[0], (Nz)*(Nx/2)*sizeof(fftw_complex));
		    }
		}
	      comp_stat2(n);

	    }
    }

  //    printf("for the adjoint solver \n\n\n");

      fprintf(fp, " averge from step:%d to %d\n",  restart_flag, tsteps);
    	for(j=0; j< qpts; j++)
	  {
	    fprintf(fp, "%f %f %f %f %f %f %f %f %f \n",  us[0][j]/nums, us[1][j]/nums, us[2][j]/nums, us[3][j]/nums, us[4][j]/nums, us[5][j]/nums, us[12][j]/nums, us[13][j]/nums, us[14][j]/nums);
	  }
	for(j=0; j< qpts; j++)
	  {
	    fprintf(fp,  " %f %f %f %f %f %f %f %f %f  \n",us[6][j]/nums, us[7][j]/nums, us[8][j]/nums, us[9][j]/nums, us[10][j]/nums, us[11][j]/nums, us[15][j]/nums, us[16][j]/nums, us[17][j]/nums);
	  }


	for (z=0; z< Nz; z++)
	  {
	    for(x=0; x< Nx; x++)
	      {
		fprintf(fp, " gradient[%d][%d]=%f+i%f\n", z, x,Re(grad[z][x]), Im(grad[z][x]));
	      }
	  }

	fclose(fp);
	  //fclose(fp);     fclose(fp2);
	// fclose(fp3);    fclose(fp4);     fclose(fp5);     fclose(fp6);
	//fclose(fp7);  fclose(fp8);  fclose(fp9);  fclose(fp10);
	fclose(afp);    fclose(afp2);
    fclose(afp3);    fclose(afp4);     fclose(afp5);    fclose(afp6);
     fclose(afp7);  fclose(afp8);  fclose(afp9);  fclose(afp10);
    /* clean up... */
    fftw_destroy_plan(pf1);   fftw_destroy_plan(pf2);
    rfftwnd_destroy_plan(pr1);  rfftwnd_destroy_plan(pr2);
    fftw_destroy_plan(Ipf1);   fftw_destroy_plan(Ipf2);

    freedMatrix(Q); freedMatrix(Qp); freedMatrix(Qpp); freedMatrix(R); 
    freedMatrix(Rp); freedMatrix(Qw); freedMatrix(Qpw); freedMatrix(Rw);
    freedMatrix(Qs); freedMatrix(Qps); freedMatrix(Qpps); freedMatrix(Rs);
    freedMatrix(Rps); freedVector(Kx); freedVector(Kz); freedMatrix(K2); 
    freec4Darray(U); freec4Darray(C); freec3Darray(CT); 
    freecMatrix(Fa); freecMatrix(Fb); freed3Darray(M); freecMatrix(TM); 
    freedMatrix(MZ); freecVector(fa); freecVector(fb); freecVector(tm);
    freedVector(cfl2);  freedVector(Rp0);  freec3Darray(ICT);
    freec4Darray(IU); freec4Darray(IC); freecMatrix(IFa); freecMatrix(IFb);
    freecMatrix(ITM); freecVector(Ifa); freecVector(Ifb); freecVector(Itm);
    freecMatrix(Uxbt); freecMatrix(Uzbt); freecMatrix(Uxb); freecMatrix(Uzb);
    freedVector(Vadd);freedVector(Vpadd); freedVector(Uadd);
    freec4Darray(AU); freec4Darray(AC);    freec4Darray(IAU); freec4Darray(IAC);
    freec5Darray(MC); freec5Darray(MIC);  freec4Darray(MU); freec4Darray(MIU);
    freec4Darray(LU);freec4Darray(LIU);
    freecMatrix(AUxb);freecMatrix(AUzb);freecMatrix(IUzb);  freecMatrix(IUxb);
    freecMatrix(IAUxb);	freecMatrix(IAUzb);	freecMatrix(grad);
    freecMatrix(GUxb);freecMatrix(GUzb);freecMatrix(GIUxb);freecMatrix(GIUzb);
    freecMatrix(HUxb);freecMatrix(HUzb);	freecMatrix(HAUxb);	freecMatrix(HAUzb);
    freecMatrix(hess);freedVector(Rpp0);freedMatrix(uu);freedMatrix(us);
    return(EXIT_SUCCESS);

}	/* end main */


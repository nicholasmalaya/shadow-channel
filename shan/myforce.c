#include <math.h>
#include <stdio.h>
#include <string.h>
#include "minChnl.h"

/**************forces used for manufacture solution 1 ***************************/
void  force0(int n, int k, int flag, double *f, double *g)
{
  extern double *Qy, **Rw;
  extern double dt,mpg, re;
  extern int dimR, qpts;
  static double e[3] = { 0., 1./2., 2./3.};

  int i, j;
  for (i = 0; i < dimR; ++i)
    {
      f[i]=0;
      g[i]=0;
      if(flag==0)
	{
        for (j = 0; j < qpts; ++j)
	{
	  //f[i] += Rw[i][j]*((1-Qy[j]*Qy[j])*cos(n*dt+e[k]*dt)
	  //  +2.*re*sin(n*dt+dt*e[k])+ mpg);
	  f[i] += Rw[i][j]*((1-Qy[j]*Qy[j])*cos(2*M_PI*(n*dt+e[k]*dt))*2*M_PI
			    +2.*re*sin(2*M_PI*(n*dt+dt*e[k]))+ mpg);
	}
	}
    }

}


void  increforce0(int n, int k, int flag, double *f, double *g)
{
  extern double *Qy, **Rw;
  extern double dt,mpg, re;
  extern int dimR, qpts;
  static  double e[3] = { 0., 1./2., 2./3.};

  int i, j;
  for (i = 0; i < dimR; ++i)
    {
      f[i]=0;
      g[i]=0;
        for (j = 0; j < qpts; ++j)
	{
	 f[i] +=  Rw[i][j]*((2-Qy[j]-Qy[j]*Qy[j])*cos(n*dt+e[k]*dt)
			    +2*re*sin(n*dt+e[k]*dt));
	}
    }

}

void increforce(int n, int k, int z, double *f, double *g)
{
  extern double *Qy, **Rw, *Kz;
  extern double dt,mpg, re;
  extern int dimR, qpts, dimQ;
  static  double e[3] = { 0., 1./2., 2./3.};
  extern int Nz;
  int i, j;

  if(z==1)
    {
      for (i = 0; i < dimR; ++i)
	{
	  f[i]=0.;
	  for (j = 0; j < qpts; ++j)
	    {
	      f[i] += Rw[i][j] * Kz[z]*(1-Qy[j])*
		( cos(n*dt+e[k]*dt)+re*Kz[z]*Kz[z]*sin(n*dt+e[k]*dt));
	    }
	}
    }
  else if(z==Nz-1)
      {
	for (i = 0; i < dimR; ++i)
	  {
	    f[i]=0.;
	    for (j = 0; j < qpts; ++j)
	      {
	        f[i] += Rw[i][j] * Kz[z]*(1-Qy[j])*
		  ( cos(n*dt+e[k]*dt)+re*Kz[z]*Kz[z]*sin(n*dt+e[k]*dt));
	      }
	  }
      }
  else
    {
      for (i = 0; i < dimR; ++i)
	{
	  f[i]=0;
	}
    }

  for (i = 0; i < dimQ; ++i)
    {
      g[i]=0.;
    }
}
/******************************end of forces for manufacture solution 1 *************/

/*******************************force used for manufacture solution 2****************/
void  force0_2(int n, int k, int flag, double *f, double *g)
{
  extern double *Qy, **Rw;
  extern double dt,mpg, re;
  extern int dimR, qpts;
    static double e[3] = { 0., 1./2, 2./3};
  //static  double e[3] = { 0., 8./15., 2./3.};
  int i, j;
  
  for (i = 0; i < dimR; ++i)
    {
      f[i]=0;
      g[i]=0;
      if(flag==0)
	{
	  for (j = 0; j < qpts; ++j)
	    {
	      f[i] +=  Rw[i][j]*sin(2.*M_PI*Qy[j])*(cos(n*dt+e[k]*dt)
						    +re*4.*M_PI*M_PI*sin(n*dt+dt*e[k]))+Rw[i][j]*mpg;
	    }
	}
    }

}

void  increforce0_2(int n, int k, int flag, double *f, double *g)
{
  extern double *Qy, **Rw;
  extern double dt,mpg, re;
  extern int dimR, qpts;
   static  double e[3] = { 0., 1./2., 2./3.};
   // static  double e[3] = { 0., 8./15., 2./3.};
  int i, j;
  for (i = 0; i < dimR; ++i)
    {
      f[i]=0;
      g[i]=0;
      if(flag==1)
	{
	  for (j = 0; j < qpts; ++j)
	    { 
	      f[i]+= Rw[i][j]*2*M_PI*sin(M_PI*Qy[j]/4+3*M_PI/4)*
		(cos(n*dt+e[k]*dt)+re*M_PI*M_PI/16*sin(n*dt+e[k]*dt));
	    }
	}
    }
  
}

void increforce_2(int n, int k, int z, double *tmp, double *tmp2)
{
  extern double *Qy, **Rw, *Kz;
  extern double dt,mpg, re;
  extern int dimR, qpts, dimQ;
   static  double h[3] = { 0., 1./2., 2./3.};
   // static  double h[3] = { 0., 8./15., 2./3.};
  extern int Nz;

  int i, j;

  for (i = 0; i < dimR; ++i)
    {
      tmp[i]=0;
    }
  for (i = 0; i < dimQ; ++i)
    {
      tmp2[i]=0;
    }
  if(z==1 ||z==Nz-1)
      {
	for (i = 0; i < dimR; ++i)
	  {
	    tmp[i]=0;
	    for (j = 0; j < qpts; ++j)
	      {
		tmp[i]+=Rw[i][j] * Kz[z]*2.*M_PI*sin(M_PI*Qy[j]/4+3.*M_PI/4)*
		  ( cos(n*dt+h[k]*dt)+re*( M_PI*M_PI/16.+Kz[z]*Kz[z])*sin(n*dt+h[k]*dt));
	      }
	  }
      }
}

void  force0_3(int n, int k, int z, double *f, double *g)
{
  extern double *Qy, **Rw;
  extern double dt,mpg, re;
  extern int dimR, qpts;
  static double e[3] = { 0., 1./2., 2./3};

  int i, j;
  for (i = 0; i < dimR; ++i)
    {
      f[i]=0;
      g[i]=0.;
        for (j = 0; j < qpts; ++j)
	{
	  f[i] += Rw[i][j]*mpg;
	}
    }

}
void  force_3(int n, int k, int z, double *f, double *g)
{
  extern double *Qy, **Rw;
  extern double dt,mpg, re;
  extern int dimR, qpts, dimQ;
  // static double h[3] = { 0., 8./15, 2./3};
   static double h[3] = { 0., 1./2., 2./3};
  double k11=2.*M_PI;
  double k22=2.*M_PI;

  int i, j;
   for (i = 0; i < dimQ; ++i)
    {
      g[i]=0;
    }
  for (i = 0; i < dimR; ++i)
    {
      f[i]=0;
      if(z==1)
	{
	  for (j = 0; j < qpts; ++j)
	      {
		f[i]+= Rw[i][j]*(k11*k11+k22*k22)*sin(M_PI*Qy[j]*2)*
		  ( cos(n*dt+h[k]*dt)+re*(4*M_PI*M_PI+k11*k11+k22*k22)*sin(n*dt+h[k]*dt));
	      }
	}
    }
}
void  increforce0_3(int n, int k, int flag, double *f, double *g)
{
  extern double *Qy, **Rw;
  extern double dt,mpg, re;
  extern int dimR, qpts,dimQ;
  //  static  double e[3] = { 0., 8./15, 2./3};
  static double e[3] = { 0., 1./2., 2./3}; 
  double k11=2.*M_PI;
  double k22=2.*M_PI;
  int i, j;

  for (i = 0; i < dimR; ++i)
    {
      f[i]=0;
      g[i]=0;
        for (j = 0; j < qpts; ++j)
	{
	  f[i]+=k22* Rw[i][j]*4*M_PI*sin(M_PI*Qy[j]/4+3*M_PI/4)*
	       (cos(n*dt+e[k]*dt)+re*M_PI*M_PI/16*sin(n*dt+e[k]*dt));
	  g[i]+=-k11* Rw[i][j]*4*M_PI*sin(M_PI*Qy[j]/4+3*M_PI/4)*
	       (cos(n*dt+e[k]*dt)+re*M_PI*M_PI/16*sin(n*dt+e[k]*dt));
	}
    }

}

void increforce_3(int n, int k, int z, double *tmp, double *tmp2)
{
  extern double *Qy, **Rw, *Kz;
  extern double dt,mpg, re;
  extern int dimR, qpts, dimQ;
  // static  double h[3] = { 0., 8./15, 2./3};
  static double h[3] = { 0., 1./2., 2./3}; 
  extern int Nz;
  double k11=2.*M_PI;
  double k22=2.*M_PI;

  int i, j;

  for (i = 0; i < dimQ; ++i)
    {
      tmp2[i]=0;
    }
  for (i = 0; i < dimR; ++i)
    {
      tmp[i]=0;
    }

  if(z==2)
      {
	for (i = 0; i < dimR; ++i)
	  {
	    tmp[i]=0;
	    for (j = 0; j < qpts; ++j)
	      {
		tmp[i]+=Rw[i][j] *4.*M_PI*sin(M_PI*Qy[j]/4+3.*M_PI/4)*(k11*k11+k22*k22)*
		  ( cos(n*dt+h[k]*dt)+re*( M_PI*M_PI/16.+4*(k11*k11+k22*k22))*sin(n*dt+h[k]*dt));
	      }
	  }
      }
}


void  force0_4(int n, int k, int z, double *f, double *g)
{
  extern double *Qy, **Rw;
  extern double dt,mpg, re;
  extern int dimR, qpts;
  //  static double e[3] = { 0., 8./15, 2./3};
   static double e[3] = { 0., 1./2., 2./3};
  double k11=2.*M_PI;
   double k22=2.*M_PI;
   int i, j;

   for (i = 0; i < dimR; ++i)
     {
       f[i]=0;
       for (j = 0; j < qpts; ++j)
	 {
	   f[i] += Rw[i][j]*mpg+ Rw[i][j]* (k22*(1-Qy[j]*Qy[j])*cos(n*dt+e[k]*dt)
					    +re*2.*k22*sin(n*dt+dt*e[k]));
	 }
     }

   for (i = 0; i < dimR; ++i)
     {
       g[i]=0;
       for (j = 0; j < qpts; ++j)
	 {
	   g[i] -=Rw[i][j]* (k11*(1-Qy[j]*Qy[j])*cos(n*dt+e[k]*dt)
					       +re*2.*k11*sin(n*dt+dt*e[k]));
	 }
     }
}

void increforce_4(int n, int k, int z, double *tmp, double *tmp2)
{
  extern double *Qy, **Rw,**Qw, *Kz;
  extern double dt,mpg, re;
  extern int dimR, qpts, dimQ;
  // static  double h[3] = { 0., 8./15, 2./3};
   static  double h[3] = { 0., 1./2., 2./3};
  extern int Nz;
  double k11=2.*M_PI;
  double k22=2.*M_PI;

  int i, j;

  for (i = 0; i < dimR; ++i)
    {
      tmp[i]=0;
    }
 for (i = 0; i < dimQ; ++i)
    {
      tmp2[i]=0;
    }
  if(z==1)
      {
	for (i = 0; i < dimQ; ++i)
	  {
	    tmp2[i]=0;
	    for (j = 0; j < qpts; ++j)
	      {
		tmp2[i]+=Qw[i][j] *k11* ((12*Qy[j]*Qy[j]-4-(k11*k11+k22*k22)*(1-Qy[j]*Qy[j])*(1-Qy[j]*Qy[j]))*cos(n*dt+h[k]*dt)
					-re*( 24-2*(k11*k11+k22*k22)*(12*Qy[j]*Qy[j]-4)+
					      (k11*k11+k22*k22)*(k11*k11+k22*k22)*(1-Qy[j]*Qy[j])*(1-Qy[j]*Qy[j]))*sin(n*dt+h[k]*dt));
	      }
	  }
 	for (i = 0; i < dimR; ++i)
	  {
	    tmp[i]=0;
	    for (j = 0; j < qpts; ++j)
	      {
		tmp[i]+=Rw[i][j] *(((k11*k11+k22*k22)* (1-Qy[j])+k22*(1-Qy[j]*Qy[j])*4*Qy[j])*cos(n*dt+h[k]*dt)
				   -re*(-k22*24*Qy[j]-(k11*k11+k22*k22)*((k11*k11+k22*k22)*(1-Qy[j])+k22*(1-Qy[j]*Qy[j])*4*Qy[j]))*sin(n*dt+h[k]*dt));
	      }
	  }
      }

   if(z==-1)
      {
	for (i = 0; i < dimR; ++i)
	  {
	    tmp[i]=0;
	    for (j = 0; j < qpts; ++j)
	      {
		tmp[i]+=Rw[i][j] *2*Qy[j]*k11*(k11*k11+k22*k22)*(1-Qy[j]*Qy[j])*(1-Qy[j]*Qy[j])*sin(n*dt+h[k]*dt)*sin(n*dt+h[k]*dt);
	      }
	  }
      }


}
void  force0_5(int n, int k, int z, double *f)
{
  extern double *Qy, **Rw;
  extern double dt,mpg, re;
  extern int dimR, qpts;
  static double e[3] = { 0., 8./15, 2./3};
  double k11=2.*M_PI;
   double k22=2.*M_PI;
  int i, j;
  if(z==0)
    {
      for (i = 0; i < dimR; ++i)
	{
	  f[i]=0;
	  for (j = 0; j < qpts; ++j)
	    {
	      f[i] += Rw[i][j]*mpg;
	    }
	}
    }
  if(z==1)
    {
      for (i = 0; i < dimR; ++i)
	{
	  f[i]=0;
	  for (j = 0; j < qpts; ++j)
	    {
	      f[i] =0;
	    }
	}
    }
}

void increforce_5(int n, int k, int z, double *tmp)
{
  extern double *Qy, **Rw,**Qw, *Kz;
  extern double dt,mpg, re;
  extern int dimR, qpts, dimQ;
  static  double h[3] = { 0., 8./15, 2./3};
  extern int Nz;
  double k11=2.*M_PI;
  double k22=2.*M_PI;

  int i, j;

  for (i = 0; i < dimR; ++i)
    {
      tmp[i]=0;
    }

  if(z==0)
      {
	for (i = 0; i < dimQ; ++i)
	  {
	    tmp[i]=0;
	    for (j = 0; j < qpts; ++j)
	      {
		tmp[i]+=Qw[i][j] *k11* ((12*Qy[j]*Qy[j]-4-(k11*k11+k22*k22)*(1-Qy[j]*Qy[j])*(1-Qy[j]*Qy[j]))*cos(n*dt+h[k]*dt)
					-re*( 24-2*(k11*k11+k22*k22)*(12*Qy[j]*Qy[j]-4)+
					      (k11*k11+k22*k22)*(k11*k11+k22*k22)*(1-Qy[j]*Qy[j])*(1-Qy[j]*Qy[j]))*sin(n*dt+h[k]*dt));
	      }
	  }
      }

  
  if(z==1)
      {
	for (i = 0; i < dimR; ++i)
	  {
	    tmp[i]=0;
	    for (j = 0; j < qpts; ++j)
	      {
		tmp[i]+=Rw[i][j] *k22 *( (1-Qy[j]*Qy[j])*4*Qy[j]*cos(n*dt+h[k]*dt)
		  +re*(24*Qy[j]+(k11*k11+k22*k22)*(1-Qy[j]*Qy[j])*4*Qy[j])*sin(n*dt+h[k]*dt));
	      }
	  }
      }

}

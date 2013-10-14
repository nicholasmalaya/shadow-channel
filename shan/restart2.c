#include "main.h"
#include "hdf5.h"
#include <stdio.h>
#include <string.h>
#include "minChnl.h"
#include "mvOps.h"
#include "arrays.h"
void restart2(int restart_flag)
{

  extern int    qpts, dimR, dimQ, Nx, Nz;
  extern mcomplex ****U, ****C;
  extern mcomplex ****IU, ****IC;
  extern double **Rw, **Rs, **MZ, *Rp0;
  extern mcomplex **Uxb, **Uzb;
  int x, y, z, i, j;
  hid_t  file_id1, dataset_a, dataset_b;   /* file identifier */
  hid_t fid1,dataset_ia, dataset_ib ;
  hid_t complex_id;
  hsize_t fdim[]={ dimR, Nz, Nx/2};
  herr_t  ret;
  char filename [50];

    /* define compound datatype for the complex number */
    typedef struct
    {
      double re;   /*real part*/
      double im;   /*imaginary part*/
    } complex_t;

    complex_id = H5Tcreate (H5T_COMPOUND, sizeof (complex_t));
    H5Tinsert (complex_id, "real", HOFFSET(complex_t,re), H5T_NATIVE_DOUBLE);
    H5Tinsert (complex_id, "imaginary", HOFFSET(complex_t,im), H5T_NATIVE_DOUBLE);

    /* define some temporal matrix to store the data to the hdf file */
    complex_t Matrix1[dimR][Nz][Nx/2];
    complex_t Matrix2[dimR][Nz][Nx/2];

    complex_t IMatrix1[dimR][Nz][Nx/2];
    complex_t IMatrix2[dimR][Nz][Nx/2];


      sprintf(filename, "data_t=%d", restart_flag);

      // open the file and dataset 
      file_id1 = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
      dataset_a = H5Dopen(file_id1,"/data_alpha");
      dataset_b = H5Dopen(file_id1,"/data_beta");

      dataset_ia = H5Dopen(file_id1,"/data_ialpha");
      dataset_ib = H5Dopen(file_id1,"/data_ibeta");

      ret = H5Dread(dataset_a, complex_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, Matrix1);
      ret = H5Dread(dataset_b, complex_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, Matrix2);

      ret = H5Dread(dataset_ia, complex_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, IMatrix1);
      ret = H5Dread(dataset_ib, complex_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, IMatrix2);


      ret=H5Dclose(dataset_a);
      ret=H5Dclose(dataset_b);

      ret=H5Dclose(dataset_ia);
      ret=H5Dclose(dataset_ib);

      ret=H5Fclose(file_id1);

      for (y=0; y<dimR; y++)
	{
	  for (z=0; z<Nz; ++z)
	    {
	      for(x=0; x<Nx/2; ++x)
		{
		  Re(C[z][ALPHA][y][x])=Matrix1[y][z][x].re;  
		  Im(C[z][ALPHA][y][x])=Matrix1[y][z][x].im;
		  Re(C[z][BETA][y][x])=Matrix2[y][z][x].re;
		  Im(C[z][BETA][y][x])=Matrix2[y][z][x].im;

		  Re(IC[z][ALPHA][y][x])=IMatrix1[y][z][x].re; 
		  Im(IC[z][ALPHA][y][x])=IMatrix1[y][z][x].im;
		  Re(IC[z][BETA][y][x])=IMatrix2[y][z][x].re;
		  Im(IC[z][BETA][y][x])=IMatrix2[y][z][x].im;
		}
	    }
	} 

      /* compute  ux hat, uy hat, uz hat given alpha, beta,. */
           initAlphaBeta2();

      /* compute the boundary condition for previous stage or time step, result is stored in Uxb and Uzb*/
	 if (increBoundary()!= NO_ERR)
	{
	  printf("increBoundary failure\n");
	}

      /* compute iux, iuy, iuz given i-alpha, i-beta */ 

       incre_initAlphaBeta2();

       memset(U[Nz/2][0][0], 0, 5*qpts*(Nx/2)*sizeof(mcomplex)); 
       memset(IU[Nz/2][0][0], 0, 5*qpts*(Nx/2)*sizeof(mcomplex)); 

       /*
	    for (i = 0; i < dimQ; ++i)
		{
		  for (z= 0; z < Nz; ++z)
		    {
		      for (x = 0; x < Nx/2; ++x)
			{
			  printf(" C[%d][ALPHA][%d][%d]=%25.16e+i%25.16e\n", z, i, x,Re(C[z][ALPHA][i][x]), Im(C[z][ALPHA][i][x]));
			  printf(" C[%d][BETA][%d][%d]=%25.16e+%25.16ei\n", z, i, x,Re(C[z][BETA][i][x]), Im(C[z][BETA][i][x]));
			}
		    }
		}
	      for (i = 0; i < dimR; ++i)
		{
		  for (z= 0; z < Nz; ++z)
		    {
		      for (x = 0; x < Nx/2; ++x)
			{
			  printf(" IC[%d][ALPHA][%d][%d]=%25.16e+%25.16ei\n", z, i, x,Re(IC[z][ALPHA][i][x]), Im(IC[z][ALPHA][i][x]));
			  printf(" IC[%d][BETA][%d][%d]=%25.16e+%25.16ei\n", z, i, x,Re(IC[z][BETA][i][x]), Im(IC[z][BETA][i][x]));
			}
		    }
		}
	      for (i = 0; i < 5; ++i)
		{
		  for (z= 0; z < Nz; ++z)
		    {
		      for (x = 0; x < Nx/2; ++x)
			{
			  for (y=0; y< qpts; y++)
			    {
			      printf(" U[%d][%i][%d][%d]=%25.16e+%25.16ei\n", z, i, y, x,Re(U[z][i][y][x]), Im(U[z][i][y][x]));
			      printf(" IU[%d][%i][%d][%d]=%25.16e+%25.16ei\n", z, i, y, x,Re(IU[z][i][y][x]), Im(IU[z][i][y][x]));
			    }
			}
		    }
		}

	      for (z=0; z< Nz; z++)
		{
		  for(x=0; x< Nx/2; x++)
		    {
		      printf("Uxb[%d][%d]=%25.16e+i%25.16e, Uzb=%25.16e+i%25.16e\n", z, x, Re(Uxb[z][x]), Im(Uxb[z][x]), Re(Uzb[z][x]), Im(Uzb[z][x]));
		    }
		}
       */
}

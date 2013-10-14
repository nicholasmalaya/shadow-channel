#include "main.h"
#include "hdf5.h"
#include <stdio.h>
#include <string.h>
#include "minChnl.h"
#include "mvOps.h"
#include "arrays.h"
int  write_data(int n)
{
    hid_t  file_id1, dataset_u, dataset_v,dataset_w;   /* file identifier */
    hid_t fid1,dataset_iu, dataset_iv,dataset_iw ;
    hid_t complex_id;
    hsize_t fdim[]={ qpts, Nz, Nx/2};
    herr_t  ret;
    char filename [50];

    extern mcomplex ****U, ****IU;
    extern int    qpts, Nx, Nz;
    int x, y,z, i;

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
    complex_t Matrix1[qpts][Nz][Nx/2];
    complex_t Matrix2[qpts][Nz][Nx/2];
    complex_t Matrix3[qpts][Nz][Nx/2];

    complex_t IMatrix1[qpts][Nz][Nx/2];
    complex_t IMatrix2[qpts][Nz][Nx/2];
    complex_t IMatrix3[qpts][Nz][Nx/2]; 
    for (y=0; y<qpts; y++)
      {
	for (z=0; z<Nz; ++z)
	  {
	    for(x=0; x<Nx/2; ++x)
	      {
		Matrix1[y][z][x].re=Re(U[z][XEL][y][x]);  /* storing solved solutions to Matrix */
		Matrix1[y][z][x].im=Im(U[z][XEL][y][x]);
		Matrix2[y][z][x].re=Re(U[z][YEL][y][x]);
		Matrix2[y][z][x].im=Im(U[z][YEL][y][x]);
		Matrix3[y][z][x].re=Re(U[z][ZEL][y][x]);
		Matrix3[y][z][x].im=Im(U[z][ZEL][y][x]);

		IMatrix1[y][z][x].re=Re(IU[z][XEL][y][x]);  /* storing solved solutions to Matrix */
		IMatrix1[y][z][x].im=Im(IU[z][XEL][y][x]);
		IMatrix2[y][z][x].re=Re(IU[z][YEL][y][x]);
		IMatrix2[y][z][x].im=Im(IU[z][YEL][y][x]);
		IMatrix3[y][z][x].re=Re(IU[z][ZEL][y][x]);
		IMatrix3[y][z][x].im=Im(IU[z][ZEL][y][x]);
	      }
	  }
      } 

	      sprintf(filename, "data(velocity)2=%d", n+1);

	      /* create File data.h5 for three data sets*/
	      file_id1 = H5Fcreate (filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

	      /* create the dataspace*/
	      fid1 = H5Screate_simple(3, fdim, NULL);
	
	      /* create datasets with name u v w*/
	      dataset_u = H5Dcreate(file_id1, "/data_u", complex_id, fid1, H5P_DEFAULT);
	      dataset_v = H5Dcreate(file_id1, "/data_v", complex_id, fid1, H5P_DEFAULT);
	      dataset_w = H5Dcreate(file_id1, "/data_w", complex_id, fid1, H5P_DEFAULT);
	      dataset_iu = H5Dcreate(file_id1, "/data_iu", complex_id, fid1, H5P_DEFAULT);
	      dataset_iv = H5Dcreate(file_id1, "/data_iv", complex_id, fid1, H5P_DEFAULT);
	      dataset_iw = H5Dcreate(file_id1, "/data_iw", complex_id, fid1, H5P_DEFAULT);

	      /* write data to corresponding datasets*/
	      ret = H5Dwrite(dataset_u,complex_id, H5S_ALL,fid1 , H5P_DEFAULT, Matrix1);
	      ret = H5Dwrite(dataset_v,complex_id, H5S_ALL,fid1 , H5P_DEFAULT, Matrix2);
	      ret = H5Dwrite(dataset_w,complex_id, H5S_ALL,fid1 , H5P_DEFAULT, Matrix3);
	      ret = H5Dwrite(dataset_iu,complex_id, H5S_ALL,fid1 , H5P_DEFAULT, IMatrix1);
	      ret = H5Dwrite(dataset_iv,complex_id, H5S_ALL,fid1 , H5P_DEFAULT, IMatrix2);
	      ret = H5Dwrite(dataset_iw,complex_id, H5S_ALL,fid1 , H5P_DEFAULT, IMatrix3);


	      /* close datasets and file */
	      ret=H5Dclose(dataset_u);
	      ret=H5Dclose(dataset_v);
	      ret=H5Dclose(dataset_w);
	      ret=H5Dclose(dataset_iu);
	      ret=H5Dclose(dataset_iv);
	      ret=H5Dclose(dataset_iw);

	      ret=H5Sclose(fid1);
	      ret=H5Fclose(file_id1);
	      /* 
	      for (i = 0; i < dimQ; ++i)
		{
		  for (z= 0; z < Nz; ++z)
		    {
		      for (x = 0; x < Nx/2; ++x)
			{
			  printf(" IC[%d][ALPHA][%d][%d]=%f+%fi\n", z, i, x,Re(IC[z][ALPHA][i][x]), Im(IC[z][ALPHA][i][x]));
			}
		    }
		}
	      for (i = 0; i < dimR; ++i)
		{
		  for (z= 0; z < Nz; ++z)
		    {
		      for (x = 0; x < Nx/2; ++x)
			{
			  printf(" IC[%d][BETA][%d][%d]=%f+%fi\n", z, i, x,Re(IC[z][BETA][i][x]), Im(IC[z][BETA][i][x]));
			}
		    }
		    }*/
  return(EXIT_SUCCESS);
}

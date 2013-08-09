#include "main.h"
#include "hdf5.h"
#include <stdio.h>
#include <string.h>
#include "minChnl.h"
#include "mvOps.h"
#include "arrays.h"
void restart(int restart_flag)
{

    extern int qpts, dimR, dimQ, Nx, Nz;
    extern mcomplex ****U, ****C;
    extern mcomplex ****IU, ****IC;

    int x, y, z, i;
    hid_t file_id1, dataset_u, dataset_v, dataset_w;    /* file identifier */
    hid_t dataset_iu, dataset_iv, dataset_iw;
    hid_t complex_id;
    herr_t ret;
    char filename[50];

    /* define compound datatype for the complex number */
    typedef struct {
        double re;              /*real part */
        double im;              /*imaginary part */
    } complex_t;

    complex_id = H5Tcreate(H5T_COMPOUND, sizeof(complex_t));
    H5Tinsert(complex_id, "real", HOFFSET(complex_t, re),
              H5T_NATIVE_DOUBLE);
    H5Tinsert(complex_id, "imaginary", HOFFSET(complex_t, im),
              H5T_NATIVE_DOUBLE);

    /* define some temporal matrix to store the data to the hdf file */
    complex_t Matrix1[qpts][Nz][Nx / 2];
    complex_t Matrix2[qpts][Nz][Nx / 2];
    complex_t Matrix3[qpts][Nz][Nx / 2];

    complex_t IMatrix1[qpts][Nz][Nx / 2];
    complex_t IMatrix2[qpts][Nz][Nx / 2];
    complex_t IMatrix3[qpts][Nz][Nx / 2];       // find the restart file 


    sprintf(filename, "data_t=%d", restart_flag);

    // open the file and dataset 
    file_id1 = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    dataset_u = H5Dopen1(file_id1, "/data_u");
    dataset_v = H5Dopen1(file_id1, "/data_v");
    dataset_w = H5Dopen1(file_id1, "/data_w");

    dataset_iu = H5Dopen1(file_id1, "/data_iu");
    dataset_iv = H5Dopen1(file_id1, "/data_iv");
    dataset_iw = H5Dopen1(file_id1, "/data_iw");

    ret =
        H5Dread(dataset_u, complex_id, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                Matrix1);
    ret =
        H5Dread(dataset_v, complex_id, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                Matrix2);
    ret =
        H5Dread(dataset_w, complex_id, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                Matrix3);
    ret =
        H5Dread(dataset_iu, complex_id, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                IMatrix1);
    ret =
        H5Dread(dataset_iv, complex_id, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                IMatrix2);
    ret =
        H5Dread(dataset_iw, complex_id, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                IMatrix3);

    ret = H5Dclose(dataset_u);
    ret = H5Dclose(dataset_v);
    ret = H5Dclose(dataset_w);
    ret = H5Dclose(dataset_iu);
    ret = H5Dclose(dataset_iv);
    ret = H5Dclose(dataset_iw);
    ret = H5Fclose(file_id1);

    for (y = 0; y < qpts; y++) {
        for (z = 0; z < Nz; ++z) {
            for (x = 0; x < Nx / 2; ++x) {
                Re(U[z][XEL][y][x]) = Matrix1[y][z][x].re;
                Im(U[z][XEL][y][x]) = Matrix1[y][z][x].im;
                Re(U[z][YEL][y][x]) = Matrix2[y][z][x].re;
                Im(U[z][YEL][y][x]) = Matrix2[y][z][x].im;
                Re(U[z][ZEL][y][x]) = Matrix3[y][z][x].re;
                Im(U[z][ZEL][y][x]) = Matrix3[y][z][x].im;

                Re(IU[z][XEL][y][x]) = IMatrix1[y][z][x].re;
                Im(IU[z][XEL][y][x]) = IMatrix1[y][z][x].im;
                Re(IU[z][YEL][y][x]) = IMatrix2[y][z][x].re;
                Im(IU[z][YEL][y][x]) = IMatrix2[y][z][x].im;
                Re(IU[z][ZEL][y][x]) = IMatrix3[y][z][x].re;
                Im(IU[z][ZEL][y][x]) = IMatrix3[y][z][x].im;
            }
        }
    }

    /* compute inital coefficients in expansion  given the value of ux hat, uy hat, uz hat. */
    initAlphaBeta();

    /* compute the boundary condition for previous stage or time step, result is stored in Uxb and Uzb */
    if (increBoundary() != NO_ERR) {
        printf("increBoundary failure\n");
    }

    /* compute initial coefficients for incremental state solver */

    incre_initAlphaBeta();
    printf(" restart computed coefficients!\n");
    for (i = 0; i < dimQ; ++i) {
        for (z = 0; z < Nz; ++z) {
            for (x = 0; x < Nx / 2; ++x) {
                printf(" IC[%d][ALPHA][%d][%d]=%f+%fi\n", z, i, x,
                       Re(IC[z][ALPHA][i][x]), Im(IC[z][ALPHA][i][x]));
            }
        }
    }
    for (i = 0; i < dimR; ++i) {
        for (z = 0; z < Nz; ++z) {
            for (x = 0; x < Nx / 2; ++x) {
                printf(" IC[%d][BETA][%d][%d]=%f+%fi\n", z, i, x,
                       Re(IC[z][BETA][i][x]), Im(IC[z][BETA][i][x]));
            }
        }
    }
    printf(" finish restart computing coefficients!\n\n");

}

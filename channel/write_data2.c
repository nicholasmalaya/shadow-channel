#include "main.h"
#include "hdf5.h"
#include <stdio.h>
#include <string.h>
#include "minChnl.h"
#include "mvOps.h"
#include "arrays.h"
int write_data2(int n)
{
    hid_t file_id1, dataset_a, dataset_b;       /* file identifier */
    hid_t fid1, dataset_ia, dataset_ib;
    hid_t fid2, dataset_Nx, dataset_Ny, dataset_Nz, dataset_dt, dataset_Re,
        dataset_mpg;
    hid_t complex_id;
    hsize_t fdim[] = { dimR, Nz, Nx / 2 };
    hsize_t fdim2[] = { 1 };
    herr_t ret;
    char filename[50];

    extern mcomplex ****C, ****IC;
    extern int qpts, Nx, Nz, dimR, dimQ;
    extern double dt, re, mpg;
    int x, y, z;
    int Ny = qpts * 2 / 3;
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
    complex_t Matrix1[dimR][Nz][Nx / 2];
    complex_t Matrix2[dimR][Nz][Nx / 2];

    complex_t IMatrix1[dimR][Nz][Nx / 2];
    complex_t IMatrix2[dimR][Nz][Nx / 2];

    for (y = 0; y < dimR; y++) {
        for (z = 0; z < Nz; ++z) {
            for (x = 0; x < Nx / 2; ++x) {
                Matrix1[y][z][x].re = Re(C[z][ALPHA][y][x]);    /* storing solved solutions to Matrix */
                Matrix1[y][z][x].im = Im(C[z][ALPHA][y][x]);
                Matrix2[y][z][x].re = Re(C[z][BETA][y][x]);
                Matrix2[y][z][x].im = Im(C[z][BETA][y][x]);

                IMatrix1[y][z][x].re = Re(IC[z][ALPHA][y][x]);  /* storing solved solutions to Matrix */
                IMatrix1[y][z][x].im = Im(IC[z][ALPHA][y][x]);
                IMatrix2[y][z][x].re = Re(IC[z][BETA][y][x]);
                IMatrix2[y][z][x].im = Im(IC[z][BETA][y][x]);
            }
        }
    }

    sprintf(filename, "data_t=%d.h5", n + 1);

    /* create File data.h5 for three data sets */
    file_id1 =
        H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    /* create the dataspace */
    fid1 = H5Screate_simple(3, fdim, NULL);
    fid2 = H5Screate_simple(1, fdim2, NULL);

    /* create datasets with name u v w */
    dataset_a =
        H5Dcreate1(file_id1, "/data_alpha", complex_id, fid1, H5P_DEFAULT);
    dataset_b =
        H5Dcreate1(file_id1, "/data_beta", complex_id, fid1, H5P_DEFAULT);

    dataset_ia =
        H5Dcreate1(file_id1, "/data_ialpha", complex_id, fid1,
                   H5P_DEFAULT);
    dataset_ib =
        H5Dcreate1(file_id1, "/data_ibeta", complex_id, fid1, H5P_DEFAULT);


    dataset_Nx =
        H5Dcreate1(file_id1, "data_Nx", H5T_STD_I32LE, fid2, H5P_DEFAULT);
    dataset_Ny =
        H5Dcreate1(file_id1, "data_Ny", H5T_STD_I32LE, fid2, H5P_DEFAULT);
    dataset_Nz =
        H5Dcreate1(file_id1, "data_Nz", H5T_STD_I32LE, fid2, H5P_DEFAULT);

    dataset_dt =
        H5Dcreate1(file_id1, "data_dt", H5T_IEEE_F64LE, fid2, H5P_DEFAULT);
    dataset_Re =
        H5Dcreate1(file_id1, "data_Re", H5T_IEEE_F64LE, fid2, H5P_DEFAULT);
    dataset_mpg =
        H5Dcreate1(file_id1, "data_mpg", H5T_IEEE_F64LE, fid2,
                   H5P_DEFAULT);

    /* write data to corresponding datasets */
    ret =
        H5Dwrite(dataset_a, complex_id, H5S_ALL, fid1, H5P_DEFAULT,
                 Matrix1);
    ret =
        H5Dwrite(dataset_b, complex_id, H5S_ALL, fid1, H5P_DEFAULT,
                 Matrix2);
    ret =
        H5Dwrite(dataset_ia, complex_id, H5S_ALL, fid1, H5P_DEFAULT,
                 IMatrix1);
    ret =
        H5Dwrite(dataset_ib, complex_id, H5S_ALL, fid1, H5P_DEFAULT,
                 IMatrix2);

    ret =
        H5Dwrite(dataset_Nx, H5T_NATIVE_INT, H5S_ALL, fid2, H5P_DEFAULT,
                 &Nx);
    ret =
        H5Dwrite(dataset_Ny, H5T_NATIVE_INT, H5S_ALL, fid2, H5P_DEFAULT,
                 &Ny);
    ret =
        H5Dwrite(dataset_Nz, H5T_NATIVE_INT, H5S_ALL, fid2, H5P_DEFAULT,
                 &Nz);
    ret =
        H5Dwrite(dataset_dt, H5T_IEEE_F64LE, H5S_ALL, fid2, H5P_DEFAULT,
                 &dt);
    ret =
        H5Dwrite(dataset_Re, H5T_IEEE_F64LE, H5S_ALL, fid2, H5P_DEFAULT,
                 &re);
    ret =
        H5Dwrite(dataset_mpg, H5T_IEEE_F64LE, H5S_ALL, fid2, H5P_DEFAULT,
                 &mpg);

    /* close datasets and file */
    ret = H5Dclose(dataset_a);
    ret = H5Dclose(dataset_b);
    ret = H5Dclose(dataset_ia);
    ret = H5Dclose(dataset_ib);
    ret = H5Dclose(dataset_Nx);
    ret = H5Dclose(dataset_Ny);
    ret = H5Dclose(dataset_Nz);
    ret = H5Dclose(dataset_dt);
    ret = H5Dclose(dataset_Re);
    ret = H5Dclose(dataset_mpg);

    ret = H5Sclose(fid1);
    ret = H5Sclose(fid2);
    ret = H5Fclose(file_id1);

    return (EXIT_SUCCESS);
}

void print_result(int n, int dctr, int flag)
{
    extern int dimQ, dimR, Nx, Nz, qpts;
    extern mcomplex ****U, ****C;
    extern mcomplex ****IU, ****IC;
    extern mcomplex **Uxb, **Uzb;
    int i, z, x, y;
    if (n >= 100) {
        printf("results comparison %d,%d, %d\n", n, dctr, flag);
        for (i = 0; i < dimQ; ++i) {
            for (z = 0; z < Nz; ++z) {
                for (x = 0; x < Nx / 2; ++x) {
                    printf(" C[%d][ALPHA][%d][%d]=%25.16e+i%25.16e\n", z,
                           i, x, Re(C[z][ALPHA][i][x]),
                           Im(C[z][ALPHA][i][x]));
                    printf(" C[%d][BETA][%d][%d]=%25.16e+%25.16ei\n", z, i,
                           x, Re(C[z][BETA][i][x]), Im(C[z][BETA][i][x]));
                }
            }
        }
        printf("results comparison %d, %d, %d\n", n, dctr, flag);
        for (i = 0; i < dimR; ++i) {
            for (z = 0; z < Nz; ++z) {
                for (x = 0; x < Nx / 2; ++x) {
                    printf(" IC[%d][ALPHA][%d][%d]=%25.16e+%25.16ei\n", z,
                           i, x, Re(IC[z][ALPHA][i][x]),
                           Im(IC[z][ALPHA][i][x]));
                    printf(" IC[%d][BETA][%d][%d]=%25.16e+%25.16ei\n", z,
                           i, x, Re(IC[z][BETA][i][x]),
                           Im(IC[z][BETA][i][x]));
                }
            }
        }
        printf("results comparison %d,%d,  %d\n", n, dctr, flag);
        for (i = 0; i < 5; ++i) {
            for (z = 0; z < Nz; ++z) {
                for (x = 0; x < Nx / 2; ++x) {
                    for (y = 0; y < qpts; y++) {
                        printf(" U[%d][%d][%d][%d]=%25.16e+%25.16ei\n", z,
                               i, y, x, Re(U[z][i][y][x]),
                               Im(U[z][i][y][x]));
                        printf(" IU[%d][%d][%d][%d]=%25.16e+%25.16ei\n", z,
                               i, y, x, Re(IU[z][i][y][x]),
                               Im(IU[z][i][y][x]));
                    }
                }
            }
        }
        printf("results comparison %d,%d,  %d\n", n, dctr, flag);
        for (z = 0; z < Nz; z++) {
            for (x = 0; x < Nx / 2; x++) {
                printf
                    ("Uxb[%d][%d]=%25.16e+i%25.16e, Uzb=%25.16e+i%25.16e\n",
                     z, x, Re(Uxb[z][x]), Im(Uxb[z][x]), Re(Uzb[z][x]),
                     Im(Uzb[z][x]));
            }
        }
    }
}

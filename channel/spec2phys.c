#include <math.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include "fftw.h"
#include "rfftw.h"
#include "main.h"
#include "minChnl.h"

void spec2phys(mcomplex * C_ptr, double * flow)
{
	/* External Variables */
	extern int Nx, Nz, qpts, dimR;
	extern mcomplex ****C, ****U;
	extern fftw_plan pf1;
	extern rfftwnd_plan pr1;
	extern fftw_complex ***CT;
	int x, y, z, i, idx;
    int XFLOW, YFLOW, ZFLOW;
	fftw_real *RT;		/* real to complex transform */
	fftw_complex *fout = NULL;
	fftw_real *rout = NULL;
	extern mcomplex ****U, ****IU;

    /* copy C_ptr to C */
    memcpy(C[0][0][0], C_ptr, Nz * 2 * dimR * (Nx / 2) * sizeof(mcomplex));

    /* Compute U from C */
    initAlphaBeta2();

    /* Compute flow from U */
	idx = (3 * Nz / 2) * (3 * Nx / 2 + 2);
	RT = (fftw_real *) CT[0][0];

    XFLOW = 0 * (3 * Nx / 2) * qpts * (3 * Nz / 2);
    YFLOW = 1 * (3 * Nx / 2) * qpts * (3 * Nz / 2);
    ZFLOW = 2 * (3 * Nx / 2) * qpts * (3 * Nz / 2);

	for (y = 0; y < qpts; ++y) {
		memset(CT[0][0], 0,
		       MAXTT * (3 * Nz / 2) * (3 * Nx / 4 +
					       1) * sizeof(fftw_complex));

		for (i = 0; i < 3; ++i) {
			for (z = 0; z < Nz / 2; ++z) {
				memcpy(CT[i][z], U[z][i][y],
				       (Nx / 2) * sizeof(fftw_complex));
				memcpy(CT[i + 3][z], IU[z][i][y],
				       (Nx / 2) * sizeof(fftw_complex));
			}
			for (z = Nz / 2 + 1; z < Nz; ++z) {
				memcpy(CT[i][z + Nz / 2], U[z][i][y],
				       (Nx / 2) * sizeof(fftw_complex));
				memcpy(CT[i + 3][z + Nz / 2],
				       IU[z][i][y],
				       (Nx / 2) * sizeof(fftw_complex));
			}
		}
		for (z = 0; z < Nz / 2; ++z) {
			memcpy(CT[6][z], U[z][DXEL][y],
			       (Nx / 2) * sizeof(fftw_complex));
		}
		for (z = Nz / 2 + 1; z < Nz; ++z) {
			memcpy(CT[6][z + Nz / 2], U[z][DXEL][y],
			       (Nx / 2) * sizeof(fftw_complex));
		}
		for (i = 0; i < 3; ++i) {
			/* Each column of CT[i] */
			fftw(pf1, Nx / 2, CT[i][0], 3 * Nx / 4 + 1, 1,
			     fout, -1, -1);

			/* Each row of CT[i] */
			rfftwnd_complex_to_real(pr1, 3 * Nz / 2,
						CT[i][0], 1,
						3 * Nx / 4 + 1, rout, -1, -1);
		}

		for (z = 0; z < (3 * Nz / 2); ++z) {
			for (x = 0; x < 3 * Nx / 2; ++x) {
                i = (x * qpts + y) * (3 * Nz / 2) + z;
				flow[XFLOW + i] =
				   RT[XEL * idx + (z * (3 * Nx / 2 + 2) + x)];
				flow[YFLOW + i] =
				   RT[YEL * idx + (z * (3 * Nx / 2 + 2) + x)];
				flow[ZFLOW + i] =
				   RT[ZEL * idx + (z * (3 * Nx / 2 + 2) + x)];
			}
		}
    }
}

#!/bin/py
#
# run the channel and create plots comparable to keefe (92)
#
import os
import sys
sys.path.append("..")
import channel 

# for plotting
import pylab as pl
import numpy as np

# reynolds number
channel.Re=2000

# invoke init
n_steps = 1
ru_steps = 0

restart_file = 'keefe_runup_stage_5.hd5'
assert os.path.exists(restart_file)
channel.init(n_steps, ru_steps, restart_file)

#
# postprocess
#
#
# void getsoln(int i_step, mcomplex ** MC_ptr,
#              int * Nz_ptr, int * Nvar_ptr, int * Ny_ptr, int * Nx_ptr);
# 
#     # us[0][y] = u1;
#     # us[1][y] = u2;
#     # us[2][y] = u3;
#     # us[3][y] = u1u2;
#     # us[4][y] = u1u3;
#     # us[5][y] = u2u3;
#     # us[12][y] = u1u1;
#     # us[13][y] = u2u2;
#     # us[14][y] = u3u3;
#     # us[18][y] = u1y;

# degree of Legendre polynomial
n = (channel.Ny - 1) * 3 / 2 + 1
y, w = numpy.polynomial.legendre.leggauss(n)

for i in range(0, n_steps + 1, 400):
    stats = channel.statistics(i)
    ax = pl.subplot(2,1,1)
    if i == 0:
        pl.plot(y 3./2 * (1 - y**2), ':k')
    pl.plot(y, stats[0], '.-')

    # calculate u_tau by finite 
    utau  = np.sqrt((1/channel.Re)*((stats[0,1] - stats[0,0]) / (y[1]-y[0])))
    retau = utau*channel.Re
    print("utau: ", utau)
    print("re_tau: ", retau)
    
    # ax.set_xscale('log')
    yp = (y[0]+1)*retau
    up = stats[0]/utau
    # pl.plot(yp, up, '.-')

    pl.subplot(2,1,2)
    # pl.plot(stats[0], np.sqrt(stats[13] - stats[1]**2)/utau, '.-')
    pl.plot(y, np.sqrt(stats[12] - stats[0]**2), '.-')


xl = pl.xlabel(r'$\mathrm{y/\delta}$')
pl.show()

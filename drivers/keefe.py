#!/bin/py
#
# run the channel and create plots comparable to keefe (92)
#
import os
import sys

sys.path.append("../channel/")
import channel 

# for plotting
import pylab as pl
import numpy as np

# try subprocess to 
# redirect all the screen outputs in C  to a buffer
#  which we can access in Python through channel.stdout
#
# http://docs.python.org/2/library/subprocess.html
#
# ------------------------------------------------------------------
# main function
# ------------------------------------------------------------------
#
# parameters
#
# number of pts. 
Nx=16
Ny=33
Nz=16

# box size
Lx=1.6
Lz=1.6

# time step
dt = .01

# run up (or spin up) time
ru_steps = 0

# steps
nsteps=1200

# restart flag: look for the largest steps in the currect directory
restart_flag = 0
for fname in os.listdir('.'):
    if fname.startswith('data_t=') and fname.endswith('.h5'):
        n = int(fname[len('data_t='):-len('.h5')])
        restart_flag = max(n, restart_flag)

restart_flag = 10000
print('restart from ', restart_flag)

# flux (mean pressure gradient)
flux = 2

# reynolds number
Re=1190

# number of time chunks
nchunk=1

# invoke init
channel.init(Nx,Ny,Nz,Lx,Lz,Re,flux,dt,ru_steps,nchunk,nsteps,restart_flag)

#
# postprocess
#
#
# void getsoln(int i_step, mcomplex ** MC_ptr,
#              int * Nz_ptr, int * Nvar_ptr, int * Ny_ptr, int * Nx_ptr);
# 
#     # us[0][y] = Qy[y];
#     # us[1+0][y] = u1;
#     # us[1+1][y] = u2;
#     # us[1+2][y] = u3;
#     # us[1+3][y] = u1u2;
#     # us[1+4][y] = u1u3;
#     # us[1+5][y] = u2u3;
#     # us[1+12][y] = u1u1;
#     # us[1+13][y] = u2u2;
#     # us[1+14][y] = u3u3;
#     # us[1+18][y] = u1y;

for i in range(0, nsteps + 1, 400):
    C = channel.getsoln(i)
    stats = channel.statistics(C)
    pl.subplot(2,1,1)
    if i == 0:
        pl.plot(stats[0], 3./2 * (1 - stats[0]**2), ':k')
    pl.plot(stats[0], stats[1], '.-')
    pl.subplot(2,1,2)
    pl.plot(stats[0], sqrt(stats[13] - stats[1]**2), '.-')

pl.show()

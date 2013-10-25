#!/bin/py
#
# run the channel and create plots comparable to keefe (92)
#
import sys
sys.path.append("../channel/")
import channel 

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

# time step and number of steps
dt = .0005
nsteps=10

# restart flag: 0 for no restart
restart_flag=0

# mpg (mean pressure gradient)
mpg=2

# reynolds number
Re=32

# run up (or spin up) time
ru_steps = 0

# number of time chunks
nchunk=1

# invoke init
channel.init(Nx,Ny,Nz,Lx,Lz,Re,mpg,dt,ru_steps,nchunk,nsteps,restart_flag)

#
# postprocess
#
#
# void getsoln(int i_step, mcomplex ** MC_ptr,
#              int * Nz_ptr, int * Nvar_ptr, int * Ny_ptr, int * Nx_ptr);
# 
for i in range(nsteps + 1):
    C = channel.getsoln(i)
    stats = channel.statistics(C)
    # us[0][y] = Qy[y];
    # us[1+0][y] = u1;
    # us[1+1][y] = u2;
    # us[1+2][y] = u3;
    # us[1+3][y] = u1u2;
    # us[1+4][y] = u1u3;
    # us[1+5][y] = u2u3;
    # us[1+12][y] = u1u1;
    # us[1+13][y] = u2u2;
    # us[1+14][y] = u3u3;
    # us[1+18][y] = u1y;
    plot(stats[0], stats[1])

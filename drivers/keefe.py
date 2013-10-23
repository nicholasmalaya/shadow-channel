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
nsteps=101

# restart flag: 0 for no restart
restart_flag=0

# mpg (mean pressure gradient)
mpg=2

# reynolds number
Re=32

# run up (or spin up) time
ru_steps = 1

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
U = channel.getsoln(0)

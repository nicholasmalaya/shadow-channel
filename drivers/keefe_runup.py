#!/bin/py
#
# run the channel and create plots comparable to keefe (92)
#
import sys
sys.path.append("..")
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
channel.Nx=16
channel.Ny=33
channel.Nz=16

# box size
channel.Lx=1.6
channel.Lz=1.6

# time step and number of steps
channel.dt = .01

# flux (mean pressure gradient)
channel.meanU = 1

# ============================================
# STAGES -- reynolds number
Re = [10000, 5000, 3000, 2500, 2000]
restart = None

for i in range(len(Re)):
    channel.Re = Re[i]
    
    # STAGE i -- time steps
    ru_steps = 4999
    n_steps = 1
    
    # STAGE i -- invoke init
    channel.init(n_steps, ru_steps, restart=restart)

    filename = 'keefe_runup_stage_{0}'.format(i+1)
    channel.save_solution(filename, channel.get_solution(0))
    channel.destroy()

    restart = filename


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
# STAGE 1 -- reynolds number
channel.Re=10000

# STAGE 1 -- time steps
ru_steps = 4999
n_steps = 1

# STAGE 1 -- invoke init
channel.init(n_steps, ru_steps, restart=None)
channel.save_solution('keefe_runup_stage_1', channel.get_solution(0))
channel.destroy()

# ============================================
# STAGE 2 -- reynolds number
channel.Re=2000

# STAGE 1 -- time steps
ru_steps = 1999
n_steps = 1

# STAGE 1 -- invoke init
channel.init(n_steps, ru_steps, restart='keefe_runup_stage_1')
channel.save_solution('keefe_runup_stage_2', channel.get_solution(0))
channel.destroy()


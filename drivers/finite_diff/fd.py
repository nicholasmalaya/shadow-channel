#!/bin/py
#
# run the channel and create plots comparable to keefe (92)
#
import os
import sys
sys.path.append("../..")
import channel 

# for plotting
from pylab import *
from numpy import *

# reynolds number
channel.Re=int(sys.argv[1])

# invoke init
n_steps = 10000

restart_file = '../keefe_runup_stage_5_qiqi.hd5'

for i_stage in range(1000000):
    next_restart = 'Re{0:04d}_{1:06d}.hd5'.format(channel.Re, i_stage)

    if not os.path.exists(next_restart):
        print('restarting from ', restart_file, ' writing to ', next_restart)
        sys.stdout.flush()

        channel.init(n_steps, 0, restart_file)
        
        u_mean, u2_mean = [], []
        for i in range(0, n_steps + 1, 400):
            U = channel.spec2phys(i)
            u_mean.append(U.mean(1).mean(2))
            u2_mean.append((U**2).mean(1).mean(2))
        
        u_mean, u2_mean = array(u_mean).mean(0), array(u2_mean).mean(0)
        
        savetxt('Re{0:04d}_{1:06d}.txt'.format(channel.Re, i_stage),
                vstack([u_mean, u2_mean]).T)
        channel.save_solution(next_restart, channel.get_solution(n_steps))

    restart_file = next_restart

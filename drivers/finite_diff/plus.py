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
channel.Re=2100

# invoke init
n_steps = 10000
ru_steps = 1000

restart_file = '../keefe_runup_stage_5_qiqi.hd5'
assert os.path.exists(restart_file)
channel.init(n_steps, ru_steps, restart_file)

y, w = channel.quad()

u_mean, u2_mean = [], []
for i in range(0, n_steps + 1, 400):
    U = channel.spec2phys(i)
    u_mean.append(U.mean(1).mean(2))
    u2_mean.append((U**2).mean(1).mean(2))

u_mean, u2_mean = array(u_mean).mean(0), array(u2_mean).mean(0)

# plot(y, u_mean.T)
# plot(y, u2_mean.T)

savetxt('Re2100.txt', vstack([u_mean, u2_mean]).T)

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
channel.init(100, 0, None)

C = channel.get_solution(0) * 0
channel.primal(C)

y, w = channel.quad()

for i in [0, 2, 10, 40, 100]:
    u = channel.spec2phys(i)
    pl.plot(y, u[0].mean(0).mean(1), '.-')


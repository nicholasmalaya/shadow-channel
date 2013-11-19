import sys
sys.path.append("..")
import channel
from numpy import *

# TANGENT WITH FORCING
Re = 200
dRe = 1E-6
ru_steps = 0
n_steps = int(sys.argv[1])
i_restart = int(sys.argv[2])
restart = 'keefe_runup_stage_5.hd5' if i_restart else None
channel.dt = 0.01
channel.meanU = 1.0
channel.Re = Re
channel.init(n_steps, ru_steps, restart=restart)

# tangent 
IC = zeros_like(channel.get_solution(0))
channel.tangent(0, n_steps, IC, 1)

from pylab import *
y, w = channel.quad()
IU = channel.spec2phys(IC)
plot(y, IU[0].mean(0).mean(1))
savefig('testTanManu{0:06d}_restart{1:1d}'.format(n_steps, restart is not None))

# channel.destroy()

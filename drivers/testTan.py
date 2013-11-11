import sys
sys.path.append("..")
import channel
from numpy import *

# UNFORCED TANGENT

n_steps = 10

channel.init(n_steps, ru_steps=100, restart=None)
# initial condition

C = channel.get_solution(0, copy=True)
C1 = channel.get_solution(n_steps, copy=True)

# perturbation
IC0 = 0.000001*(random.rand(C.shape[0],C.shape[1],C.shape[2],C.shape[3]) + 1j * random.rand(C.shape[0],C.shape[1],C.shape[2],C.shape[3]))
print(linalg.norm(IC0))
# run perturbed primal
channel.primal(C + IC0)

C1p = channel.get_solution(n_steps, copy=True)
print(linalg.norm(C1p-C1))
# run tangent
#ICh = zeros(C.shape,complex)
#channel.tangent(0, n_steps, ICh,0)
#print linalg.norm(ICh)

IC = copy(IC0)
channel.tangent(0, n_steps, IC,0)
print(linalg.norm(IC))
print(linalg.norm(C1p-C1-IC))

channel.destroy()

'''
# TANGENT WITH FORCING
Re = 2000
dRe = 1E-6
ru_steps = 0
n_steps = 10
restart = 'keefe_runup_stage_5.hd5'
# restart = None
channel.dt = 0.01
channel.meanU = 1.0
# positive perturbation
channel.Re = Re + dRe
channel.init(n_steps, ru_steps, restart=restart)
Cp = channel.get_solution(n_steps,copy=True)
channel.destroy()

# negative perturbation
channel.Re = Re - dRe
channel.init(n_steps, ru_steps, restart=restart)
Cm = channel.get_solution(n_steps,copy=True)
channel.destroy()

# tangent 
channel.Re = Re
channel.init(n_steps, ru_steps, restart=restart)
IC = zeros(Cm.shape,complex)
channel.tangent(0, n_steps, IC,1)

dC = (1./(2*dRe))*(Cp - Cm)

print linalg.norm(dC)
print linalg.norm(IC)
print linalg.norm(dC-IC)
# channel.destroy()

import numpy
n = int((channel.Ny - 1) * 3 / 2 + 1)
y, w = numpy.polynomial.legendre.leggauss(n)

ds = channel.statistics(dC)
incre_s = channel.statistics(IC)
from pylab import *
plot(y, ds[0])
plot(y, incre_s[0])
s = channel.statistics(Cp)
figure(); plot(y, s[0])
show()
'''

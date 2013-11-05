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
print linalg.norm(IC0)
# run perturbed primal
channel.primal(C + IC0)

C1p = channel.get_solution(n_steps, copy=True)
print linalg.norm(C1p-C1)
# run tangent
#ICh = zeros(C.shape,complex)
#channel.tangent(0, n_steps, ICh,0)
#print linalg.norm(ICh)

IC = copy(IC0)
channel.tangent(0, n_steps, IC,0)
print linalg.norm(IC)
print linalg.norm(C1p-C1-IC)

channel.destroy()
'''
# TANGENT WITH FORCING
Re = 5000
dRe = 1E-5
ru_steps = 0
n_steps = 1
# positive perturbation
channel.Re = Re + dRe
channel.init(n_steps, ru_steps, restart=None)
Cp = channel.get_solution(0)
channel.destroy()

# negative perturbation
channel.Re = Re - dRe
channel.init(n_steps, ru_steps, restart=None)
Cm = channel.get_solution(0)
channel.destroy()

# tangent 
channel.Re = Re
channel.init(n_steps, ru_steps, restart=None)
IC = zeros(Cm.shape,complex)
channel.tangent(0, n_steps, IC,1)
channel.destroy()
'''

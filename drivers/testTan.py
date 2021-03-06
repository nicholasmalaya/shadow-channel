import sys
sys.path.append("..")
import channel
from numpy import *

'''
# UNFORCED TANGENT

channel.Re = 20000
n_steps = 1000

channel.init(n_steps, ru_steps=1000, restart=None)
# initial condition

C = channel.get_solution(0, copy=True)
C1 = channel.get_solution(n_steps, copy=True)

U = channel.spec2phys(C)
U1 = channel.spec2phys(C1)

n = int((channel.Ny - 1) * 3 / 2 + 1)
y, w = numpy.polynomial.legendre.leggauss(n)

plot(y, U[0].mean(0).mean(1))
plot(y, U1[0].mean(0).mean(1))

print("Primal")
for i in range(0,n_steps):
    Ctmp = channel.get_solution(i, copy=True)
    print(linalg.norm(Ctmp))

# perturbation
random.seed(0)
IC0 = 0.000001*(random.rand(C.shape[0],C.shape[1],C.shape[2],C.shape[3]) + 1j * random.rand(C.shape[0],C.shape[1],C.shape[2],C.shape[3]))
print("Perturbation")
print(linalg.norm(IC0))

channel.primal(C + IC0)
#channel.save_solution('perturbed_C.hd5', C + IC0)
#channel.destroy()

# run perturbed primal
# channel.init(n_steps, ru_steps=0, restart='perturbed_C.hd5')

print("Perturbed Primal")
for i in range(0,n_steps):
    Ctmp = channel.get_solution(i, copy=True)
    print(linalg.norm(Ctmp))


C1p = channel.get_solution(n_steps, copy=True)
print("Norm of difference")
print(linalg.norm(C1p-C1))
# run tangent
#ICh = zeros(C.shape,complex)
#channel.tangent(0, n_steps, ICh,0)
#print linalg.norm(ICh)

IC = copy(IC0)
channel.tangent(0, n_steps, IC,0)
print("Norm of tangent")
print(linalg.norm(IC))
print("Difference of Norms")
print(linalg.norm(C1p-C1-IC))

channel.destroy()

'''
# TANGENT WITH FORCING
Re = 2000
dRe = 1E-6
ru_steps = 0
n_steps = 1000
restart = 'keefe_runup_stage_5.hd5'
restart = None
channel.dt = 0.01
channel.meanU = 1.0
# positive perturbation
channel.Re = Re + dRe
channel.init(n_steps, ru_steps, restart=restart)
# channel.primal(channel.get_solution(0) * 0)
Cp = channel.get_solution(n_steps,copy=True)
channel.destroy()

# negative perturbation
channel.Re = Re - dRe
channel.init(n_steps, ru_steps, restart=restart)
# channel.primal(channel.get_solution(0) * 0)
Cm = channel.get_solution(n_steps,copy=True)
channel.destroy()

# tangent 
channel.Re = Re
channel.init(n_steps, ru_steps, restart=restart)
# channel.primal(channel.get_solution(0) * 0)
IC = zeros(Cm.shape,complex)
channel.tangent(0, n_steps, IC, 1)

dC = (1./(2*dRe))*(Cp - Cm)

print(linalg.norm(dC))
print(linalg.norm(IC))
print(linalg.norm(dC-IC))
# channel.destroy()
'''
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

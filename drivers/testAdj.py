import os
import sys
sys.path.append("..")
import channel
from numpy import *

restart_file = 'keefe_runup_stage_5_qiqi.hd5'
if not os.path.exists(restart_file):
    url = 'https://www.dropbox.com/s/8nqihkqf8sd4w8z/keefe_runup_stage_5.hd5'
    print('Please download\n{0}\nand name it\n{1}\n'.format(url, restart_file))

# Primal
n_steps = 100

channel.Re = 2000
channel.init(n_steps, ru_steps=0, restart=restart_file) # turbulent

# channel.Re = 100
# channel.init(n_steps, ru_steps=0, restart=None) # laminar
C = channel.get_solution(0)

# perturbation
IC = zeros(C.shape, complex)
IC[0,0,:,0] = random.normal(size=C.shape[2])
IC = random.normal(size=C.shape) + 0j
channel.tangent(0, 1, IC, 0)
IC0 = IC.copy()

# run tangent
print('tangent')
channel.tangent(1, n_steps-1, IC, 0)
IC_tmp = IC.copy()
channel.tangent(n_steps-1, n_steps, IC_tmp, 0)
print('done')

# linear functional
AC = zeros(C.shape, complex)
AC[0,0,:,0] = random.normal(size=C.shape[2])
AC = random.normal(size=C.shape) + 0j
channel.adjoint(n_steps, n_steps-1, AC, 0, strength=0.0)
AC0 = AC.copy()

# run adjoint
print('adjoint')
channel.adjoint(n_steps-1, 1, AC, 0, strength=0.0)
print('done')

y, w = channel.quad()

Au = channel.spec2phys(AC)
Au0 = channel.spec2phys(AC0)
Iu = channel.spec2phys(IC)
Iu0 = channel.spec2phys(IC0)

print((Au0 * Iu * w[:,newaxis]).sum())
print((Au * Iu0 * w[:,newaxis]).sum())


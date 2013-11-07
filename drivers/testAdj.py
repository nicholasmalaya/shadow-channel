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
n_steps = 10

# channel.Re = 2000
# channel.init(n_steps, ru_steps=0, restart=restart_file) # turbulent

channel.Re = 100
channel.init(n_steps, ru_steps=0, restart=None) # laminar
C = channel.get_solution(0)

# perturbation
IC = zeros(C.shape, complex)
IC[0,0,:,0] = random.normal(size=C.shape[2]) \
            + 1j * random.normal(size=C.shape[2])
# IC = random.normal(size=C.shape) + 1j * random.normal(size=C.shape)
IC0 = IC.copy()

# run tangent
channel.tangent(0, n_steps, IC, 0)

# linear functional
AC = zeros(C.shape, complex)
AC[0,0,:,0] = random.normal(size=C.shape[2]) \
            + 1j * random.normal(size=C.shape[2])
# AC = random.normal(size=C.shape) + 1j * random.normal(size=C.shape)
AC0 = AC.copy()

# run adjoint
channel.adjoint(n_steps, 0, AC, 0)

Au = channel.statistics(AC)[1]
Au0 = channel.statistics(AC0)[1]
Iu = channel.statistics(IC)[1]
Iu0 = channel.statistics(IC0)[1]


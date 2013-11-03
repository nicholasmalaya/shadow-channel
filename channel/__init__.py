import numpy as np
import channel.c_channel
#
# default parameters
#
# number of pts. 
Nx=16
Ny=33
Nz=16

# box size
Lx=1.6
Lz=1.6

# time step and number of steps
dt = .01

# mean velocity
meanU = 1.

# reynolds number
Re=5000

def save_solution(filename, C):
    if not filename.endswith('.hd5'):
        filename += '.hd5'
    assert C.shape == (c_channel.c_Nz(), 2, c_channel.c_dimR(),
                       c_channel.c_Nx() / 2)
    c_channel.c_save_solution(filename, C)

def read_solution(filename):
    if not filename.endswith('.hd5'):
        filename += '.hd5'
    shp = (c_channel.c_Nz(), 2, c_channel.c_dimR(), c_channel.c_Nx() / 2)
    C = np.zeros(shp, complex)
    c_channel.c_read_solution(filename, C)
    return C

def get_solution(i_step, copy=True):
    C = c_channel.c_getsoln(i_step)
    assert C.shape == (c_channel.c_Nz(), 2, c_channel.c_dimR(),
                       c_channel.c_Nx() / 2)
    if copy:
        C = C.copy()
    return C

def init(n_steps, ru_steps=0, restart=None):
    assert n_steps >= 0 and ru_steps >= 0
    c_channel.c_init(Nx,Ny,Nz,Lx,Lz,Re,meanU*2,dt,n_steps)
    if restart is None:
        c_channel.c_primal(ru_steps, np.zeros([0,0,0,0], complex))
    else:
        C = read_solution(restart)
        assert C.shape == (c_channel.c_Nz(), 2, c_channel.c_dimR(),
                           c_channel.c_Nx() / 2)
        c_channel.c_primal(ru_steps, C)

def destroy():
    c_channel.c_destroy()
    
def primal(C_init):
    assert C_init.shape == (c_channel.c_Nz(), 2, c_channel.c_dimR(),
                            c_channel.c_Nx() / 2)
    c_channel.c_primal(0, C_init)

def statistics(solution):
    if isinstance(solution, str):
        solution = read_solution(solution)
    if isinstance(solution, int):
        solution = get_solution(solution, copy=False)

    assert solution.dtype == complex
    assert solution.shape == (c_channel.c_Nz(), 2, c_channel.c_dimR(),
                              c_channel.c_Nx() / 2)

    stats = np.zeros([21, c_channel.c_qpts()], np.float64)
    c_channel.c_statistics(solution, stats)
    return stats


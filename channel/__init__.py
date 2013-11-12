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

def quad():
    y, w = c_channel.c_quad()
    return y.copy(), w.copy()

def save_solution(filename, C):
    '''
    Save the solution C to a .hd5 file.  Must be called after init()
    The solution must be a complex array of dimension (Nz, 2, dimR, Nx/2)
    The .hd5 extension can be appended automatically.  No return value
    '''
    if not filename.endswith('.hd5'):
        filename += '.hd5'
    assert C.shape == (c_channel.c_Nz(), 2, c_channel.c_dimR(),
                       c_channel.c_Nx() / 2)
    c_channel.c_save_solution(filename, C)

def read_solution(filename):
    '''
    Read a .hd5 file and return a solution.  Must be called after init()
    Returns a complex array of dimension (Nz, 2, dimR, Nx/2)
    The .hd5 extension can be appended automatically.
    '''
    if not filename.endswith('.hd5'):
        filename += '.hd5'
    shp = (c_channel.c_Nz(), 2, c_channel.c_dimR(), c_channel.c_Nx() / 2)
    C = np.zeros(shp, complex)
    c_channel.c_read_solution(filename, C)
    return C

def get_solution(i_step, copy=True):
    '''
    Get the primal solution at the i_step time step (not including run-up).
    Returns a complex array of dimension (Nz, 2, dimR, Nx/2)
    If copy is True (default), make a copy of the solution.
                               otherwise return a view of the solution (danger)
    '''
    C = c_channel.c_getsoln(i_step)
    assert C.shape == (c_channel.c_Nz(), 2, c_channel.c_dimR(),
                       c_channel.c_Nx() / 2)
    if copy:
        C = C.copy()
    return C

def init(n_steps, ru_steps=0, restart=None):
    '''
    Must be called before any other calls in this module.
    n_steps: the number of time steps to run and save
    ru_steps: the number of run up time steps (not including n_steps).
              These steps are not saved.
    restart: If None (default), start from a perturbed laminar solution.
             If a filename, read the .hd5 file as the restart file.
    '''
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
    '''
    Must be called before calling init again to initialize another channel
    flow simulation.  Cleans up things.
    '''
    c_channel.c_destroy()
    
def primal(C_init):
    '''
    Solve the primal equation starting from C_init for nsteps time steps,
    where nsteps is specified in calling init().
    After this function call, C_init stores the solution after
    nsteps time steps.
    '''
    assert C_init.shape == (c_channel.c_Nz(), 2, c_channel.c_dimR(),
                            c_channel.c_Nx() / 2)
    c_channel.c_primal(0, C_init)

def tangent(start_step, end_step, IC_init, inhomo):
    '''
    Solve the tangent equation starting from IC_init for
    end_step - start_step time steps, starting from the start_step time step.
    After this function call, IC_init stores the solution after
    nsteps time steps.
    '''
    assert IC_init.shape == (c_channel.c_Nz(), 2, c_channel.c_dimR(),
                             c_channel.c_Nx() / 2)
    assert 0 <= start_step <= end_step <= c_channel.c_nsteps()
    c_channel.c_tangent(start_step, end_step, IC_init, inhomo)

def adjoint(start_step, end_step, AC_init, inhomo, strength):
    assert AC_init.shape == (c_channel.c_Nz(), 2, c_channel.c_dimR(),
                             c_channel.c_Nx() / 2)
    assert 0 <= end_step <= start_step <= c_channel.c_nsteps()
    c_channel.c_adjoint(start_step, end_step, AC_init, inhomo, strength)

"""
def statistics(solution):
    if isinstance(solution, str):
        solution = read_solution(solution)
    if isinstance(solution, int):
        solution = get_solution(solution, copy=False)

    assert solution.dtype == complex
    assert solution.shape == (c_channel.c_Nz(), 2, c_channel.c_dimR(),
                              c_channel.c_Nx() / 2)

    stats = np.zeros([20, c_channel.c_qpts()], np.float64)
    c_channel.c_statistics(solution, stats)
    return stats
"""

def spec2phys(solution):
    if isinstance(solution, str):
        solution = read_solution(solution)
    if isinstance(solution, int):
        solution = get_solution(solution, copy=False)

    assert solution.dtype == complex
    assert solution.shape == (c_channel.c_Nz(), 2, c_channel.c_dimR(),
                              c_channel.c_Nx() / 2)

    flow = np.zeros([3, int(c_channel.c_Nx() * 3 / 2), c_channel.c_qpts(),
                     int(c_channel.c_Nz() * 3 / 2)], np.float64)
    c_channel.c_spec2phys(solution, flow)
    return flow


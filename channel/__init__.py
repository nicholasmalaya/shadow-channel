import numpy as np
import channel.c_channel
#
# global lock (per process)
#
__is_initialized__ = False
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
    c_channel.c_read_solution(str(filename), C)
    return C

def get_solution(i_step, copy=True):
    '''
    Get the primal solution at the i_step time step (not including run-up).
    Returns a complex array of dimension (Nz, 2, dimR, Nx/2)
    If copy is True (default), make a copy of the solution.
                               otherwise return a view of the solution (danger)
    '''
    C = c_channel.c_getsoln(int(i_step))
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
    assert not __is_initialized__
    __is_initialized__ = True

    assert n_steps >= 0 and ru_steps >= 0
    c_channel.c_init(int(Nx),int(Ny),int(Nz),
                     float(Lx),float(Lz),float(Re),
                     float(meanU*2),float(dt),int(n_steps))
    if restart is None:
        c_channel.c_primal(int(ru_steps), np.zeros([0,0,0,0], complex))
    else:
        C = read_solution(restart)
        assert C.shape == (c_channel.c_Nz(), 2, c_channel.c_dimR(),
                           c_channel.c_Nx() / 2)
        c_channel.c_primal(int(ru_steps), C)

def para_init(n_steps, C):
    '''
    Like init, but takes a complex array of dimension (Nz, 2, dimR, Nx/2) 
    as an input for the initial condition instead of a restart file.  
    Used in parallel LSS.
    '''
    assert n_steps >= 0 
    c_channel.c_init(int(Nx),int(Ny),int(Nz),
                     float(Lx),float(Lz),float(Re),
                     float(meanU*2),float(dt),int(n_steps))
    assert C.shape == (c_channel.c_Nz(), 2, c_channel.c_dimR(),
                           c_channel.c_Nx() / 2)
    c_channel.c_primal(0, C)

def destroy():
    '''
    Must be called before calling init again to initialize another channel
    flow simulation.  Cleans up things.
    '''
    assert __is_initialized__
    __is_initialized__ = False
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
    c_channel.c_tangent(int(start_step), int(end_step), IC_init, int(inhomo))

def adjoint(start_step, end_step, AC_init, inhomo, strength):
    assert AC_init.shape == (c_channel.c_Nz(), 2, c_channel.c_dimR(),
                             c_channel.c_Nx() / 2)
    assert 0 <= end_step <= start_step <= c_channel.c_nsteps()
    c_channel.c_adjoint(int(start_step), int(end_step), AC_init, int(inhomo),
                        float(strength))

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

#def ddt_project(i_step,v):
#    '''
#    returns the magnitude of the projection of the array v onto the primal 
#    time derivative dudt at time step i_step.  
#    Makes v orthogonal to dudt at time step i_step
#    '''
#
#    if i_step == 0:
#        dudt = (get_solution(i_step+1, copy=True) - get_solution(i_step, copy=True)) / dt
#    else:
#        dudt = (get_solution(i_step, copy=True) - get_solution(i_step-1, copy=True)) / dt
#    #print "Cp = ", np.linalg.norm(get_solution(i_step, copy=True)) 
#    #print "Cm = ",  np.linalg.norm(get_solution(i_step-1, copy=True))
#    #print "dCdt= ", np.linalg.norm(dudt)
#    dudt_p = spec2phys(dudt)
#    v_p = spec2phys(v)
#    print "dudt= ", np.linalg.norm(dudt_p)
#    print "v= ", np.linalg.norm(v_p)
#
#    y, w = quad()
#
#    num = np.sum(dudt_p * v_p * w[np.newaxis,np.newaxis,:,np.newaxis])
#    den = np.sum(dudt_p * dudt_p * w[np.newaxis,np.newaxis,:,np.newaxis])
#
#    v -= dudt * (num/den)
#    print "num= ", num
#    print "den= ", den
#
#    return (num/den)

def ddt_project(i_step,v):
    '''
    returns the magnitude of the projection of the array v onto the primal 
    time derivative dudt at time step i_step.  
    Makes v orthogonal to dudt at time step i_step
    '''
    return c_channel.c_ddt_project(int(i_step),v)




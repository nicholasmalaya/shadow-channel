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

def uxBar(i_step,profile=True):
    '''
    returns the mean x-velocity profile (or centerline x-velocity) 
    '''
    C = get_solution(i_step, copy=True)
    ux = spec2phys(C)
    ux = ux[0]
    uxbar = (ux.mean(2)).mean(0)
   
    if not(profile):
        uxbar = uxbar[c_channel.c_qpts()/2 + 1]
 
    return uxbar

def uxBarAvg(start_step,end_step,T,profile=True):
    '''
    returns the time averaged mean x-velocity profile
    '''
    if profile:
        uxbar = np.zeros(c_channel.c_qpts())
    else:
        uxbar = 0.0

    assert 0 <= start_step <= end_step <= c_channel.c_nsteps() 
    for i in range(start_step,end_step):
        uxbar = uxbar + (dt/T) * uxBar(i,profile) 

    return uxbar



def vxBar(i_step,profile=True):
    '''
    returns the mean x-velocity sensitivity
    '''
    IC = c_channel.c_getincresoln(int(i_step))
    IC = IC.copy()
    assert IC.shape == (c_channel.c_Nz(), 2, c_channel.c_dimR(),
                       c_channel.c_Nx() / 2)
    vx = spec2phys(IC)
    vx = vx[0]
    vxbar = (vx.mean(2)).mean(0)
   
    if not(profile):
        vxbar = vxbar[c_channel.c_qpts()/2 + 1]
 
    return vxbar

def vxBarAvg(start_step,end_step,T,profile=True):
    '''
    returns the time averaged mean x-velocity sensitivity profile
    '''
    if profile:
        vxbar = np.zeros(c_channel.c_qpts())
    else:
        vxbar = 0.0

   
    assert 0 <= start_step <= end_step <= c_channel.c_nsteps()
    for i in range(start_step,end_step):
        vxbar = vxbar + (dt/T) * vxBar(i,profile) 

    return vxbar

def uxBarGrad(n_chunk,chunk_bounds,T,uxbar_avg,zeta,profile=True):
    '''
    returns the time averaged sensitivity of the mean x-velocity profile
    to the reynolds number.  Inputs are the number of time chunks, the 
    starting and ending step index of each chunk (chunk_bounds), 
    total simulation time T, the time averaged mean x-velocity profile
    uxbar_avg, and time dilation terms zeta (length n_chunk vector)
    if profile = 1, sensitivity of the entire profile is compute, otherwise
    only the sesitivity of the centerline velocity is computed
    '''
    if profile:
        grad = np.zeros(c_channel.c_qpts())
    else:
        grad = 0.0    

    for i in range(n_chunk):
        start_step = chunk_bounds[i]
        end_step = chunk_bounds[i+1]
        vxbar_avg = vxBarAvg(start_step,end_step,T,profile)
        grad = grad + vxbar_avg + (1./T) * zeta[i]*(uxBar(end_step,profile)-uxbar_avg)

    return grad


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
    global __is_initialized__
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
    global __is_initialized__
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

def z_vel(i_step):
    '''
    returns the z-velocity profile  
    '''
    C = get_solution(i_step, copy=True)
    uz = spec2phys(C)
    uz = uz[0] # 2
 
    return uz

def z_vel2Avg(start_step,end_step,T,profile=True):
    '''
    returns the time averaged mean square z-velocity profile
    '''
    if profile:
        uz2avg = np.zeros(c_channel.c_qpts())
    else:
        uz2avg = 0.0
        y, w = quad()

    assert 0 <= start_step <= end_step <= c_channel.c_nsteps() 
    for i in range(start_step,end_step):
        uz = z_vel(i)
        uz2 = uz # uz2 = uz * uz
        uz2 = (uz2.mean(2)).mean(0)
        if not(profile): uz2 = (w * uz2).sum()
        uz2avg = uz2avg + (dt/T) * uz2  

    return uz2avg



def Iz_vel(i_step,project=False):
    '''
    returns the incremental z-velocity 
    '''
    IC = c_channel.c_getincresoln(int(i_step))
    IC = IC.copy()
    assert IC.shape == (c_channel.c_Nz(), 2, c_channel.c_dimR(),
                       c_channel.c_Nx() / 2)
    if project: 
        ddt_project(i_step,IC)

    vz = spec2phys(IC)
    vz = vz[0] #2
 
    return vz

def dz_vel2Avg(start_step,end_step,T,profile=True):
    '''
    returns the time averaged mean squared z-velocity sensitivity profile
    '''
    if profile:
        d_uz2avg = np.zeros(c_channel.c_qpts())
    else:
        d_uz2avg = 0.0
        y, w = quad()


   
    assert 0 <= start_step <= end_step <= c_channel.c_nsteps()
    for i in range(start_step,end_step):
        uz = z_vel(i)
        vz = Iz_vel(i,project = False)
        d_uz2 = vz # d_uz2 = 2 * uz * vz
        d_uz2 = (d_uz2.mean(2)).mean(0)
        if not(profile): d_uz2 = (w * d_uz2).sum()
        d_uz2avg = d_uz2avg + (dt/T) * d_uz2

    return d_uz2avg

def uz2BarGrad(n_chunk,chunk_bounds,T,uz2_avg,zeta,profile=True):
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
        y, w = quad()   

    for i in range(n_chunk):
        start_step = chunk_bounds[i]
        end_step = chunk_bounds[i+1]
        d_uz2Avg = dz_vel2Avg(start_step,end_step,T,profile)
        uzEnd = (z_vel(end_step).mean(2)).mean(0)
        uz2End = uzEnd # uz2End = uzEnd * uzEnd
        if not(profile): uz2End = (w * uz2End).sum()
        grad = grad + d_uz2Avg + (1./T) * zeta[i]*(uz2End-uz2_avg)

    return grad


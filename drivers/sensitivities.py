import sys
sys.path.append("..")
from scipy import *
import numpy as np
from channel import *

def vel(dirc,i_step):
    '''
    returns the velocity in direction dirc 
    '''
    C = get_solution(i_step, copy=True)
    u = spec2phys(C)
    u = u[dirc] 
 
    return u

def velAvg(dirc,start_step,end_step,T):
    '''
    returns the time averaged mean velocity profile in the dirc direction
    '''
    uavg = np.zeros(c_channel.c_qpts())

    assert 0 <= start_step <= end_step <= c_channel.c_nsteps() 
    for i in range(start_step,end_step):
        u = vel(dirc,i)
        u = (u.mean(2)).mean(0)
        uavg = uavg + (dt/T) * u  

    return uavg



def vel2Avg(dirc,start_step,end_step,T):
    '''
    returns the time averaged mean square z-velocity profile
    '''
    u2avg = np.zeros(c_channel.c_qpts())

    assert 0 <= start_step <= end_step <= c_channel.c_nsteps() 
    for i in range(start_step,end_step):
        u = vel(dirc,i)
        u2 = u * u
        u2 = (u2.mean(2)).mean(0)
        u2avg = u2avg + (dt/T) * u2  

    return u2avg



def I_vel(dirc,i_step,project=False):
    '''
    returns the incremental velocity in direction dirc
    '''
    IC = c_channel.c_getincresoln(int(i_step))
    IC = IC.copy()
    assert IC.shape == (c_channel.c_Nz(), 2, c_channel.c_dimR(),
                       c_channel.c_Nx() / 2)
    if project: 
        ddt_project(i_step,IC)

    v = spec2phys(IC)
    v = (1./(Re * Re)) * v[dirc] 
 
    return v

def dvelAvg(dirc,start_step,end_step,T):
    '''
    returns the time averaged mean squared velocity sensitivity profile
    '''
    vavg = np.zeros(c_channel.c_qpts())
   
    assert 0 <= start_step <= end_step <= c_channel.c_nsteps()
    for i in range(start_step,end_step):
        u = vel(dirc,i)
        v = I_vel(dirc,i,project = False)
        v = (v.mean(2)).mean(0)
        vavg = vavg + (dt/T) * v

    return vavg


def dvel2Avg(dirc,start_step,end_step,T):
    '''
    returns the time averaged mean squared velocity sensitivity profile
    '''
    du2avg = np.zeros(c_channel.c_qpts())


   
    assert 0 <= start_step <= end_step <= c_channel.c_nsteps()
    for i in range(start_step,end_step):
        u = vel(dirc,i)
        v = I_vel(dirc,i,project = False)
        du2 = 2 * u * v
        du2 = (du2.mean(2)).mean(0)
        du2avg = du2avg + (dt/T) * du2

    return du2avg

def uBarGrad(dirc,n_chunk,chunk_bounds,T,uavg,zeta):
    '''
    returns the time averaged sensitivity of the mean squared velocity profile
    to the reynolds number.  Inputs are the number of time chunks, the 
    starting and ending step index of each chunk (chunk_bounds), 
    total simulation time T, the time averaged mean x-velocity profile
    uxbar_avg, and time dilation terms zeta (length n_chunk vector)
    '''
    grad = np.zeros(c_channel.c_qpts())

    for i in range(n_chunk):
        start_step = chunk_bounds[i]
        end_step = chunk_bounds[i+1]
        duAvg = dvelAvg(dirc,start_step,end_step,T)
        uEnd = (vel(dirc,end_step).mean(2)).mean(0)
        grad = grad + duAvg + (1./T) * (1./(Re * Re)) * zeta[i]*(uEnd-uavg)

    return grad

def u2BarGrad(dirc,n_chunk,chunk_bounds,T,u2avg,zeta):
    '''
    returns the time averaged sensitivity of the mean squared velocity profile
    to the reynolds number.  Inputs are the number of time chunks, the 
    starting and ending step index of each chunk (chunk_bounds), 
    total simulation time T, the time averaged mean x-velocity profile
    uxbar_avg, and time dilation terms zeta (length n_chunk vector)
    '''
    grad = np.zeros(c_channel.c_qpts())

    for i in range(n_chunk):
        start_step = chunk_bounds[i]
        end_step = chunk_bounds[i+1]
        du2Avg = dvel2Avg(dirc,start_step,end_step,T)
        uEnd = (vel(dirc,end_step).mean(2)).mean(0)
        u2End = uEnd * uEnd
        grad = grad + du2Avg + (1./T) * (1./(Re * Re)) * zeta[i] * (u2End-u2avg)

    return grad



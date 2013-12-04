import sys
import os
os.environ[ 'MPLCONFIGDIR' ] = '/tmp/'
sys.path.append("..")
from pylab import *
from numpy import *
from mpi4py import MPI
from pariter import *

class Wrapper(object):
    def __init__(self, Nx, Ny, Nz, n, nb, dT, tan, adj, proj):
        self.Nx, self.Ny, self.Nz = int(Nx), int(Ny), int(Nz)
        self.m, self.n, self.dT = int(Nx*(Ny-2)*Nz), int(n), float(dT)
        self.nb = nb
        self._tan, self._adj, self._proj = tan, adj, proj

    # ------------------ interface ------------------
    def forward(self, i, v0, inhomo):
        v1 = v0.copy()
        start_step = self.nb[i]
        end_step = self.nb[i+1]
        self._proj(start_step, v1)
        self._tan(start_step, end_step, v1, inhomo)
        eta = self._proj(end_step, v1)
        return v1, eta

    def backward(self, i, w0, strength):
        w1 = w0.copy()
        start_step = self.nb[i+1]
        end_step = self.nb[i]
        self._proj(start_step, w1)
        self._adj(start_step, end_step, w1, 0, strength)
        self._proj(end_step, w1)
        return w1

    # ------------------ KKT matvec ----------------------
    def matvec(self, x, inhomo=0):
        mpi_comm = MPI.COMM_WORLD
        mpi_rank = mpi_comm.Get_rank()
        mpi_size = mpi_comm.Get_size()
        tag_w, tag_v = 800, 801
        n_v = self.m * self.n

        # print mpi_rank, mpi_size, x.shape, self.m, n_v
        x = x[::2] + 1j * x[1::2]
        v = x[:n_v].reshape([-1,self.Nz,2,self.Ny-2,self.Nx/2])
        w = x[n_v:].reshape([-1,self.Nz,2,self.Ny-2,self.Nx/2])
        # first initiate nonblocking communication of w
        w_requests = []
        if mpi_rank > 0:
            w_requests.append(mpi_comm.Isend(w[0], mpi_rank-1, tag_w))
            w = vstack([w, zeros([1,self.Nz,2,self.Ny-2,self.Nx/2],complex)])
        else:
            w = vstack([zeros([1,self.Nz,2,self.Ny-2,self.Nx/2],complex), w, zeros([1,self.Nz,2,self.Ny-2,self.Nx/2],complex)])
        assert w.shape[0] == self.n + 1
        if mpi_rank < mpi_size - 1:
            w_requests.append(mpi_comm.Irecv(w[-1], mpi_rank+1, tag_w))

        # do computation of v while w is going over the network
        R_w = zeros([self.n,self.Nz,2,self.Ny-2,self.Nx/2],complex)
        self.zeta = zeros(self.n)
        for i in range(self.n):
            vip, eta = self.forward(i, v[i], inhomo)
            self.zeta[i] = -eta #/ self.dT
            if i < self.n - 1:
                R_w[i+1] = v[i+1] - vip

        # make sure w is ready, then initiating sending the last vip
        MPI.Request.waitall(w_requests)

        v_requests = []
        if mpi_rank < mpi_size - 1:
            v_requests.append(mpi_comm.Isend(vip, mpi_rank+1, tag_v))
        if mpi_rank > 0:
            v0m = zeros([self.Nz,2,self.Ny-2,self.Nx/2],complex)
            v_requests.append(mpi_comm.Irecv(v0m, mpi_rank-1, tag_v))

        # do computation of w while the last vip is going over the network
        R_v = zeros([self.n, self.Nz,2,self.Ny-2,self.Nx/2],complex) 
        for i in range(self.n):
            wim = self.backward(i, w[i+1], strength=1.0)
            R_v[i] = w[i] - wim
        
        # match the last vip
        MPI.Request.waitall(v_requests)

        if mpi_rank > 0:
            R_w[0] = v[0] - v0m
        else:
            R_w = R_w[1:]

        return hstack([ravel(R_v), ravel(R_w)]).view(dtype=float)



import channel

mpi_comm = MPI.COMM_WORLD
mpi_rank = mpi_comm.Get_rank()
mpi_size = mpi_comm.Get_size()

# Directory for restart files
use_restart = False
file_address = '/home/blonigan/LSS/shadow-channel/'
results_dir = 'results/Lam1000/'
MR_iter = 10

# Key Parameters
# Re, {Nx, Ny, Nz} ru_steps, n_chunk, t_chunk, dt
channel.Re = 500
channel.Nx = 16
channel.Ny = 33
channel.Nz = 16
channel.dt = 0.01

n_chunk = 1
chunk_steps = 5
n_steps = n_chunk * chunk_steps + 1

# defined vector with indicies of start and end steps for each time chunk
chunk_bounds = chunk_steps * arange(0,n_chunk + 1)
# Lx, Lz, meanU
channel.Lx = 1.6
channel.Lz = 1.6
channel.meanU = 1.0

restart = None #  "keefe_runup_stage_5"

ru_steps = 0 if mpi_rank == 0 else 0

# initial condition
u0 = zeros([channel.Nz,2,channel.Ny-2,channel.Nx/2],complex) 

# solve portion of primal on each processor
if mpi_rank == 0:
    print "Starting Primal ..."

if mpi_rank > 0:
    mpi_comm.Recv(u0, mpi_rank - 1, 1)
    channel.para_init(n_steps, u0)
else:
    channel.init(n_steps,ru_steps, restart = restart)

if mpi_rank < mpi_size - 1:
    u0 = channel.get_solution(n_steps, copy=True)
    mpi_comm.Send(u0, mpi_rank + 1, 1)

if mpi_rank == (mpi_size-1):
    print "Primal Complete!"

# Compute time averaged objective function over all time (and processors), Jbar
# Compute Gradient
T = ((mpi_size - 1) * (n_steps-1) + n_steps)*channel.dt
if mpi_rank == (mpi_size-1): 
    Jbar = channel.velAvg(0,n_steps,T,profile=False)
else:
    Jbar = channel.velAvg(0,n_steps-1,T,profile=False)

Jbar = mpi_comm.allreduce(Jbar)

if mpi_rank == 0: print "Total sim time:", T, "Jbar:", Jbar


pde = Wrapper(channel.Nx, channel.Ny, channel.Nz, n_chunk, chunk_bounds,
              channel.dt * chunk_steps,
              channel.tangent, channel.adjoint, channel.ddt_project)

# construct matrix rhs
nvw = 2 * pde.n - (1 if mpi_rank == 0 else 0)
x = zeros(2 * pde.m * nvw)
rhs = -pde.matvec(x, 1) - pde.matvec(x, 0)

# solve
from scipy import sparse
import scipy.sparse.linalg as splinalg


# INITIAL GUESS HERE!
if use_restart:
    if mpi_rank == 0: print 'Reading Restart File(s)'
    filename = file_address + results_dir + 'vw_array{0}.npy'.format(mpi_rank)
    vw = load(filename)
    assert (vw.size == rhs.size)
else:
    vw = zeros(rhs.size)

oper = splinalg.LinearOperator((vw.size, vw.size), matvec=pde.matvec,
                               dtype=float)

par_dot = lambda vw1, vw2: mpi_comm.allreduce(dot(vw1, vw2))
# def par_dot(vw1, vw2):
#     return mpi_comm.allreduce(dot(vw1, vw2))
par_norm = lambda vw : sqrt(par_dot(vw, vw))

class Callback:
    'Convergence monitor'
    def __init__(self, pde, T, Jbar):
        self.n = 0
        self.pde = pde
        self.T = T
        self.Jbar = Jbar
        self.hist = []

    def __call__(self, x):
        self.n += 1
        if self.n == 1 or self.n % 10 == 0:
            mpi_rank = mpi_comm.Get_rank()
            mpi_size = mpi_comm.Get_size()
            resnorm = par_norm(self.pde.matvec(x, 1))
            nb = self.pde.nb.copy()
            if mpi_rank < (mpi_size - 1):
                nb[-1] = nb[-1] - 1
            # Compute Gradient
            grad = channel.uavgGrad(self.pde.n,nb,self.T,self.Jbar,self.pde.zeta,profile=False)
            grad = mpi_comm.allreduce(grad)
            if mpi_rank == 0: print 'iter ', self.n, resnorm, grad
            self.hist.append([self.n, resnorm, grad])
        sys.stdout.flush()

# --- solve with minres (if cg converges this should converge -#
callback = Callback(pde, T, Jbar)
callback(vw)
vw, info = par_minres(oper, rhs, vw, par_dot, maxiter=MR_iter, tol=1E-6,
                      callback=callback)

pde.matvec(vw, 1)


# save vw array to file with mpi_rank in name...
filename = file_address + results_dir + 'vw_array{0}.npy'.format(mpi_rank)
save(filename,vw)

# gradient plots
import sensitivities as sens
dirc = 0

if mpi_rank == (mpi_size-1):
    uxbar_avg = sens.velAvg(dirc,0,n_steps,T)
else:
    uxbar_avg = sens.velAvg(dirc,0,n_steps-1,T)

len = uxbar_avg.shape[0]
for i in range(len):
    uxbar_avg[i] = mpi_comm.allreduce(uxbar_avg[i])

nb = chunk_bounds


grad = sens.uBarGrad(dirc,n_chunk,nb,T,uxbar_avg,pde.zeta)
for i in range(len):
    grad[i] = mpi_comm.allreduce(grad[i])

if mpi_rank == 0:
    y,w = channel.quad()

    figure()
    plot(uxbar_avg,y)
    figure()
    plot(grad,y)
    # save data
    filename = file_address + results_dir + 'grad_profile.npy'
    save(filename,grad)
    
# v and w plots

vxhist = []
vxprojhist = []
strt = nb[0]
if (mpi_rank == mpi_size - 1):
    fnsh = nb[-1]
else:
    fnsh = nb[-1]-1

for i in range(strt,fnsh+1):
    vx = (sens.I_vel(dirc,i,project=False).mean(2)).mean(0)
    vxproj = (sens.I_vel(dirc,i,project=True).mean(2)).mean(0)
    vxhist.append(vx)   
    vxprojhist.append(vxproj)   

vxhist = array(vxhist)
vxprojhist = array(vxprojhist)
zetas = pde.zeta

# send data to process 0
tag_vx = 501
tag_vxproj = 502
tag_zeta = 503
vx_requests = []

if mpi_rank > 0:
    vx_requests.append(mpi_comm.Send(vxhist, 0,tag_vx))    
    vx_requests.append(mpi_comm.Send(vxprojhist, 0,tag_vxproj))  
    vx_requests.append(mpi_comm.Send(zetas, 0,tag_zeta))  


if mpi_rank == 0:
    # receive data from other processes, append
    for i in range(1,mpi_size):
        if i == (mpi_size - 1):
            vxtmp = zeros([n_steps, len])
            vxprojtmp = zeros([n_steps, len])
        else:
            vxtmp = zeros([n_steps-1, len])
            vxprojtmp = zeros([n_steps-1, len])
 
        zetatmp = zeros(pde.zeta.size)
        vx_requests.append(mpi_comm.Recv(vxtmp, i, tag_vx))
        vx_requests.append(mpi_comm.Recv(vxprojtmp, i, tag_vxproj))
        vx_requests.append(mpi_comm.Recv(zetatmp, i, tag_zeta))

        vxhist = vstack([vxhist,vxtmp]) 
        vxprojhist = vstack([vxprojhist,vxprojtmp]) 
        zetas = vstack([zetas,zetatmp])
    
    # plot!
    t = channel.dt * arange(int(T/channel.dt))
    y,w = channel.quad()

    figure()
    contourf(y, t, vxprojhist, 100); colorbar()
    axis([-1.0, 1.0, t[0], t[-1]])
    show()
    

    # save data
    filename = file_address + results_dir + 'tan_hist.npy'
    save(filename,vxhist)
    filename = file_address + results_dir + 'tan_proj_hist.npy'
    save(filename,vxprojhist)
    filename = file_address + results_dir + 'zetas.npy'            
    save(filename,zetas)


        

# clean memory
channel.destroy()

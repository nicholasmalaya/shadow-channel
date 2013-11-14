import sys
sys.path.append("..")
from pylab import *
from numpy import *
#from scipy.interpolate import interp1d
#from scipy.integrate import ode
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
        #print "fwd", i, norm(v1)
        return v1, eta

    def backward(self, i, w0, strength):
        w1 = w0.copy()
        start_step = self.nb[i+1]
        end_step = self.nb[i]
        self._proj(start_step, w1)
        self._adj(start_step, end_step, w1, 0, strength)
        self._proj(end_step, w1)
        #print "bwk", i, norm(w1)
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
            self.zeta[i] = eta / self.dT
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

# Key Parameters
# Re, {Nx, Ny, Nz} ru_steps, n_chunk, t_chunk, dt
channel.Re = 2000
channel.Nx = 16
channel.Ny = 33
channel.Nz = 16
channel.dt = 0.01

n_chunk = 2
chunk_steps = 5
n_steps = n_chunk * chunk_steps + 1

# defined vector with indicies of start and end steps for each time chunk
chunk_bounds = chunk_steps * arange(0,n_chunk + 1)
# Lx, Lz, meanU
channel.Lx = 1.6
channel.Lz = 1.6
channel.meanU = 1.0

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
    restart = "keefe_runup_stage_5"
    channel.init(n_steps,ru_steps, restart = restart)

if mpi_rank < mpi_size - 1:
    u0 = channel.get_solution(n_steps, copy=True)
    mpi_comm.Send(u0, mpi_rank + 1, 1)

if mpi_rank == (mpi_size-1):
    print "Primal Complete!"

# Compute time averaged objective function over all time (and processors), Jbar
# Jbar = mpi_comm.allreduce(kuramoto.cvar.JBAR)
# kuramoto.c_assignJBAR(Jbar)

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

vw = zeros(rhs.size)
oper = splinalg.LinearOperator((vw.size, vw.size), matvec=pde.matvec,
                               dtype=float)

par_dot = lambda vw1, vw2: mpi_comm.allreduce(dot(vw1, vw2))
# def par_dot(vw1, vw2):
#     return mpi_comm.allreduce(dot(vw1, vw2))
par_norm = lambda vw : sqrt(par_dot(vw, vw))

class Callback:
    'Convergence monitor'
    def __init__(self, pde):
        self.n = 0
        self.pde = pde
        self.hist = []

    def __call__(self, x):
        self.n += 1
        if self.n == 1 or self.n % 10 == 0:
            resnorm = par_norm(self.pde.matvec(x, 1))
            # gradient = mpi_comm.allreduce(kuramoto.c_grad())
            if mpi_rank == 0: print 'iter ', self.n, resnorm
            self.hist.append([self.n, resnorm])
        sys.stdout.flush()

# --- solve with minres (if cg converges this should converge -#
callback = Callback(pde)
callback(rhs * 0)
vw, info = par_minres(oper, rhs, vw, par_dot, maxiter=100, tol=1E-6,
                      callback=callback)

pde.matvec(vw, 1)

# output data to file with mpi_rank in name...

# clean memory
channel.destroy()

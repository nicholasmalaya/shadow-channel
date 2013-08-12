import sys
from pylab import *
from numpy import *
from scipy.interpolate import interp1d
from scipy.integrate import ode
from mpi4py import MPI
from pariter import *

class Wrapper(object):
    def __init__(self, m, n, dT, tan, adj, proj):
        self.m, self.n, self.dT = int(m), int(n), float(dT)
        self._tan, self._adj, self._proj = tan, adj, proj

    # ------------------ interface ------------------
    def forward(self, i, v0, inhomo):
        v1 = v0.copy()
        self._tan(i, v1, inhomo)
        eta = self._proj(i + 1, 0, v1)
        return v1, eta
        
    def backward(self, i, w0, strength):
        w1 = w0.copy()
        self._adj(i, w1, strength)
        self._proj(i, 0, w1)
        return w1

    # ------------------ KKT matvec ----------------------
    def matvec(self, x, inhomo=0):
        mpi_comm = MPI.COMM_WORLD
        mpi_rank = mpi_comm.Get_rank()
        mpi_size = mpi_comm.Get_size()
        tag_w, tag_v = 800, 801

        n_v = self.m * self.n

        # print mpi_rank, mpi_size, x.shape, self.m, n_v
        v = x[:n_v].reshape([-1, self.m])
        w = x[n_v:].reshape([-1, self.m])

        # first initiate nonblocking communication of w
        w_requests = []
        if mpi_rank > 0:
            w_requests.append(mpi_comm.Isend(w[0], mpi_rank-1, tag_w))
            w = vstack([w, zeros([1,self.m])])
        else:
            w = vstack([zeros([1,self.m]), w, zeros([1,self.m])])
        assert w.shape[0] == self.n + 1
        if mpi_rank < mpi_size - 1:
            w_requests.append(mpi_comm.Irecv(w[-1], mpi_rank+1, tag_w))

        # do computation of v while w is going over the network
        R_w = zeros([self.n, self.m])
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
            v0m = zeros(self.m)
            v_requests.append(mpi_comm.Irecv(v0m, mpi_rank-1, tag_v))

        # do computation of w while the last vip is going over the network
        R_v = zeros([self.n, self.m])
        for i in range(self.n):
            wim = self.backward(i, w[i+1], strength=.01)
            R_v[i] = w[i] - wim

        # match the last vip
        MPI.Request.waitall(v_requests)

        if mpi_rank > 0:
            R_w[0] = v[0] - v0m
        else:
            R_w = R_w[1:]

        return hstack([ravel(R_v), ravel(R_w)])



import kuramoto
# c, n_grid, T0, n_chunk, t_chunk, dt_max

mpi_comm = MPI.COMM_WORLD
mpi_rank = mpi_comm.Get_rank()
mpi_size = mpi_comm.Get_size()

T0 = 500 if mpi_rank == 0 else 0

u0 = rand(127)
if mpi_rank > 0:
    mpi_comm.Recv(u0, mpi_rank - 1, 1)

kuramoto.c_init(0.5, u0, 500, 50, 2, 0.2);

if mpi_rank < mpi_size - 1:
    mpi_comm.Send(u0, mpi_rank + 1, 1)

pde = Wrapper(kuramoto.cvar.N_GRID, kuramoto.cvar.N_CHUNK,
              kuramoto.cvar.DT_STEP * kuramoto.cvar.N_STEP,
              kuramoto.c_tangent, kuramoto.c_adjoint, kuramoto.c_project_ddt)

# construct matrix rhs
nvw = 2 * pde.n - (1 if mpi_rank == 0 else 0)
x = zeros(pde.m * nvw)
rhs = pde.matvec(x, -1) - pde.matvec(x, 0)

# solve
from scipy import sparse
import scipy.sparse.linalg as splinalg

vw = zeros(rhs.size)
oper = splinalg.LinearOperator((vw.size, vw.size), matvec=pde.matvec,
                               dtype=float)

par_dot = lambda vw1, vw2: mpi_comm.allreduce(dot(vw1, vw2))
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
            if mpi_rank == 0: print 'iter ', self.n, resnorm
            self.hist.append([self.n, resnorm])
        sys.stdout.flush()

# --- solve with minres (if cg converges this should converge -#
callback = Callback(pde)
callback(rhs * 0)
vw, info = par_minres(oper, rhs, vw, par_dot, maxiter=100, tol=1E-6,
                      callback=callback)

pde.matvec(vw, 1)
'''
u, v, v0, uEnd = [], [], [], []
for i in range(kuramoto.cvar.N_CHUNK):
    for j in range(kuramoto.cvar.N_STEP):
        u.append(kuramoto.c_u(i, j))
        v.append(kuramoto.c_v(i, j))
        v0.append(v[-1].copy())
        kuramoto.c_project_ddt(i, j, v0[-1])
    uEnd.append(kuramoto.c_u(i, kuramoto.cvar.N_STEP))

u, v, v0, uEnd = array(u), array(v), array(v0), array(uEnd)

subplot(1,2,1); contourf(u, 100); colorbar()
subplot(1,2,2); contourf(v0, 100); colorbar()

du = v.mean(0) + (pde.zeta[:,newaxis] * (uEnd - u.mean(0))).mean(0)
figure()
plot(u.mean(0) - 0.1 * du)
plot(u.mean(0) + 0.1 * du)
'''

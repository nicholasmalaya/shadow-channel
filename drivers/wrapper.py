import sys
sys.path.append("..")
from pylab import *
from numpy import *
#from scipy.interpolate import interp1d
#from scipy.integrate import ode
#from lssode import *

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
        
    def backward(self, i, w0, strength, inhomo):
        w1 = w0.copy()
        start_step = self.nb[i+1]
        end_step = self.nb[i]
        self._proj(start_step, w1)
        self._adj(start_step, end_step, w1, inhomo, strength)
        self._proj(end_step, w1)
        return w1

    # ------------------ KKT matvec ----------------------
    def matvec(self, x, inhomo=0):
        print "matvec called, inhomo =", inhomo
        n_v = self.m * self.n
        n_w = self.m * (self.n - 1)
        assert x.shape == (n_v + n_w,)
        v = x[:n_v].reshape([-1,self.Nz,2,self.Ny-2,self.Nx/2])
        w = vstack([x[n_v:].reshape([-1,self.Nz,2,self.Ny-2,self.Nx/2]), zeros([1,self.Nz,2,self.Ny-2,self.Nx/2],complex)])

        R_w = zeros([self.n - 1, self.Nz,2,self.Ny-2,self.Nx/2],complex)
        R_v = zeros([self.n, self.Nz,2,self.Ny-2,self.Nx/2],complex)
        self.zeta = zeros(self.n)
        for i in range(self.n):
            vip, eta = self.forward(i, v[i], inhomo)
            self.zeta[i] = eta / self.dT
            wim = self.backward(i, w[i], strength=1.0, inhomo=0)
            if i < self.n - 1:
                R_w[i] = v[i+1] - vip
            if i > 0:
                R_v[i] = w[i-1] - wim
            else:
                R_v[i] = - wim
        return hstack([ravel(R_v), ravel(R_w)])


import channel 
# Key Parameters
# Re, {Nx, Ny, Nz} ru_steps, n_chunk, t_chunk, dt
channel.Re = 2000
channel.Nx = 16
channel.Ny = 33
channel.Nz = 16
channel.dt = 0.01

ru_steps = 0
n_chunk = 3
chunk_steps = 10
n_steps = n_chunk * chunk_steps + 1 

# defined vector with indicies of start and end steps for each time chunk
chunk_bounds = chunk_steps * arange(0,n_chunk + 1)
# Lx, Lz, meanU
channel.Lx = 1.6
channel.Lz = 1.6
channel.meanU = 1.0

# Initial Condition
restart = "keefe_runup_stage_5"

channel.init(n_steps,ru_steps, restart = restart)

pde = Wrapper(channel.Nx, channel.Ny, channel.Nz, n_chunk, chunk_bounds, 
              channel.dt * chunk_steps,
              channel.tangent, channel.adjoint, channel.ddt_project)

# construct matrix rhs
x = zeros(pde.m * (2 * pde.n - 1),complex)
rhs = -pde.matvec(x, 1) - pde.matvec(x, 0)

# solve
from scipy import sparse
import scipy.sparse.linalg as splinalg

w = zeros(rhs.size)
oper = splinalg.LinearOperator((w.size, w.size), matvec=pde.matvec, dtype=float)

class Callback:
    'Convergence monitor'
    def __init__(self, pde):
        self.n = 0
        self.pde = pde
        self.hist = []

    def __call__(self, x):
        self.n += 1
        if self.n == 1 or self.n % 10 == 0:
            resnorm = norm(self.pde.matvec(x, 1))
            print('iter ', self.n, resnorm)
            self.hist.append([self.n, resnorm])
        sys.stdout.flush()

# --- solve with minres (if cg converges this should converge -#
callback = Callback(pde)
callback(rhs * 0)
vw, info = splinalg.minres(oper, rhs, maxiter=100, tol=1E-6,
                           callback=callback)

#pde.matvec(vw, 1)
#channel.destroy()

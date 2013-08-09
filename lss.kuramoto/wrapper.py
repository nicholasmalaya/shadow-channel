import sys
sys.path.append('/home/qiqi/git/lssode')
from pylab import *
from numpy import *
from scipy.interpolate import interp1d
from scipy.integrate import ode
from lssode import *

class Wrapper(object):
    def __init__(self, m, n, tan, adj, proj):
        self.m, self.n = int(m), int(n)
        self._tan, self._adj, self._proj = tan, adj, proj

    # ------------------ interface ------------------
    def forward(self, i, v0, inhomo):
        v1 = v0.copy()
        self._tan(i, v1, inhomo)
        eta = self._proj(i + 1, v1)
        return v1, eta
        
    def backward(self, i, w0, strength):
        w1 = w0.copy()
        self._adj(i, w1, strength)
        self._proj(i, w1)
        return w1

    # ------------------ KKT matvec ----------------------

    def matvec(self, x, inhomo=0):
        n_v = self.m * self.n
        n_w = self.m * (self.n - 1)
        assert x.shape == (n_v + n_w,)
        v = x[:n_v].reshape([-1, self.m])
        w = vstack([x[n_v:].reshape([-1, self.m]), zeros([1,self.m])])

        R_w = zeros([self.n - 1, self.m])
        R_v = zeros([self.n, self.m])
        for i in range(self.n):
            vip, eta = self.forward(i, v[i], inhomo)
            wim = self.backward(i, w[i], strength=.01)
            if i < self.n - 1:
                R_w[i] = v[i+1] - vip
            if i > 0:
                R_v[i] = w[i-1] - wim
            else:
                R_v[i] = - wim
        return hstack([ravel(R_v), ravel(R_w)])


import kuramoto
# c, n_grid, T0, n_chunk, t_chunk, dt_max
kuramoto.c_init(0.5, 127, 500, 100, 2, 0.2);

pde = Wrapper(kuramoto.cvar.N_GRID, kuramoto.cvar.N_CHUNK,
              kuramoto.c_tangent, kuramoto.c_adjoint, kuramoto.c_project_ddt)

# construct matrix rhs
x = zeros(pde.m * (2 * pde.n - 1))
rhs = pde.matvec(x, -1) - pde.matvec(x, 0)

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
            print 'iter ', self.n, resnorm
            self.hist.append([self.n, resnorm])
        sys.stdout.flush()

# --- solve with minres (if cg converges this should converge -#
callback = Callback(pde)
callback(rhs * 0)
vw, info = splinalg.minres(oper, rhs, maxiter=100, tol=1E-6,
                           callback=callback)

pde.matvec(vw, -1)

u, v = [], []
for i in range(kuramoto.cvar.N_CHUNK):
    for j in range(kuramoto.cvar.N_STEP):
        u.append(kuramoto.c_u(i, j))
        v.append(kuramoto.c_v(i, j))

u = array(u)
v = array(v)

subplot(1,2,1); contourf(u, 100); colorbar()
subplot(1,2,2); contourf(v, 100); colorbar()
show()

import sys
sys.path.append('/home/qiqi/git/lssode')
from pylab import *
from numpy import *
from scipy.interpolate import interp1d
from scipy.integrate import ode
from lssode import *

class Wrapper(object):
    def __init__(self, m, n, dT, tan, adj, proj):
        self.m, self.n, self.dT = int(m), int(n), float(dT)
        self._tan, self._adj, self._proj = tan, adj, proj

    # ------------------ interface ------------------
    def forward(self, i, v0, inhomo):
        v1 = v0.copy()
        self._proj(i, 0, v1)
        self._tan(i, v1, inhomo)
        eta = self._proj(i + 1, 0, v1)
        return v1, eta
        
    def backward(self, i, w0, strength, inhomo):
        w1 = w0.copy()
        self._proj(i + 1, 0, w1)
        self._adj(i, w1, strength, inhomo)
        self._proj(i, 0, w1)
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
        self.zeta = zeros(self.n)
        for i in range(self.n):
            vip, eta = self.forward(i, v[i], inhomo=0)
            self.zeta[i] = eta / self.dT
            wim = self.backward(i, w[i], 0.01, inhomo)
            if i < self.n - 1:
                R_w[i] = v[i+1] - vip
            if i > 0:
                R_v[i] = w[i-1] - wim
            else:
                R_v[i] = - wim
        return hstack([ravel(R_v), ravel(R_w)])


import kuramoto
# c, n_grid, T0, n_chunk, t_chunk, dt_max
# u0 = rand(127)
u0 = zeros(127)
u0[64] = 1

kuramoto.c_init(0.5, u0, 500, 25, 4, 0.2);

pde = Wrapper(kuramoto.cvar.N_GRID, kuramoto.cvar.N_CHUNK,
              kuramoto.cvar.DT_STEP * kuramoto.cvar.N_STEP,
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

pde.matvec(vw, 1)

u, v, w, v0, uEnd = [], [], [], [], []
for i in range(kuramoto.cvar.N_CHUNK):
    for j in range(kuramoto.cvar.N_STEP):
        u.append(kuramoto.c_u(i, j))
        v.append(kuramoto.c_v(i, j))
        w.append(kuramoto.c_w(i, j))
        v0.append(v[-1].copy())
        kuramoto.c_project_ddt(i, j, v0[-1])
    uEnd.append(kuramoto.c_u(i, kuramoto.cvar.N_STEP))

u, v, w, v0, uEnd = array(u), array(v), array(w),  array(v0), array(uEnd)

# Plots

m = u.shape[0]
u = hstack([zeros([m,1]), u, zeros([m,1])])
v0 = hstack([zeros([m,1]), v0, zeros([m,1])])
w = hstack([zeros([m,1]), w, zeros([m,1])])

x = linspace(0,128,129)
dt = kuramoto.cvar.DT_STEP
T1 = kuramoto.cvar.N_CHUNK * kuramoto.cvar.N_STEP * dt
t = dt * arange(int(T1 / dt))

#subplot(2,2,1)
title('$u(t)$')
#contourf(x, t, u, 100); colorbar()
contour(x,t,u,[-3, 4]);
xticks( array([0., 32., 64., 96., 128.]  ))
xlabel('$x$')
ylabel('$t$')
axis('scaled')
axis([x[0], x[-1], t[0], t[-1]])

figure()
#subplot(2,2,3)
title('$\hat{v}(t)$')
contourf(x, t, v0, 100); colorbar()
contour(x,t,u,[0.3*u.min() , 0.3*u.max()]);
xticks( array([0., 32., 64., 96., 128.]  ))
xlabel('$x$')
ylabel('$t$')
axis('scaled')
axis([x[0], x[-1], t[0], t[-1]])

figure()
#subplot(2,2,4)
title('$\hat{w}(t)$')
contourf(x, t, w, 100); colorbar()
contour(x,t,u,[0.3*u.min() , 0.3*u.max()]);
xticks( array([0., 32., 64., 96., 128.]  ))
xlabel('$x$')
ylabel('$t$')
axis('scaled')
axis([x[0], x[-1], t[0], t[-1]])


#subplot(1,2,1); contourf(u, 100); colorbar()
#subplot(1,2,2); contourf(v0, 100); colorbar()
#figure()
#subplot(1,2,1); contourf(u, 100); colorbar()
#subplot(1,2,2); contourf(w, 100); colorbar()

#du = v.mean(0) + (pde.zeta[:,newaxis] * (uEnd - u.mean(0))).mean(0)
#figure()
#plot(u.mean(0) - 0.1 * du)
#plot(u.mean(0) + 0.1 * du)

show()

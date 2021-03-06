import sys
sys.path.append("..")
from pylab import *
from numpy import *

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
        n_v = self.m * self.n
        n_w = self.m * (self.n - 1)
        assert x.shape == (2*n_v + 2*n_w,)
        x = x[::2] + 1j * x[1::2]

        v = x[:n_v].reshape([-1,self.Nz,2,self.Ny-2,self.Nx/2])
        w = vstack([x[n_v:].reshape([-1,self.Nz,2,self.Ny-2,self.Nx/2]), zeros([1,self.Nz,2,self.Ny-2,self.Nx/2],complex)])
        R_w = zeros([self.n - 1, self.Nz,2,self.Ny-2,self.Nx/2],complex)
        R_v = zeros([self.n, self.Nz,2,self.Ny-2,self.Nx/2],complex)
        self.zeta = zeros(self.n)
        for i in range(self.n):
            vip, eta = self.forward(i, v[i], inhomo)
            self.zeta[i] = -eta #/ self.dT
            wim = self.backward(i, w[i], strength=1.0)
            if i < self.n - 1:
                R_w[i] = v[i+1] - vip
            if i > 0:
                R_v[i] = w[i-1] - wim
            else:
                R_v[i] = - wim
        return hstack([ravel(R_v), ravel(R_w)]).view(dtype=float)


import channel 
# Key Parameters
# Re, {Nx, Ny, Nz} ru_steps, n_chunk, t_chunk, dt
channel.Re = 500
channel.Nx = 16
channel.Ny = 33
channel.Nz = 16
channel.dt = 0.01

ru_steps = 0
n_chunk = 2
chunk_steps = 5
n_steps = n_chunk * chunk_steps + 1 

# defined vector with indicies of start and end steps for each time chunk
chunk_bounds = chunk_steps * arange(0,n_chunk + 1)
# Lx, Lz, meanU
channel.Lx = 1.6
channel.Lz = 1.6
channel.meanU = 1.0

# Initial Condition
#restart = "keefe_runup_stage_5"
restart = None

channel.init(n_steps,ru_steps, restart = restart)


# compute time averaged objective function Jbar
T = n_steps * channel.dt 
Jbar = channel.velAvg(0,n_steps,T,profile=False)

print "Total sim time:", T, "Jbar:", Jbar

pde = Wrapper(channel.Nx, channel.Ny, channel.Nz, n_chunk, chunk_bounds, 
              channel.dt * chunk_steps,
              channel.tangent, channel.adjoint, channel.ddt_project)

# construct matrix rhs
x = zeros(2 * pde.m * (2 * pde.n - 1))
rhs = - pde.matvec(x,1) - pde.matvec(x, 0)
print "RHS norm: ", norm(rhs)
# solve
from scipy import sparse
import scipy.sparse.linalg as splinalg

w = zeros(rhs.size)
oper = splinalg.LinearOperator((w.size, w.size), matvec=pde.matvec, dtype=float)

class Callback:
    'Convergence monitor'
    def __init__(self, pde, rhs, T, Jbar):
        self.n = 0
        self.pde = pde
        self.rhs = rhs
        self.T = T
        self.Jbar = Jbar
        self.hist = []

    def __call__(self, x):
        self.n += 1
        if self.n == 1 or self.n % 10 == 0:
            resnorm1 = norm(self.pde.matvec(x,0)-self.rhs)
            resnorm = norm(self.pde.matvec(x, 1))
            # Compute Gradient
            grad = channel.uavgGrad(self.pde.n,self.pde.nb,self.T,self.Jbar,self.pde.zeta,profile=False)
            print('iter ', self.n, resnorm, resnorm1, grad)
            self.hist.append([self.n, resnorm, grad])
        sys.stdout.flush()

# --- solve with minres (if cg converges this should converge -#

# read in restart

# vw0 = load("vw_array.npy")
vw0 = 0 * rhs

callback = Callback(pde,rhs,T,Jbar)
callback(vw0.copy())


vw, info = splinalg.minres(oper, rhs, x0=vw0, maxiter=10, tol=1E-6,
                           callback=callback)


pde.matvec(vw, 1)
#channel.destroy()
# v and w plots
vxhist = []
for i in range(0,n_steps):
    vx = (channel.Ivel(i,project=True).mean(2)).mean(0)
    vxhist.append(vx)

t = channel.dt * arange(n_steps)
y,w = channel.quad()

figure()
contourf(y, t, vxhist, 100); colorbar()
axis([-1.0, 1.0, t[0], t[-1]])

# gradient plots
uxavg = channel.velAvg(0,n_steps,T,profile=True)
grad = channel.uavgGrad(n_chunk,chunk_bounds,T,uxavg,pde.zeta,profile=True)

figure()
mag = 1
#plot(uz2bar_avg,y,uz2bar_avg + mag*grad,y,uz2bar_avg - mag*grad,y)
plot(grad,y)
show()
'''
import sensitivities as sens
dirc = 0
y,w = channel.quad()
# gradient plots
ux2bar_avg = sens.vel2Avg(dirc,0,n_steps,T)
gradu2 = sens.u2BarGrad(dirc,n_chunk,chunk_bounds,T,ux2bar_avg,pde.zeta)

uxbar_avg = sens.velAvg(dirc,0,n_steps,T)
gradu = sens.uBarGrad(dirc,n_chunk,chunk_bounds,T,uxbar_avg,pde.zeta)
 
        
figure()
plot(gradu2,y,gradu,y)

# v plot
vxhist = []
for i in range(0,n_steps):
    vx = (sens.I_vel(0,i,project=True).mean(2)).mean(0)
    vxhist.append(vx)
 
vxhist = array(vxhist)
t = channel.dt * arange(n_steps)

figure()
title('U')
contourf(y, t, vxhist, 100); colorbar()
axis([-1.0, 1.0, t[0], t[-1]])
'''


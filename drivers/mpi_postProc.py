import sys
import os # delete
os.environ[ 'MPLCONFIGDIR' ] = '/tmp/' #delete
sys.path.append("..")
from pylab import *  # change to scipy
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
        MPI.Request.waitall(w_requests) # capitalize

        v_requests = []
        if mpi_rank < mpi_size - 1:
            v_requests.append(mpi_comm.Isend(vip, mpi_rank+1, tag_v))
        if mpi_rank > 0:
            v0m = zeros([self.Nz,2,self.Ny-2,self.Nx/2],complex)
            v_requests.append(mpi_comm.Irecv(v0m, mpi_rank-1, tag_v))

        # do computation of w while the last vip is going over the network
        R_v = zeros([self.n, self.Nz,2,self.Ny-2,self.Nx/2],complex) 
        for i in range(self.n):
            wim = self.backward(i, w[i+1], strength=0.01)
            R_v[i] = w[i] - wim
        
        # match the last vip
        MPI.Request.waitall(v_requests) # capitalize

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
file_address = '/home/blonigan/LSS/'
results_dir = 'shadow-channel/'

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

pde = Wrapper(channel.Nx, channel.Ny, channel.Nz, n_chunk, chunk_bounds,
              channel.dt * chunk_steps,
              channel.tangent, channel.adjoint, channel.ddt_project)

# construct matrix rhs
nvw = 2 * pde.n - (1 if mpi_rank == 0 else 0)
x = zeros(2 * pde.m * nvw)


# INITIAL GUESS HERE!
if mpi_rank == 0: print 'Reading Restart File(s)'
filename = file_address + results_dir + 'vw_array{0}.npy'.format(mpi_rank)
vw = load(filename)
assert (vw.size == x.size)

pde.matvec(vw, 1)


import sensitivities as sens
T = ((mpi_size - 1) * (n_steps-1) + n_steps)*channel.dt

# gradient plots
if mpi_rank == (mpi_size-1):
    uxbar_avg = sens.velAvg(0,0,n_steps,T)
    uybar_avg = sens.velAvg(1,0,n_steps,T)
    uzbar_avg = sens.velAvg(2,0,n_steps,T)
    ux2bar_avg = sens.vel2Avg(0,0,n_steps,T)
    uy2bar_avg = sens.vel2Avg(1,0,n_steps,T)
    uz2bar_avg = sens.vel2Avg(2,0,n_steps,T)
else:
    uxbar_avg = sens.velAvg(0,0,n_steps-1,T)
    uybar_avg = sens.velAvg(1,0,n_steps-1,T)
    uzbar_avg = sens.velAvg(2,0,n_steps-1,T)
    ux2bar_avg = sens.vel2Avg(0,0,n_steps-1,T)
    uy2bar_avg = sens.vel2Avg(1,0,n_steps-1,T)
    uz2bar_avg = sens.vel2Avg(2,0,n_steps-1,T)

len = uz2bar_avg.shape[0]
for i in range(len):
    uxbar_avg[i] = mpi_comm.allreduce(uxbar_avg[i])
    uybar_avg[i] = mpi_comm.allreduce(uybar_avg[i])
    uzbar_avg[i] = mpi_comm.allreduce(uzbar_avg[i])
    ux2bar_avg[i] = mpi_comm.allreduce(ux2bar_avg[i])
    uy2bar_avg[i] = mpi_comm.allreduce(uy2bar_avg[i])
    uz2bar_avg[i] = mpi_comm.allreduce(uz2bar_avg[i])

nb = chunk_bounds

if mpi_rank < (mpi_size - 1):
    nb[-1] = nb[-1] - 1

gradux = sens.uBarGrad(0,n_chunk,nb,T,uxbar_avg,pde.zeta)
graduy = sens.uBarGrad(1,n_chunk,nb,T,uybar_avg,pde.zeta)
graduz = sens.uBarGrad(2,n_chunk,nb,T,uzbar_avg,pde.zeta)
gradux2 = sens.u2BarGrad(0,n_chunk,nb,T,ux2bar_avg,pde.zeta)
graduy2 = sens.u2BarGrad(1,n_chunk,nb,T,uy2bar_avg,pde.zeta)
graduz2 = sens.u2BarGrad(2,n_chunk,nb,T,uz2bar_avg,pde.zeta)
for i in range(len):
    gradux[i] = mpi_comm.allreduce(gradux[i])
    graduy[i] = mpi_comm.allreduce(graduy[i])
    graduz[i] = mpi_comm.allreduce(graduz[i])
    gradux2[i] = mpi_comm.allreduce(gradux2[i])
    graduy2[i] = mpi_comm.allreduce(graduy2[i])
    graduz2[i] = mpi_comm.allreduce(graduz2[i])

if mpi_rank == 0:
    y,w = channel.quad()

    figure()
    plot(uxbar_avg,y)
    figure()
    plot(gradux,y)
    # save data
    filename = file_address + results_dir + 'gradux_profile.npy'
    save(filename,gradux)
    filename = file_address + results_dir + 'graduy_profile.npy'
    save(filename,graduy)
    filename = file_address + results_dir + 'graduz_profile.npy'
    save(filename,graduz)
    filename = file_address + results_dir + 'gradux2_profile.npy'
    save(filename,gradux2)
    filename = file_address + results_dir + 'graduy2_profile.npy'
    save(filename,graduy2)
    filename = file_address + results_dir + 'graduz2_profile.npy'
    save(filename,graduz2)
    
# v and w plots

vxhist = []
vyhist = []
vzhist = []
for i in range(nb[0],nb[-1]+1):
    vx = (sens.I_vel(0,i,project=True).mean(2)).mean(0)
    vy = (sens.I_vel(1,i,project=True).mean(2)).mean(0)
    vz = (sens.I_vel(2,i,project=True).mean(2)).mean(0)
    vxhist.append(vx)
    vyhist.append(vy)
    vzhist.append(vz)   

vxhist = array(vxhist)
vyhist = array(vyhist)
vzhist = array(vzhist)

# send data to process 0
tag_vx = 501
tag_vy = 502
tag_vz = 503
vx_requests = []
vy_requests = []
vz_requests = []

if mpi_rank > 0:
    vx_requests.append(mpi_comm.Send(vxhist, 0,tag_vx))    
    vy_requests.append(mpi_comm.Send(vyhist, 0,tag_vy))
    vz_requests.append(mpi_comm.Send(vzhist, 0,tag_vz))

if mpi_rank == 0:
    # receive data from other processes, append
    for i in range(1,mpi_size):
        if i == (mpi_size - 1):
            vxtmp = zeros([n_steps, len])
            vytmp = zeros([n_steps, len])
            vztmp = zeros([n_steps, len])
        else:
            vxtmp = zeros([n_steps-1, len])
            vytmp = zeros([n_steps-1, len])
            vztmp = zeros([n_steps-1, len])
        vx_requests.append(mpi_comm.Recv(vxtmp, i, tag_vx))
        vy_requests.append(mpi_comm.Recv(vytmp, i, tag_vy))
        vz_requests.append(mpi_comm.Recv(vztmp, i, tag_vz))

        vxhist = vstack([vxhist,vxtmp]) 
        vyhist = vstack([vyhist,vytmp])       
        vzhist = vstack([vzhist,vztmp]) 
 
    # plot!
    t = channel.dt * arange(int(T/channel.dt))
    y,w = channel.quad()

    figure()
    contourf(y, t, vxhist, 100); colorbar()
    axis([y[0], y[-1], t[0], t[-1]])
    show()
    
    # save data
    filename = file_address + results_dir + 'tan_ux_hist.npy'
    save(filename,vxhist)
    filename = file_address + results_dir + 'tan_uy_hist.npy'
    save(filename,vyhist)
    filename = file_address + results_dir + 'tan_uz_hist.npy'
    save(filename,vzhist)
 

# clean memory
channel.destroy()

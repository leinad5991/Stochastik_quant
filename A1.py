import numpy as np
import matplotlib.pyplot as pl

tmax=0.1
xmax=10.
Nt=10000
Nx=100
dt=tmax/Nt
dx=xmax/Nx

Q=np.zeros((Nt,Nx))
Q2=np.zeros((Nt,Nx))

def step(W,Q,dx,dt,nt,N):

    for i in range(1,N-1):
        dw=np.sqrt(dt)*np.random.randn()
        Q[nt][i]=Q[nt-1][i]+(Q[nt-1][i+1]-2*Q[nt-1][i]+2*Q[nt-1][i-1])/dx**2 *dt - Q[nt-1][i]*dt + W[i]*dw
    return Q
	
def braun(N,dx):
    W=np.array([0])
    for i in range(N):
        dw=np.sqrt(dx)*np.random.randn()
        W=np.append(W,W[-1]+dw)
    return W

for i in range(Nt-1):
    W=braun(Nx,dx)
    Q=step(W,Q,dx,dt,i+1,Nx)



pl.plot(Q[-1])

pl.show()

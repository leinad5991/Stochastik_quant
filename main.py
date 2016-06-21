import numpy as np
import matplotlib.pyplot as pl

dt=0.0001
dx=0.1
Nt=100000
Nx=100
tmax = dt*Nt
xmax = dx*Nx

Q=np.zeros((Nt,Nx))

def step(Q,dx,dt,nt,N):

    for i in range(1,N-1):
        dw=np.sqrt(2*dt/dx)*np.random.randn()
        Q[nt][i]=Q[nt-1][i]+ (Q[nt-1][i+1]-2*Q[nt-1][i]+Q[nt-1][i-1])/dx**2 * dt - Q[nt-1][i]*dt + dw
    return Q

def Energy(q,dx,xmax):
    v=(q[0:-2]-q[1:-1])/dx
    E=1/2*v**2+1./2.*q[1:-1]**2
    return sum(E)/xmax

E=0
for i in range(1,Nt):
    Q=step(Q,dx,dt,i,Nx)
    E=E+Energy(Q[i],dx,xmax)

print E/Nt
pl.plot(sum(Q)/Nt)
pl.show()
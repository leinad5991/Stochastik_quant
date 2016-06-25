import numpy as np
import matplotlib.pyplot as pl

dt=0.00002
dx=0.1
Nt=50000
Nx=300
tmax = dt*Nt
xmax = dx*Nx

w=1.
m=1.
k=w**2*m

def step(Q,dx,dt,N):
    Qn=np.zeros(N)

    #Periodische Randbedingung
    dw=np.sqrt(2*dt/dx)*np.random.randn()
    Qn[0]=Q[0]+ m*(Q[1]-2*Q[0]+Q[-1])/dx**2 * dt -k**Q[0]*dt + dw
    dw=np.sqrt(2*dt/dx)*np.random.randn()
    Qn[-1]=Q[-1]+ m*(Q[0]-2*Q[-1]+Q[-2])/dx**2 * dt - k*Q[-1]*dt + dw

    for i in range(1,N-1):
        dw=np.sqrt(2*dt/dx)*np.random.randn()
        Qn[i]=Q[i]+ m*(Q[i+1]-2*Q[i]+Q[i-1])/dx**2 * dt - k*Q[i]*dt + dw
    return Qn

def histogram(Q):
    global density_axes
    hist, density_axes = np.histogram(Q, bins=Nx , range=(-5,5),density=True)
    return hist

def Energy(q,dx,xmax):
    v=(q[0:-2]-q[1:-1])/dx
    E=1/2*m*v**2+1./2.*k*q[1:-1]**2
    return sum(E)/xmax

def Energy2(q):
    global density_axes
    min1=density_axes[0]
    max1=density_axes[-1]
    int=(max1-min1)
    dx=int/Nx
    w=np.sqrt(q)

    Ekin=-1/2.*w[1:-2]*(w[0:-3]+w[2:-1]-2*w[1:-2])/dx**2
    Epot=1/2.*density_axes[0:-1]**2*q

    return (sum(Epot)+sum(Ekin))*dx

E=0
Qm=np.zeros(Nx)
q=0

Q=np.zeros(Nx)
for i in range(1,Nt):
    Q=step(Q,dx,dt,Nx)

    Qm=histogram(Q)+Qm
    if i%1000==0:
        print Energy2(Qm/i)

print Energy2(Qm/Nt)
pl.plot(density_axes[0:-1],Qm/Nt)
pl.show()
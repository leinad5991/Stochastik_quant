import numpy as np
import matplotlib.pyplot as pl

dt=0.001
dx=0.1
Nt=100000
Nx=50
tmax = dt*Nt
xmax = 2*np.pi

dx=xmax/Nx


w=1.
m=1.
k=w**2*m

def step(Q,Qh,dx,dt,N,h):
    Qn=np.zeros(N)
    Qhn=np.zeros(N)
    #Periodische Randbedingung
    dw=np.sqrt(2*dt/dx)*np.random.randn()
    Qn[0]=Q[0]+ m*(Q[1]-2*Q[0]+Q[-1])/dx**2 * dt -V(Q[0])*dt + dw
    Qhn[0]=Qh[0]+ m*(Qh[1]-2*Qh[0]+Qh[-1])/dx**2 * dt -V(Qh[0])*dt + dw+ h*dt


    dw=np.sqrt(2*dt/dx)*np.random.randn()
    Qn[-1]=Q[-1]+ m*(Q[0]-2*Q[-1]+Q[-2])/dx**2 * dt - V(Q[-1])*dt + dw
    Qhn[-1]=Qh[-1]+ m*(Qh[0]-2*Qh[-1]+Qh[-2])/dx**2 * dt - V(Qh[-1])*dt + dw +h*dt

    for i in range(1,N-1):
        dw=np.sqrt(2*dt/dx)*np.random.randn()
        Qn[i]=Q[i]+ m*(Q[i+1]-2*Q[i]+Q[i-1])/dx**2 * dt - V(Q[i])*dt + dw
        Qhn[i]=Qh[i]+ m*(Qh[i+1]-2*Qh[i]+Qh[i-1])/dx**2 * dt - V(Qh[i])*dt + dw
    return (Qn,Qhn)

def V(x):
    return x

def histogram(Q):
    global density_axes
    hist, density_axes = np.histogram(Q, bins=Nx , range=(-5,5),density=True)
    return hist

def Energy(q,dx,xmax):
    v=(q[0:-2]-q[1:-1])/dx
    E=1/2*m*v**2+1./2.*k*q[1:-1]**2
    return sum(E)/xmax

def Propagator(Q,Nx):
    P=np.zeros(Nx)
    for i in range(Nx):
            P[i]=Q[i]*Q[0]
    return P

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


Q=np.zeros(Nx)
Qh=np.zeros(Nx)

Qs=np.zeros(Nx)
Qhs=np.zeros(Nx)

h=0.01
for i in range(1000):
    (Q,Qh)=step(Q,Qh,dx,dt,Nx,h)

for i in range(1,Nt):
    (Q,Qh)=step(Q,Qh,dx,dt,Nx,h)
    Qs=Qs+Q
    Qhs=Qhs+Qh
    if i%10000==0:
        pl.plot(((Qhs-Qs)/i/h+Qs[0]*Qs/i**2).real)
        pl.show()

#print Energy2(Qm/Nt)
pl.plot((Qhs-Qs)/Nt/h+Qs[0]*Qs/Nt**2)
pl.show()

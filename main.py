import numpy as np
import matplotlib.pyplot as pl
import matplotlib.animation as anim

dt=0.1


Nx=10
dx1=1.
xmax=Nx*dx1

xaxes=np.arange(0,xmax,dx1)
dx=np.zeros(Nx)+dx1



w=1
m=1.
k=w**2*m




def step(Q,Qh,dx,dt,N,h):

    dw=np.sqrt(2*dt/dx)*np.random.randn(N)

    Qn=Q+ m*(np.roll(Q,1)-2*Q+np.roll(Q,-1))/dx**2 * dt -V(Q)*dt + dw
    Qhn=Qh+ m*(np.roll(Qh,1)-2*Qh + np.roll(Qh,-1))/dx**2 * dt -V(Qh)*dt + dw

    Qhn[0]=Qhn[0]+ h*dt

    return (Qn,Qhn)


def Vcosh(x):
    b=0.1**(1/2.)
    V0=1.
    return 2*b*V0/np.cosh(b*x)**2* np.tanh(b* x)

def potcosh(x):
    b=0.1**(1/2.)
    V0=1.
    return -V0/np.cosh(b*x)**2


def Vdopp(x):
    x0=3.
    return 1/8.*x0**2.* (4.*x* (x**2./x0-1.))/x0

def potdopp(x):
    x0=3.
    return 1/8.*x0**2.*(x**2./x0**2. - 1.)**2.


def histogram(Q):
    global density_axes
    hist, density_axes = np.histogram(Q, bins=Nx , range=(-10,10),density=True)
    return hist

def Propagator(Q,Nx):
    P=np.zeros(Nx)
    for i in range(Nx):
            P[i]=Q[i]*Q[0]
    return P
def steigung(x,dx):
    return (-x+np.roll(x,-1))/dx

def Energy(q):
    global density_axes,PP
    min1=density_axes[0]
    max1=density_axes[-1]
    int=(max1-min1)
    dx=int/Nx
    w=np.sqrt(q)

    Ekin=-1/2.*w[1:-2]*(w[:-3]+w[2:-1]-2*w[1:-2])/dx**2
    n1=0
    n2=Nx-1
    while w[n1]==0 and n1<Nx-1: n1=n1+1
    while w[n2]==0 and n2>0:
        n2=n2-1
    #if n1!=0: print Ekin[n1:n2]

    Epot=pot(density_axes[0:-1])*q

    return (sum(Epot)+sum(Ekin))*dx

Q=np.zeros(Nx)
Qh=Q


Qs=np.zeros(Nx)
Qhs=Qs
PP=np.zeros(Nx-1)
P=0
P0=0

V=Vcosh
pot=potcosh


h=0.0000000001
for i in range(10000):
    (Q,Qh)=step(Q,Qh,dx,dt,Nx,h)


def multistep(n):
    global Q,Qh,P,Qs,Qhs,P0,dt
    P0=0
    Qfiktiv=np.zeros(n)
    print Qfiktiv
    for i in range(n):
        (Q,Qh)=step(Q,Qh,dx,dt,Nx,h)
        Qfiktiv[i]=sum((Qhs-Qs)/h/i)*dx1

        #P0=P0+histogram(Q)

        Qs=Qs+Q
        Qhs=Qhs+Qh

    P=P+P0
    pl.plot(np.log(Qfiktiv))
    print steigung(np.log(Qfiktiv),dt)
    pl.show()

def extract(Q):
    global PP,xaxes
    st=steigung(np.log(Q),dx)
    v1=sum(st[-10:])/10

    norm1=Q/np.exp(v1*xaxes)


def animate(i):
    global Q,Qh,P,P0,Qs,Qhs,xaxes,density_axes,h,P0,vplot,PP
    mn=1000
    multistep(mn)

    st=steigung(np.log(  (Qhs-Qs)/i/mn/h   ),dx)
    Energy(P)
    #extract((Qhs-Qs)/i/mn/h)
    #line.set_data(xaxes,st)
    line.set_data(xaxes,st)

    print max(st)
    print -min(st)
    print ""
    return line,

def init():
    line.set_data([], [])
    return line,
fig = pl.figure()
#ax = pl.axes(xlim=(0, 10000), ylim=(0, 0.5))
line, = pl.plot([], [])
histogram(Q)


multistep(6000)
#pl.plot(density_axes[1:],pot(density_axes[1:])/20.)


#anima = anim.FuncAnimation(fig, animate, init_func=init, blit=True)
#pl.show()


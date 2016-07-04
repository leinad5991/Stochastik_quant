import numpy as np
import matplotlib.pyplot as pl
import matplotlib.animation as anim

dt=0.01
dx=1.
Nt=1000000
Nx=100
tmax = dt*Nt
#xmax = 8*np.pi

dx=1
xmax=dx*Nx

w=1
m=1.
k=w**2*m



def step(Q,Qh,dx,dt,N,h):
    Qn=np.zeros(N)
    Qhn=np.zeros(N)
    #Periodische Randbedingung
    dw=np.sqrt(2*dt/dx)*np.random.randn()
    Qn[0]=Q[0]+ m*(Q[1]-2*Q[0]+Q[-1])/dx**2 * dt -V(Q[0])*dt + dw
    Qhn[0]=Qh[0]+ m*(Qh[1]-2*Qh[0]+Qh[-1])/dx**2 * dt -V(Qh[0])*dt + dw + h*dt


    dw=np.sqrt(2*dt/dx)*np.random.randn()
    Qn[-1]=Q[-1]+ m*(Q[0]-2*Q[-1]+Q[-2])/dx**2 * dt - V(Q[-1])*dt + dw
    Qhn[-1]=Qh[-1]+ m*(Qh[0]-2*Qh[-1]+Qh[-2])/dx**2 * dt - V(Qh[-1])*dt + dw

    for i in range(1,N-1):
        dw=np.sqrt(2*dt/dx)*np.random.randn()
        mhm=0
        if i==N/2:mhm=h*dt
        Qn[i]=Q[i]+ m*(Q[i+1]-2*Q[i]+Q[i-1])/dx**2 * dt - V(Q[i])*dt + dw
        Qhn[i]=Qh[i]+ m*(Qh[i+1]-2*Qh[i]+Qh[i-1])/dx**2 * dt - V(Qh[i])*dt + dw

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
    x0=5.
    return 1/8.*x0**2.* (4.*x* (x**2./x0-1.))/x0

def potdopp(x):
    x0=5.
    return 1/8.*x0**2.*(x**2./x0**2. - 1.)**2.


def histogram(Q):
    global density_axes
    hist, density_axes = np.histogram(Q, bins=Nx , range=(-10,10),density=True)
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
def steigung(x,dx):

    return (x[1:]-x[:-1])/dx

def Energy2(q):
    global density_axes
    min1=density_axes[0]
    max1=density_axes[-1]
    int=(max1-min1)
    dx=int/Nx
    w=np.sqrt(q)

    Ekin=-1/2.*w[1:-2]*(w[0:-3]+w[2:-1]-2*w[1:-2])/dx**2
    Epot=pot(density_axes[0:-1])*q

    return (sum(Epot)+sum(Ekin))*dx


Q=np.zeros(Nx)
Qh=Q

Qs=np.zeros(Nx)
Qhs=Qs

P=0
P0=0

xaxes=np.arange(0,xmax,dx)

V=Vdopp
pot=potdopp


h=0.01
for i in range(1000):
    (Q,Qh)=step(Q,Qh,dx,dt,Nx,h)


def multistep(n):
    global Q,Qh,P,Qs,Qhs,P0
    P0=0
    for i in range(n):
        (Q,Qh)=step(Q,Qh,dx,dt,Nx,h)

        P0=P0+histogram(Q)

        Qs=Qs+Q
        Qhs=Qhs+Qh
    P=P+P0
def animate(i):
    global Q,Qh,P,P0,Qs,Qhs,xaxes,density_axes,h,P0,vplot
    mn=1000
    multistep(mn)

    line.set_data(density_axes[1:],P0/mn)
    print (Energy2(P/i/mn) ,steigung(np.log(  (Qhs-Qs)/i/mn/h ),dx)[-1])
    vplot.set_data(density_axes[1:],P0/mn)
    #line.set_data(xaxes[1:],steigung(np.log((Qhs-Qs)/i/mn/h),dx))
    return line,

def init():
    line.set_data([], [])
    return line,

fig = pl.figure()
ax = pl.axes(xlim=(-10, 10), ylim=(0, 0.6))
line, = pl.plot([], [])

histogram(Q)
pl.plot(density_axes[1:],pot(density_axes[1:])/8.)
vplot,=pl.plot([],[])

anima = anim.FuncAnimation(fig, animate, init_func=init, blit=True)
pl.show()


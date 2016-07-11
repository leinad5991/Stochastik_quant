import numpy as np
import matplotlib.pyplot as pl
import matplotlib.animation as anim

Nt=1000000
Nx=100
#xmax = 8*np.pi


xmax=10.

dx=xmax/Nx
w=1
m=1.
k=w**2*m

def QQ(x,E):
    x0=5.
    return 2.*(E-1/8.*x0**2.*(x**2./x0**2. - 1.)**2.)

def sig(E):
    for i in range(2, Nx):
        W[i] = 2 * W[i - 1] -  QQ((i - 1) * dx, E) * W[i - 1] * (dx) ** 2 - W[i - 2]
    return np.sign(W[Nx-1])

E0 = 0.489205
E1 = 1.420638
E = E0
#E = E1

W = np.zeros(Nx)
W[0]=0
W[1]=0.001


for i in range(2,Nx):
    W[i]=2*W[i-1]-QQ((i-1)*dx,E)*W[i-1]*(dx)**2-W[i-2]


'''
E=0.02 #Vorrausgesetzt, dass E0>0.02

while E<2.:
    F = E-0.01
    if sig(E) == sig(F):
        pass
    else:
        while F<E:
            G = F+0.000001
            if sig(F)==sig(G):
                pass
            else:
                print F
            F+=0.000001
    E += 0.01
'''
pl.plot(W)
pl.show()
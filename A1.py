import numpy as np
import matplotlib.pyplot as pl

tmax=1.
xmax=1.
N=100
dt=tmax/N
dx=xmax/N
X=np.array([1])
Q=np.zeros((N,N))
W=np.array([0])

t=np.arange(N+1)*dt

def step(W,Q,dx,nt,N):
    dw=np.sqrt(dt)*np.random.randn()
    for i in range(1,N-1):

        Q[nt][i]=Q[nt][i-1]+(Q[nt-1][i+1]-2*Q[nt-1][i]+2*Q[nt-1][i-1])/dx**2 *dt - Q[nt-1][i]*dt + W[i]*dw
        print "QQQQ"
        print (Q[nt-1][i+1]-2*Q[nt-1][i]+2*Q[nt-1][i-1])/dx**2 *dt
        print Q[nt-1][i]*dt
        print W[i]*dw
    return Q
	
def braun(N,dx):
    W=np.array([0])
    for i in range(N):
        dw=np.sqrt(dx)*np.random.randn()
        W=np.append(W,W[-1]+dw)
    return W

W=braun(N,dx)
for i in range(N-1):
    W=braun(N,dx)
    Q=step(W,Q,dx,i+1,N)

pl.plot(Q[100-2])

pl.show()

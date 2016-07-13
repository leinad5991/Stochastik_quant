import numpy as np
import matplotlib.pyplot as pl
import matplotlib.animation as anim

Nx = 100

xmax = 10.

dx = xmax / Nx

w = 1
m = 1.
k = w ** 2 * m


def QQ(x, E):
    x0 = 5.
    return 2. * (E - 1 / 8. * x0 ** 2. * (x ** 2. / x0 ** 2. - 1.) ** 2.)


def sig(E, N=10000000, dx=0.000001):
    W = np.zeros(N)
    W[0] = 0
    W[1] = 0.001
    i = 1
    while abs(W[i]) < 100 / dx and i < N - 1:
        i = i + 1
        W[i] = 2 * W[i - 1] - QQ((i - 1) * dx, E) * W[i - 1] * (dx) ** 2 - W[i - 2]
    return np.sign(W[i])


def findEnergy(E2, N=10000000, dx=0.000001):
    E2 = np.longdouble(E2)
    acuracy = 0.0000000000001
    s1 = sig(E2, N=1000, dx=0.01)
    dE = np.longdouble(0.1)
    while sig(E2, N, dx) == s1:
        E1 = E2
        E2 = E2 + dE

    while abs(E2 - E1) > acuracy:
        mean = (E2 + E1) / 2.
        if sig(mean, N, dx) == s1:
            E1 = mean
        else:
            E2 = mean
    return mean


print findEnergy(0.5, N=1000000, dx=0.001)

# Es funktioniert jetzt man kann die Energien jetzt viel exakter ermitteln


E0 = 0.4892043
E0 = 0.489205
# E1 = 1.420638
# E = E0
E = E0

W = np.zeros(Nx)
W[0] = 0
W[1] = 0.001
WM = np.zeros(Nx)
WM[1] = 0
WM[0] = 0.001
# Das Problem war das man die afangsbedingugen auch spiegeln musst :)
# wenn du noch comiten kannst währe nett
for i in range(2, Nx):
    W[i] = 2 * W[i - 1] - QQ((i - 1) * dx, E) * W[i - 1] * (dx) ** 2 - W[i - 2]

for i in range(2, Nx):
    WM[i] = 2 * WM[i - 1] - QQ((i - 3) * dx, E) * WM[i - 1] * (dx) ** 2 - WM[i - 2]

'''

E=0.4 #Vorrausgesetzt, dass E0>0.02

while E<2.:
    F = E-0.01
    if sig(E) == sig(F):
        pass
    else:
        while F<E:
            G = F+0.0000001
            if sig(F)==sig(G):
                pass
            else:
                print F
            F+=0.0000001
    E += 0.01
'''
xaxes = np.arange(-Nx + 2, Nx - 2, 1.) * dx
pl.plot(xaxes, -QQ(xaxes, E0))
pl.plot(xaxes, np.concatenate((WM[2:][::-1], W[:-2]), axis=0))
pl.show()

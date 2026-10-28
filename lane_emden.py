import numpy as np
import matplotlib.pyplot as plt

def solveLaneEmden(xi0, th0, y0, h, n):
    coords = {"xi": [], "th": []}
    while th0 > 0:
        dy = lambda xi, th, y: -th**n - 2/xi*y if xi != 0 else -1/3
        k1th = h*y0
        k1y = h*dy(xi0, th0, y0)
        k2th = h*(y0+k1y/2)
        k2y = h*dy(xi0+h/2, th0+k1th/2, y0+k1y/2)
        k3th = h*(y0+k2y/2)
        k3y = h*dy(xi0+h/2, th0+k2th/2, y0+k2y/2)
        k4th = h*(y0+k3y)
        k4y = h*dy(xi0+h, th0+k3th, y0+k3y)
        th1 = th0 + 1/6 * (k1th + 2*k2th + 2*k3th + k4th)
        y1 = y0 + 1/6 * (k1y + 2*k2y + 2*k3y + k4y)
        coords["xi"].append(xi0+h)
        coords["th"].append(th1)
        xi0 = xi0+h
        th0 = th1
        y0 = y1
    return coords


c0 = solveLaneEmden(0,1,0,0.001,0)
c1 = solveLaneEmden(0,1,0,0.001,1)
c2 = solveLaneEmden(0,1,0,0.001,2)
c3 = solveLaneEmden(0,1,0,0.001,3)
c4 = solveLaneEmden(0,1,0,0.001,4)

plt.plot(c0["xi"], c0["th"])
plt.plot(c1["xi"], c1["th"])
plt.plot(c2["xi"], c2["th"])
plt.plot(c3["xi"], c3["th"])
plt.plot(c4["xi"], c4["th"])
plt.show()

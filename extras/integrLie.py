import numpy as np
from numpy.linalg import norm
from scipy.linalg import expm
from scipy.optimize import fsolve

def skw(x):
    return np.array([[0,-x[2],x[1]],[x[2],0,-x[0]],[-x[1],x[0],0]])


def skwinv(x):
    return np.array([x[2,1],x[0,2],x[1,0]])


def dexpinvso3(a,b):
    Theta = norm(a)
    theta = Theta * 0.5
    aa = skw(a)
    bb = skw(b)
    operator = np.eye(3)-0.5*aa-(Theta*(np.cos(theta)/np.sin(theta))-2)/(2*Theta**2) * aa@aa
    return operator@bb


def implieeul(A, f, y0, h):
    res = lambda x, x0, dt : - x + skwinv(dexpinvso3(x, dt*A(expm(skw(x))@x0)))
    F1 = fsolve(res, f, args=(y0, h))
    return expm(skw(F1))@y0


def implliemidp(A, f, y0, h):
    res = lambda x, x0, dt : - x + skwinv(dexpinvso3(0.5*x, dt*A(expm(skw(0.5*x))@x0)))
    F1 = fsolve(res, f, args=(y0, h))
    return expm(skw(F1))@y0
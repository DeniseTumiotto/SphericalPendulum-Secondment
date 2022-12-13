import numpy as np
from numpy.linalg import norm
from scipy.linalg import expm
from scipy.optimize import fsolve

def skw(x):
    return np.array([[0,-x[2],x[1]],[x[2],0,-x[0]],[-x[1],x[0],0]])

def skwinv(x):
    return np.array([x[2,1],x[0,2],x[1,0]])

def expso3(x):
    theta = norm(x)
    if theta > 1e-16:
        fun1 = np.sin(theta)/theta
        fun2 = (1-np.cos(theta))/theta**2
    else:
        fun1 = 1.0 - theta**2/6
        fun2 = 1.0/2 - theta**2/24 + theta**4/720
    return np.eye(3) + fun1*skw(x) + fun2*skw(x)@skw(x)

def dexpinvso3(a,b):
    Theta = norm(a)
    theta = Theta * 0.5
    aa = skw(a)
    if Theta > 1e-16:
        fun = (1-theta*(np.cos(theta)/np.sin(theta)))/(Theta**2)
    else:
        fun =  1/12 + Theta**2/720 + Theta**4/(30240)
    operator = np.eye(3) - 0.5*aa + fun*aa@aa
    return operator@b

def tantrso3(x):
    theta = norm(x)
    if theta > 1e-16:
        fun1 = (np.cos(theta)-1)/theta**2
        fun2 = (1-(np.sin(theta)/theta))/theta**2
    else:
        fun1 = - 1.0/2 + theta**2/24 - theta**4/720
        fun2 = 1/6 - theta**2/120 + theta**4/5040
    rslt = np.eye(3) + fun1*skw(x) + fun2*skw(x)@skw(x)
    return rslt.transpose()

def dexpinvse3(a,b):
    # scrivi qualcosa qui!!
    return

def expse3(x):
    u = np.copy(x[:3])
    v = np.copy(x[3:])
    A11 = expso3(u) # 3x3
    A12 = tantrso3(u)@v # 3x1
    return np.array([A11, A12])

def implieeul(A, f, y0, h):
    if f.shape == (3,) or f.shape == (3,1):
        dexp = dexpinvso3
        exp = lambda x: expm(skw(x))
    else:
        dexp = dexpinvse3
        exp = expse3
    res = lambda x, x0, dt : - x + dexp(x, dt*A(exp(x)@x0))
    F1 = fsolve(res, f, args=(y0, h))
    return expm(skw(F1))@y0

def implliemidp(A, f, y0, h):
    if f.shape == (3,) or f.shape == (3,1):
        dexp = dexpinvso3
        exp = lambda x: expm(skw(x))
    else:
        dexp = dexpinvse3
        exp = expse3
    res = lambda x, x0, dt : - x + dexp(0.5*x, dt*A(exp(0.5*x)@x0))
    F1 = fsolve(res, f, args=(y0, h))
    return expm(skw(F1))@y0

def impllietrap(A, f, y0, h):
    if f.shape == (3,) or f.shape == (3,1):
        dexp = dexpinvso3
        exp = lambda x: expm(skw(x))
    else:
        dexp = dexpinvse3
        exp = expse3
    res = lambda x, x0, dt : - x + dexp(0.5*x+0.5*x0, dt*A(exp(0.5*x+0.5*x0)@x0))
    F1 = fsolve(res, f, args=(y0, h))
    return expm(0.5*skw(F1)+0.5*skw(y0))@y0
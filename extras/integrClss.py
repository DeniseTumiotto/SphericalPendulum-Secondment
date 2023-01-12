from numpy.linalg import norm
from scipy.optimize import fsolve

def midptrule(f, y, y0, h):
    res = lambda x, x0, dt: - x + x0 + dt*f((x0+x)/2)
    y = fsolve(res, y, args=(y0, h))
    return y

def sphmidpt(f,y,y0,h):
    res = lambda x, x0, dt: - x + x0 + dt*f((x0+x)/norm(x0+x))
    y = fsolve(res, y, args=(y0, h))
    return y
import numpy as np
from numpy.linalg import norm
import integrLie as lie

# Helpful functions
def normalize(x):
    return x/norm(x)

# Define problem


# Integration parameters
t0 = 0
te = 1
nTime = 100
timeVec = np.linspace(t0,te,nTime)
dt = timeVec[1]-timeVec[0]

# Problem parameters
m0 = 1
n = normalize(-2*np.ones((3,1))+4*np.random.rand(3,1))
inertia = m0 * np.eye(3)

# initialize
alpha0 = np.ones((3,1))
w0 = np.random.rand(3,1)

YY = np.zeros((6,nTime+1))

YY[:3,0] = np.copy(alpha0)
YY[3:,0] = np.copy(w0)


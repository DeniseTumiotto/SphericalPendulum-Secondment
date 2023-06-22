import numpy as np
from scipy.linalg import expm, logm, sqrtm
from scipy.optimize import fsolve
from alive_progress import alive_bar
import matplotlib.pyplot as plt

def is_pos_def(x):
    return np.all(np.linalg.eigvals(x) > 0)

def metric(x, u, v):
    inv_x = np.linalg.inv(x)
    return np.trace(inv_x @ u @ inv_x @ v)

def my_grad(x, y):
    n = int(np.sqrt(len(x)))
    x = x.reshape((n, n))
    x1_2 = sqrtm(x)
    xminus1_2 = np.linalg.inv(sqrtm(x))
    return - x1_2 @ logm(xminus1_2 @ y @ xminus1_2) @ x1_2

def gile(y, y0, h, data_Y):
    n = y.size
    res = lambda x, x0, dt, yy : (- x0 + exponential(x, dt, grad_V(x, yy))).reshape(n)
    rslt = fsolve(res, y, args=(y0, h, data_Y))
    return rslt

def exponential(x, t, v):
    n = int(np.sqrt(len(x)))
    x = x.reshape((n, n))
    x1_2 = sqrtm(x)
    xminus1_2 = np.linalg.inv(sqrtm(x))
    return x1_2 @ expm(t * xminus1_2 @ v @ xminus1_2) @ x1_2

def logarithm(x, y):
    n = int(np.sqrt(len(x)))
    x = x.reshape((n, n))
    x1_2 = sqrtm(x)
    xminus1_2 = np.linalg.inv(sqrtm(x))
    return x1_2 @ logm(xminus1_2 @ y @ xminus1_2) @ x1_2

def distance(x, y):
    A = np.linalg.inv(x) @ y
    my_eig = np.log(np.linalg.eigvals(A))
    return np.linalg.norm(my_eig)

def grad_V(x,y):
    k = len(y)
    n = int(np.sqrt(len(x)))
    rslt = np.zeros((n, n))
    for i in range(k):
        rslt = rslt + my_grad(x, y[i])
    return rslt

k = 5
n = 3
epsilon = 1e-5
Y = []

for i in range(k):
    A = 2 * np.random.rand(n, n) - 1
    Y.append(A.transpose() @ A + epsilon * np.eye(n))

t0 = 0
te = 0.5
dt = 1e-3
my_time = np.arange(t0, te, dt)
sol1 = []
sol2 = []
dist = []

A0 = 2 * np.random.rand(n, n) - 1
Y0 = A0.transpose() @ A0 + epsilon * np.eye(n)
sol1.append(Y0)
A0 = 2 * np.random.rand(n, n) - 1
Y0 = A0.transpose() @ A0 + epsilon * np.eye(n)
sol2.append(Y0)
dist.append(distance(sol1[0], sol2[0]))

# evaluate flow and distance

with alive_bar(my_time.size, bar = 'smooth') as bar:
    for counter, the_time in enumerate(my_time):
        # crr_sol1 = gile(sol1[counter], sol1[0], the_time, Y)
        # crr_sol2 = gile(sol2[counter], sol2[0], the_time, Y)
        crr_sol1 = gile(sol1[counter], sol1[0], dt, Y)
        crr_sol2 = gile(sol2[counter], sol2[0], dt, Y)
        if not is_pos_def(crr_sol1.reshape((n, n))):
            print('Not a pos def matrix!')
        sol1.append(crr_sol1.reshape((n, n)))
        sol2.append(crr_sol2.reshape((n, n)))
        dist.append(distance(sol1[counter+1], sol2[counter+1]))
        bar()

# plot distance

plt.rcParams["figure.figsize"] = [10.00, 6.00]

fig1 = plt.figure()
plt.plot(np.append(my_time, te), dist, label='flow distance')
plt.plot(np.append(my_time, te), dist[0]*np.ones((my_time.size + 1, )), label='initial distance')
plt.legend()
plt.grid()

plt.show()

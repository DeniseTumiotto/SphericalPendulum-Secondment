import numpy as np
from scipy.linalg import expm, logm, sqrtm
from scipy.optimize import fsolve
from alive_progress import alive_bar
import matplotlib.pyplot as plt
import datetime
import pickle
import time
import winsound
import seaborn as sns

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

def gie(y, y0, h, data_Y, pp):
    n = y.size
    res = lambda x, x0, dt, yy, p : (- x0 + exponential(x, dt, grad_V(x, yy, p))).reshape(n)
    rslt = fsolve(res, y, args=(y0, h, data_Y, pp))
    return rslt

def expInt(y, y0, h, data_Y, pp):
    n = y.size
    rslt = exponential(y0.reshape((n,)), h, grad_V(y.reshape((n,)), data_Y, pp)).reshape(n)
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
    xminus1_2 = np.linalg.inv(sqrtm(x))
    # A = np.linalg.inv(x) @ y
    A = xminus1_2 @ y @ xminus1_2
    my_eig = np.log(np.linalg.eigvals(A))
    return np.linalg.norm(my_eig)

def grad_V(x,y, param):
    k = len(y)
    n = int(np.sqrt(len(x)))
    rslt = np.zeros((n, n))
    for i in range(k):
        rslt = rslt + param * 0.5 * my_grad(x, y[i])
    return rslt

k = 5
n = 2
epsilon = 1e-5
Y = []

for i in range(k):
    A = 2 * np.random.rand(n, n) - 1
    Y.append(A.transpose() @ A + epsilon * np.eye(n))

t0 = 0
dt = 1e-1
te = 2 + dt
my_time = np.arange(t0, te, dt)
time_tosave = my_time.tolist()
sol1 = []
sol2 = []
sol_expl1 = []
sol_expl2 = []
dist = []
dist_expl = []

A0 = 2 * np.random.rand(n, n) - 1
Y0 = A0.transpose() @ A0 + epsilon * np.eye(n)
sol1.append(Y0)
sol_expl1.append(Y0)
A0 = 2 * np.random.rand(n, n) - 1
Y0 = A0.transpose() @ A0 + epsilon * np.eye(n)
sol2.append(Y0)
sol_expl2.append(Y0)
# dist.append(distance(sol1[0], sol2[0]))

my_param = 1.0

# evaluate flow and distance

with alive_bar(my_time.size, bar = 'smooth') as bar:
    for counter, the_time in enumerate(my_time):
        crr_sol1 = gie(sol1[counter], sol1[0], the_time, Y, my_param)
        crr_sol2 = gie(sol2[counter], sol2[0], the_time, Y, my_param)
        crr_sol_expl1 = expInt(sol1[0], sol1[0], the_time, Y, my_param)
        crr_sol_expl2 = expInt(sol2[0], sol2[0], the_time, Y, my_param)
        if not is_pos_def(crr_sol1.reshape((n, n))):
            print('Not a pos def matrix!')
        sol1.append(crr_sol1.reshape((n, n)))
        sol2.append(crr_sol2.reshape((n, n)))
        sol_expl1.append(crr_sol_expl1.reshape((n, n)))
        sol_expl2.append(crr_sol_expl2.reshape((n, n)))
        dist.append(distance(sol1[counter+1], sol2[counter+1]))
        dist_expl.append(distance(sol_expl1[counter+1], sol_expl2[counter+1]))
        bar()

# now = datetime.datetime.now().strftime('%Y%m%dT%H%M%S')

# fileName = open('C:/Users/denis/Desktop/PHD/code/SphericalPendulum-Secondment/out/'+now+'dist.pkl', 'wb')
# pickle.dump(dist, fileName, -1)
# fileName.close()
# fileName = open('C:/Users/denis/Desktop/PHD/code/SphericalPendulum-Secondment/out/'+now+'time.pkl', 'wb')
# pickle.dump(time_tosave, fileName, -1)
# fileName.close()

# plot distance

plt.rcParams["figure.figsize"] = [10.00, 6.00]

fig1 = plt.figure(figsize=(8,4))
sns.set_style("darkgrid")
# plt.plot(my_time, dist, label='flow distance')
# plt.plot(my_time, dist[0]*np.ones((my_time.size, )), label='initial distance')
plt.plot(my_time, dist_expl, linewidth=2, label='explicit integrator')
plt.plot(my_time, dist, linewidth=2, label='GIE')
plt.legend()
plt.grid(True)
plt.xlabel('h')
plt.ylabel('d')
plt.savefig('SPD_example.pdf')

winsound.PlaySound("SystemExit", winsound.SND_ALIAS)

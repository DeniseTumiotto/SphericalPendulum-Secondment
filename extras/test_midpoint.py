import integrLie as lie
import integrClss as cl
import numpy as np
import matplotlib.pyplot as plt

def rDistance(y,z):
    return 2*np.arcsin(np.linalg.norm(y-z,2)/2)


e3 = np.array([0,0,1])
gGrav = 9.81
d = 0.0
k = 100.0

tol = 1.0e-10

t0 = 0.0
te = 10.0
h = np.arange(0.1,te,0.1)
N = np.size(h)

phi0 = -np.pi/2+(np.pi/4)*np.random.rand()
q00 = np.array([np.cos(phi0)*np.cos(0.0), np.cos(phi0)*np.sin(0.0),np.sin(phi0)]).reshape(3,1)
q01 = np.array([np.cos(phi0+0.01)*np.cos(0.0), np.cos(phi0+0.01)*np.sin(0.0),np.sin(phi0+0.01)]).reshape(3,1)
# q0rand = - np.ones((3,1)) + 2 * np.random.rand(3,1)
# q0 = q0rand / np.linalg.norm(q0rand)
w0 = np.zeros((3,1))

sol = np.zeros((2,6,N+1))
sol[0,:,0] = np.concatenate((q00,w0)).reshape(6,)
sol[1,:,0] = np.copy(sol[0,:,0])
# sol[0,:,0] = np.concatenate((q00,w0)).reshape(6,)
# sol[1,:,0] = np.concatenate((q01,w0)).reshape(6,)

qNorm = np.zeros((2,N+1))
qNorm[0,0] = np.linalg.norm(sol[0,:3,0])-1.0
qNorm[1,0] = np.linalg.norm(sol[1,:3,0])-1.0
wq = np.zeros((2,N+1))
wq[0,0] = np.dot(sol[0,3:,0],sol[0,:3,0])
wq[1,0] = np.dot(sol[1,3:,0],sol[1,:3,0])
RD = np.zeros((N+1,))
RD[0] = rDistance(sol[0,:3,0],sol[1,:3,0])


def A(x):
    q = x[:3]
    w = x[3:]
    v = gGrav*np.cross(e3,q) - d*w + k*(np.transpose(w)@q)*q
    return np.concatenate((w,v))

def vecField(x):
    Y = A(x)
    A1 = np.cross(x[:3],Y[:3])
    A2 = np.cross(x[3:],Y[:3]) + np.cross(x[:3],Y[3:])
    return np.concatenate((A1,A2))

for j in range(N):
    sol[0,:,j+1] = lie.impliemidp(A, sol[0,:,j+1], sol[0,:,0], h[j])
    sol[1,:,j+1] = cl.midptrule(vecField, sol[1,:,j], sol[1,:,0], h[j])
    # sol[1,:,j+1] = lie.impliemidp(A, sol[1,:,j+1], sol[1,:,j], h[0])
    # sol[1,:,j+1] = lie.implieeul(A, sol[1,:,j+1], sol[1,:,0], h[j])
    qNorm[0,j+1] = np.linalg.norm(sol[0,:3,j+1])-1.0
    qNorm[1,j+1] = np.linalg.norm(sol[1,:3,j+1])-1.0
    wq[0,j+1] = np.dot(sol[0,3:,j+1],sol[0,:3,j+1])
    wq[1,j+1] = np.dot(sol[1,3:,j+1],sol[1,:3,j+1])
    RD[j+1] = rDistance(sol[0,:3,j+1],sol[1,:3,j+1])



fig = plt.figure("Trajectory")
ax = fig.add_subplot(projection = '3d')
u, v = np.mgrid[0: 2 * np.pi: 30j, 0: np.pi: 20j]
x = np.cos(u) * np.sin(v)
y = np.sin(u) * np.sin(v)
z = np.cos(v)
ax.plot_surface(x, y, z, cmap = plt.cm.YlGnBu_r,alpha=0.3)
ax.plot(sol[0,0,:],sol[0,1,:],sol[0,2,:],'blue')
ax.plot(sol[1,0,:],sol[1,1,:],sol[1,2,:],'red')


fig = plt.figure("Norm")
plt.plot(h,qNorm[0,1:],'blue')
plt.plot(h,qNorm[1,1:],'red')
plt.grid()

fig = plt.figure("Orthogonal")
plt.plot(h,wq[0,1:],'blue')
plt.plot(h,wq[1,1:],'red')
plt.grid()

fig = plt.figure("Distance")
plt.plot(h,RD[1:])
plt.grid()

plt.show()

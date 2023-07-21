import numpy as np
from numpy import cos, sin
from numpy.linalg import norm
from scipy.optimize import fsolve

def mat2ang(x):
    x = x.reshape((3,3))
    theta = np.arccos(0.5 * (x[0,0]+x[1,1]+x[2,2]-1))
    return theta

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

def logso3(x):
    x = x.reshape((3,3))
    theta = mat2ang(x)
    if theta > 1e-8:
        rslt = theta/(2 * sin(theta)) * skwinv(x - np.transpose(x))
    else:
        rslt = 0.5 * (1 + 1/6 * theta**2 + 7/360 * theta ** 4) * skwinv(x - np.transpose(x))
    return rslt

def dexpinvso3(a,b):
    Theta = norm(a)
    theta = Theta * 0.5
    aa = skw(a)
    if Theta > 1e-8:
        fun = (1-theta*(np.cos(theta)/np.sin(theta)))/(Theta**2)
    else:
        fun =  1/12 + Theta**2/720 + Theta**4/30240
    operator = np.eye(3) - 0.5*aa + fun*aa@aa
    return operator@b

def tantrso3(x):
    theta = norm(x)
    if theta > 1e-14:
        fun1 = (np.cos(theta)-1)/theta**2
        fun2 = (1-(np.sin(theta)/theta))/theta**2
    else:
        fun1 = - 1.0/2 + theta**2/24 - theta**4/720
        fun2 = 1/6 - theta**2/120 + theta**4/5040
    rslt = np.eye(3) - fun1*skw(x) + fun2*skw(x)@skw(x)
    return rslt.transpose()

def taninvso3(x):
    Theta = norm(x)
    theta = Theta * 0.5
    if Theta > 1e-8:
        fun = (1-theta*(np.cos(theta)/np.sin(theta)))/(Theta**2)
    else:
        fun =  1/12 + Theta**2/720 + Theta**4/30240
    return np.eye(3) + 0.5*skw(x) + fun*skw(x)@skw(x)

def C2(x):
    a = np.copy(x[:3])
    b = np.copy(x[3:])
    Theta = norm(a)
    theta = Theta * 0.5
    if Theta > 1e-8:
        fun = (1-theta*(np.cos(theta)/np.sin(theta)))/(Theta**2)
    else:
        fun =  1/12 + Theta**2/720 + Theta**4/30240
    return 0.5*skw(b) + fun * (skw(b)@skw(a) + skw(a)@skw(b))

def dexpinvse3(a,b):
    A = taninvso3(a[:3]) # 3,3
    B = C2(a)            # 3,3
    return np.block([[A, np.zeros((3,3))],[B, A]])@b

def expse3(x):
    u = np.copy(x[:3])
    v = np.copy(x[3:])
    A11 = expso3(u)         # 3,3
    A12 = tantrso3(u)@v     # 3,
    return np.column_stack([A11, A12])

class Lie_RK:
    def __init__(self, a, b, c, stages, dt):
        self.a = a
        self.b = b
        self.c = c
        self.s = stages
        self.h = dt
    
    def explicit_integrate(self, vec_field, initial_guess):
        theta  = np.zeros(( 3, self.s))
        om     = np.zeros(( 3, self.s))
        r      = np.zeros(( 9, self.s))
        dOm    = np.zeros(( 3, self.s))
        dTheta = np.zeros(( 3, self.s))

        # first stage
        # theta[:, 0] = 0
        om    [:, 0] = np.copy(initial_guess['om0'])
        r     [:, 0] = (expso3(theta[:, 0]) @ initial_guess['r0']).reshape((9, ))
        dOm   [:, 0] = vec_field(r[:, 0], om[:, 0])
        dTheta[:, 0] = taninvso3(theta[:, 0]) @ om[:, 0]
        theta_plus   =                        self.h * self.b[0] * dTheta[:, 0]
        om_plus      = initial_guess['om0'] + self.h * self.b[0] * dOm   [:, 0]

        # following stages
        for i in range(1, self.s):
            om [:, i] = np.copy(initial_guess['om0'])
            for j in range(i):
                theta[:, i] = theta[:, i] + self.h * self.a[i,j] * dTheta[:, j]
                om   [:, i] = om   [:, i] + self.h * self.a[i,j] * dOm   [:, j]
            r     [:, i] = (expso3(theta[:, i]) @ initial_guess['r0']).reshape((9, ))
            dOm   [:, i] = vec_field(r[:, i], om[:, i])
            dTheta[:, i] = taninvso3(theta[:, i]) @ om[:, i]
            theta_plus   = theta_plus + self.h * self.b[i] * dTheta[:, i]
            om_plus      = om_plus    + self.h * self.b[i] * dOm   [:, i]

        rslt_R  = expso3(theta_plus) @ initial_guess['r0']

        return rslt_R.reshape((9, )), om_plus
    
    def stage(self, om, theta, r0, om0, vec_field):
        rslt = np.zeros(( (3+3)*self.s))
        for i in range(self.s):
            rslt[3*i:3*i+3] = np.copy(om0)
            for j in range(self.s):
                rslt[3*i:3*i+3]                   = rslt[3*i:3*i+3]                   + self.h * self.a[i,j] * vec_field(expso3(theta[3*i:3*i+3])@r0.reshape((3,3)), om[3*i:3*i+3])
                rslt[3*self.s+3*i:3*self.s+3*i+3] = rslt[3*self.s+3*i:3*self.s+3*i+3] + self.h * self.a[i,j] * (taninvso3(theta[3*i:3*i+3]) @ om[3*i:3*i+3])
        return rslt

    def implicit_integrate(self, vec_field, initial_guess):
        om_theta   = np.zeros(( (3+3)*self.s, ))
        om_plus    = np.copy(initial_guess['om0'])
        theta_plus = np.zeros(( 3, ))

        res = lambda x: - x + self.stage(x[:3*self.s], x[3*self.s:], initial_guess['r0'], initial_guess['om0'], vec_field)
        om_theta = fsolve(res, om_theta)

        for i in range(self.s):
            theta_plus = theta_plus + self.h * self.b[i] * om_theta[3*self.s+3*i:3*self.s+3*i+3]
            om_plus    = om_plus    + self.h * self.b[i] * om_theta[3*i:3*i+3]

        rslt_R  = expso3(theta_plus) @ initial_guess['r0']

        return rslt_R.reshape((9, )), om_plus
    
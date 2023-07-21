import numpy as np
import datetime
import pickle
import integrLie as lie
import time
import winsound
from alive_progress import alive_bar

# define problem 
def rhs(rotMat, om, par):
    # rslt = np.zeros((9+3, 1))
    # rslt[:9] = np.reshape(rotMat @ lie.skw(om), (9, 1))
    # rslt[9:] = - np.linalg.inv(par['inertia']) @ (par['damping'] * om + par['stiffness'] * lie.logso3(rotMat))
    rslt = - np.linalg.inv(par['inertia']) @ (par['damping'] * om + par['stiffness'] * lie.logso3(rotMat))
    return rslt

def main(input):
    # initialize problem
    parameter = dict()
    parameter['mass']      = 15
    parameter['inertia']   = parameter['mass'] * np.eye(3) # inertia of a sphere
    parameter['damping']   = 0.5
    parameter['stiffness'] = 1
    parameter['direction'] = np.array([0., 0., 1.]) # normalized axis of rotation
    parameter['alpha']     = 0
    parameter['velocity']  = 10.
    parameter['r0']        = lie.expso3(parameter['alpha'] * parameter['direction'])
    parameter['om0']       = parameter['velocity'] * parameter['direction']

    my_IC = dict()
    my_IC['r0']     = lie.expso3(parameter['alpha'] * parameter['direction'])
    my_IC['om0']    = parameter['velocity'] * parameter['direction']
    # my_IC['om0']    = np.array([np.sqrt(2)/2, np.sqrt(2)/2, 0.])

    # time interval
    parameter['t0'] = 0.
    parameter['te'] = 1.
    parameter['dt'] = input
    time_span = np.arange(parameter['t0'], parameter['te'], parameter['dt'])
    n_steps   = len(time_span)

    # numerical parameters
    parameter['stages'] = 1
    parameter['a'] = np.array([[0.5]])
    parameter['b'] = np.array([1.])
    parameter['c'] = np.sum(parameter['a'], 1)
    # parameter['stages'] = 2
    # parameter['a'] = np.array([[0., 0.], [0.5, 0.]])
    # parameter['b'] = np.array([0., 1.])
    # parameter['c'] = np.sum(parameter['a'], 1)
    # parameter['stages'] = 3
    # parameter['a'] = np.array([[0., 0., 0.], [0.5, 0., 0.], [-1., 2., 0.]])
    # parameter['b'] = np.array([1/6, 2/3, 1/6])
    # parameter['c'] = np.sum(parameter['a'], 1)

    my_integrator = lie.Lie_RK(parameter['a'], parameter['b'], parameter['c'], parameter['stages'], parameter['dt']) # object of class Lie_RK

    # solution
    r    = np.zeros((9, n_steps))
    om   = np.zeros((3, n_steps))
    r    [:, 0] = np.copy(my_IC['r0'].reshape((9, )))
    om   [:, 0] = np.copy(my_IC['om0'])

    # time integration
    with alive_bar(n_steps-1, bar = 'smooth') as bar:
        for i in range(1, n_steps):
            r[:, i], om[:, i]= my_integrator.implicit_integrate(lambda x, y: rhs(x, y, parameter), my_IC)
            my_IC['r0']      = np.copy(r[:, i].reshape((3,3)))
            my_IC['om0']     = np.copy(om[:, i])
            bar()

    now = datetime.datetime.now().strftime('%Y%m%dT%H%M%S')

    fileName = open('extras/test_eq/out/'+now+'_par.pkl', 'wb')
    pickle.dump(parameter, fileName, -1)
    fileName.close()

    np.save('extras/test_eq/out/'+now+'_sols', np.block([[r], [om]]))

if __name__ == '__main__':
    for counter,timeStepSize in enumerate(np.logspace(-4,-3,6)):
        print('Start simulation '+str(counter+1)+'.')
        start = time.time()
        main(timeStepSize)
        end = time.time()
        time.sleep(1)
        total_time = float(time.time() - start)
        if total_time > 60:
            seconds = int(total_time % 60)
            minutes = total_time // 60
            total_time = str(int(minutes)) + ' minute(s) ' + str(int(seconds)) + ' second(s)'
        else:
            total_time = str(int(total_time)) + ' second(s)'
        print('Completed simulation '+str(counter+1)+' in '+total_time)
    winsound.PlaySound("SystemExit", winsound.SND_ALIAS)
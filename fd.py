import numpy as np
import matplotlib.pyplot as plt
import random

# Parameters
  # Data from https://doi.org/10.1016/0895-7177(93)90025-T
a = 0.1
b = 0.9
gamma = 2
d = 10

Nx = 15
Nt = 10000
L = 4
T = 5


x = np.linspace(0, L, Nx+1) # mesh points in space
dx = x[1] - x[0]
t = np.linspace(0, T, Nt+1) # mesh points in time
dt = t[1] - t[0]

mu = dt / dx**2

if (mu > 0.5):
    print(mu)
    exit


u = np.zeros((Nx+1, Nt+1)) # unknown u
v = np.zeros((Nx+1, Nt+1)) # unknown v

def f(u, v): # define f
    return a - u + u**2*v 

def g(u, v): # define g
    return b - u**2*v

def I(x):
    return 0.5*random.random()

def J(x):
    return 0.75*random.random()

# Print iterations progress
def printProgressBar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█', printEnd = "\r"):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
        printEnd    - Optional  : end character (e.g. "\r", "\r\n") (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print(f'\r{prefix} |{bar}| {percent}% {suffix}', end = printEnd)
    # Print New Line on Complete
    if iteration == total: 
        print()

printProgressBar(0, (Nt-1), prefix = 'Progress:', suffix = 'Complete', length = 50)
# Set initial condition u(x,0) = I(x) / v(x, 0) = J(x)
for i in range(0, Nx+1):
    u[i, 0] = I(x[i])
    v[i, 0] = J(x[i])
    
for m in range(0, Nt-1):
    # Insert boundary conditions
    u[0,m] = u[1,m]; u[Nx,m]= u[Nx-1,m]
    v[0, m] = v[1,m]; v[Nx, m] = v[Nx-1,m]
    
    # Compute u at inner mesh points
    for j in range(1, Nx):
        u[j, m+1] = u[j,m] + mu*(u[j-1,m] - 2*u[j,m] + u[j+1,m]) + dt*gamma*f(u[j,m], v[j,m])             
        # forever let it be known, Luna forgot that it was dt * gamma * f
        v[j, m+1] = v[j,m] + d*mu*(v[j-1,m] - 2*v[j,m] + v[j+1,m]) + dt*gamma*g(u[j,m], v[j,m])
    printProgressBar(m, (Nt-1), prefix = 'Progress:', suffix = 'Complete', length = 50)
                
# print(u)
# print(v)

# plot
fig = plt.figure()
plt.imshow(u.T, extent=[-4,4,-1,1], aspect=4)


fig2 = plt.figure()
plt.imshow(v.T, extent=[-4,4,-1,1], aspect=4)

fig3 = plt.figure()
# print(u[:,Nt])
plt.plot(t, u[Nx,:])

plt.show()

# fig = plt.figure()
# ax2 = plt.axes()
# ax2.plot(X, T)
# plt.show()
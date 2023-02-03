import numpy as np
import matplotlib.pyplot as plt

# Parameters
a = 0.1
b = 0.9
Nx = 100
Nt = 50
gamma = 0.5
d = 1
L = 1
T = 1


x = np.linspace(0, L, Nx+1) # mesh points in space
dx = x[1] - x[0]
t = np.linspace(0, T, Nt+1) # mesh points in time
dt = t[1] - t[0]

mu = dt / dx**2

u = np.zeros((Nx+1, Nt+1)) # unknown u
v = np.zeros((Nx+1, Nt+1)) # unknown v

def f(u, v):
    return a - u + u**2*v 

def g(u, v):
    return b - u**2*v

def I(x):
    return 1

def J(x):
    return 1

# Set initial condition u(x,0) = I(x) / v(x, 0) = J(x)
for i in range(0, Nx+1):
    u[i, 0] = I(x[i])
    v[i, 0] = J(x[i])
    
for m in range(0, Nt-1):
    # Insert boundary conditions
    u[0,m] = 0; u[Nx,m]= 0
    v[0, m] = 0; v[Nx, m] = 0
    
    # Compute u at inner mesh points
    for j in range(1, Nx):
        u[j, m+1] = u[j,m] + mu*(u[j-1,m] - 2*u[j,m] + u[j+1,m]) + dt*gamma*f(u[j,m], v[j,m]) 
        # forever let it be known, Luna forgot that it was dt * gamma * f
        v[j, m+1] = v[j,m] + d*mu*(u[j-1,m] - 2*u[j,m] + u[j+1,m]) + dt*gamma*g(u[j,m], v[j,m])
        print(u[j,m])
            
# print(u)

# plot
fig = plt.figure()
plt.imshow(u, interpolation='nearest')
plt.show()

# fig = plt.figure()
# ax2 = plt.axes()
# ax2.plot(X, T)
# plt.show()
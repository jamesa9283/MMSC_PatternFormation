import numpy as np
import matplotlib.pyplot as plt

# Parameters
M = 5
a = 0.1
b = 0.9
Nx = 40
Nt = 100
gamma = 0.1
d = 1
L = 1
T = 0.06


x = np.linspace(0, L, Nx+1) # mesh points in space
dx = x[1] - x[0]
t = np.linspace(0, T, Nt+1) # mesh points in time
dt = t[1] - t[0]

mu = dt / dx**2

u = np.zeros((Nx+1, Nt+1)) # unknown u

def f(u, x, t):
    return 0

def I(x):
    return 0.1

# Set initial condition u(x,0) = I(x)
for i in range(0, Nx+1):
    u[i, 0] = I(x[i])
    
for m in range(0, Nt-1):
    # Insert boundary conditions
    u[0,m] = 0; u[Nx,m]=0
    
    # Compute u at inner mesh points
    for j in range(1, Nx):
        u[j, m+1] = u[j,m] + mu*(u[j-1,m] - 2*u[j,m] + u[j+1,m]) + gamma*f(u[j,m], x[j], t[m])
        print(u[j,m])
            
# print(u)

# plot
fig = plt.figure()
ax = plt.axes(projection='3d')

X, T = np.meshgrid(x, t)
ax.plot_surface(T, X, u.T, linewidth=2.0)

fig = plt.figure()
ax2 = plt.axes()
ax2.plot(X, T)
plt.show()
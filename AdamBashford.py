import numpy as np 
import matplotlib.pyplot as plt
import math
from mpl_toolkits.mplot3d import Axes3D

### Using AB2 BM2 numerical methods
# Parameters
# Data from https://doi.org/10.1016/0895-7177(93)90025-T
a = 0.1 #also in function
b = 0.9 #also in function
Nx = 150
Nt = 1000
gamma = 2
d = 10
L = 20
T = 6.4

#create meshgrid
x = np.linspace(0, L, Nx+1) # mesh points in space
dx = x[1] - x[0]
t = np.linspace(0, T, Nt+1) # mesh points in time
dt = t[1] - t[0]

mu = dt / dx**2

#and empty arrays for U and V 
U = np.zeros([Nx+1,1]) # unknown u
V = np.zeros([Nx+1,1]) # unknown v

def f(u, v): # define f
    return a - u + (u**2)*v 

def g(u, v): # define g
    return b - (u**2)*v

#creating tridiag matrices
Identity_BCs = np.diagflat([-1]+[0 for i in range(Nx+1-2)],1) +\
np.diagflat([0 for i in range(Nx+1-2)]+[-1],-1) +\
np.diagflat([1 for i in range(Nx+1)]) # WHY?!?!?!?

I = np.diagflat([1 for i in range(Nx+1)])


d2x = np.diagflat([0]+[-2 for i in range(Nx+1-2)]+[0]) +\
        np.diagflat([0]+[1 for i in range(Nx+1-2)],1) +\
        np.diagflat([1 for i in range(Nx+1-2)]+[0],-1)

#and setting the left hand sides of the equations
LHS_U = Identity_BCs-mu*d2x
LHS_V = Identity_BCs-d*mu*d2x

#setting initial values of U and V 
for i in range(0, Nx+1):
    U[i,0] = 0.75*np.random.rand(1)
    V[i,0] = 0.7*np.random.rand(1)

U_record = np.zeros((Nt+1, Nx+1)) # initialise U
V_record = np.zeros((Nt+1, Nx+1)) # intiialise V

U_record[0] = U[:,0] 
V_record[0] = V[:,0]

def f_vec(U,V):
    if len(U) and len(V) != Nx+1:
        print("Size Error")
    vec = np.zeros([Nx+1,1])
    for i in range(0,len(U)):
        vec[i] = f(U[i],V[i])
    return vec

def g_vec(U,V):
    if len(U) and len(V) != Nx+1:
        print("Size Error")
    vec = np.zeros([Nx+1,1])
    for i in range(0,len(V)):
        vec[i] = g(U[i],V[i])
    return vec

for timestep in range(1,3):
    RHS_U = U+dt*gamma*(f_vec(U,V))
    RHS_V = V+dt*gamma*(g_vec(U,V))
    RHS_U[0] = 0
    RHS_U[-1] = 0
    RHS_V[0] = 0
    RHS_V[-1] = 0
    U_new = np.linalg.solve(LHS_U, RHS_U)
    V_new = np.linalg.solve(LHS_V, RHS_V)
    
    U = U_new
    V = V_new
    
    U_record[timestep] = U[:,0]
    V_record[timestep] = V[:,0]

#using ab2am2 methods 
ab2_LHS_U = 2*Identity_BCs-mu*d2x
ab2_LHS_V = 2*Identity_BCs-d*mu*d2x

def ab2_f_vec(U1,U2,V1,V2):
    if len(U1) and len(V1) != Nx+1:
        print("Size Error")
    vec = np.zeros([Nx+1,1])
    for i in range(0,len(U1)):
        vec[i] = 3*f(U1[i],V1[i])-f(U2[i],V2[i])
    return vec

def ab2_g_vec(U1,U2,V1,V2):
    if len(U1) and len(V1) != Nx+1:
        print("Size Error")
    vec = np.zeros([Nx+1,1])
    for i in range(0,len(U1)):
        vec[i] = 3*g(U1[i],V1[i])-g(U2[i],V2[i])
    return vec

for i in range(3,Nt+1):
    # TODO Correct Formula 2
    #using formula 1 
    #ab2_RHS_U = (mu*d2x + 2*Identity_BCs)*U_record[i-1].reshape(len(U_record[0]),1) +\
    #            dt*gamma*(ab2_f_vec(U_record[i-1],U_record[i-2],V_record[i-1],V_record[i-2]))
    
    #using formula 2
    ab2_RHS_U = (mu*(U_record[i-1].reshape(len(U_record[0]),1) - 2*U_record[i-2].reshape(len(U_record[0]),1) + U_record[i-3].reshape(len(U_record[0]),1)) +\
                2*Identity_BCs)* U_record[i-1].reshape(len(U_record[0]),1) +\
                dt*gamma*(ab2_f_vec(U_record[i-1],U_record[i-2],V_record[i-1],V_record[i-2]))
                
    ab2_RHS_U2 = (2*Identity_BCs + (1 - mu)*U_record[i-1].reshape(len(U_record[0]),1) + (2*mu-1)*U_record[i-2].reshape(len(U_record[0]),1) - mu*U_record[i-3].reshape(len(U_record[0]),1)) * U_record[i-1].reshape(len(U_record[0]),1) +\
                dt*gamma*(ab2_f_vec(U_record[i-1],U_record[i-2],V_record[i-1],V_record[i-2]))

    ab2_LHS_U2 = (2*Identity_BCs - (1 - mu)*U_record[i-1].reshape(len(U_record[0]),1) - (2*mu-1)*U_record[i-2].reshape(len(U_record[0]),1) + mu*U_record[i-3].reshape(len(U_record[0]),1))

    ab2_RHS_U[0] = 0
    ab2_RHS_U[-1] = 0
    U_new = np.linalg.solve(ab2_LHS_U2, ab2_RHS_U2)

    #using formula 1
    #ab2_RHS_V = (mu*d2x + 2*Identity_BCs)*V_record[i-1].reshape(len(V_record[0]),1) +\
    #            dt*d*gamma*(ab2_g_vec(U_record[i-1],U_record[i-2],V_record[i-1],V_record[i-2]))

    #using formula 2
    ab2_RHS_V = (mu*(V_record[i-1].reshape(len(V_record[0]),1) - 2*V_record[i-2].reshape(len(V_record[0]),1) + V_record[i-3].reshape(len(V_record[0]),1)) +\
                2*Identity_BCs)* V_record[i-1].reshape(len(V_record[0]),1) +\
                dt*gamma*(ab2_g_vec(U_record[i-1],U_record[i-2],V_record[i-1],V_record[i-2]))
                
    ab2_RHS_V2 = (2*Identity_BCs + (1 - mu)*U_record[i-1].reshape(len(U_record[0]),1) + (2*mu-1)*U_record[i-2].reshape(len(U_record[0]),1) - mu*U_record[i-3].reshape(len(U_record[0]),1)) * U_record[i-1].reshape(len(U_record[0]),1) +\
                dt*d*gamma*(ab2_g_vec(U_record[i-1],U_record[i-2],V_record[i-1],V_record[i-2]))


    ab2_RHS_V[0] = 0
    ab2_RHS_V[-1] = 0
    V_new = np.linalg.solve(ab2_LHS_U2, ab2_RHS_V2)

    U_record[i] = U_new[:,0]
    V_record[i] = V_new[:,0]

#plotting results 
#2D
plt.title("U")
plt.plot(x,V_record[-1])

#3D
fig = plt.figure()
ax = plt.axes(projection='3d')

X, T = np.meshgrid(x, t)
ax.plot_surface(T,X, V_record, linewidth=2.0)

plt.show()

#heat map
#plot the 2D colour map of U 
plt.figure()
plt.imshow((U_record).T, cmap='plasma', origin='lower', extent=[0, 15, 0, 6.4], aspect='auto')
plt.colorbar()
plt.title("U colour map")
plt.xlabel("x")
plt.ylabel("t")
plt.show()

 
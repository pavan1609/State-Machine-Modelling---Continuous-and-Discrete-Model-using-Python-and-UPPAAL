 #Python program to model flow in pipe, where viscosity of the fluid is temperature dependent, and the Temperature varies with Radius as given µ(T) = µ0 (1 + α (T − T0)) T(r) = T0 + βr2


import numpy as np
import matplotlib.pyplot as plt

#Parameters
mu = 0.01   #Viscosity of fluid
R = 0.05    #Radius of Pipe
N = 100     #Grid Points
delta_r = R / N
mu0 = 0.05   #Base Viscosity
alpha = 0.08 #Temperature Expansion Coefficient 
beta = 5     #constant 
T0 = 273     # Reference temperature

#Initialisation for Velocity Vector u_z, Cofficeint Matrix A for equation, Constant Vector b
u_z = np.zeros(N)
A = np.zeros((N, N))
b = np.zeros(N)

#Creating radius values into array at different points
r_values = np.linspace(0, R, N)
#r_values = np.zeros(N)
#for i in range(0, N):
#    r_values[i] = i * delta_r
    
# Initial Boundary conditions
A[0, 0] = 1.0   # Boundary condition at r = 0
A[0, 1] = -1.0  # Boundary condition at r = 0
b[0] = 0        # Boundary condition at r = 0

# Final Boundary Conditions
A[-1, -1] = 1  #Boundary Condition at r = R
b[-1] = 0      #Boundary Condition at r = R

#Forming Coefficient matrix from the equation using for loop
for i in range(1, N-1):
    T = T0 + (beta * r_values[i]**2)
    mu = mu0 * (1 + (alpha * (T - T0)))
    a = mu / (delta_r**2)
             
    A[i, i-1] = (mu / (delta_r**2)) * r_values[i]
    A[i, i] = -(mu / (delta_r**2)) * ((2 * r_values[i]) + delta_r)
    A[i, i+1] = (mu / (delta_r**2)) * (r_values[i] + delta_r)
    
#Forming Constant Vector matrix from equation using for loop
for i in range(1, N-1):
    b[i] = -100 * r_values[i]
    
#Function to solve the linear equations
u_z = np.linalg.solve(A, b)
#print(u_z)
#print(r_values)

#plotting the graph of velocity values at each radius point
plt.plot(r_values,u_z)
plt.xlabel('radius of pipe (in meters)',fontsize=12)
plt.ylabel('Velocity of fluid  U_z (in m/s)',fontsize=12)
plt.title('Velocity profile at each Radius Point',fontsize=12)
plt.grid()
plt.show()

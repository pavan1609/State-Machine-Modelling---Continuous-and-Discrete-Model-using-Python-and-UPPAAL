#Python program to solve the system of linear equations obtained by converting the Equation into difference equation. Create a 3D plot of velocity profile as a function of both radial (r) and axial (z) coordinates.



import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

Nr = 100
Nz = 100
R = 0.05
L = 1.0
delta_r = R / (Nr-1)
delta_z = L / (Nz-1)
print (delta_r)
print (delta_z)
dpdz = -100
A = np.zeros((Nr * Nz, Nr * Nz))  # A matrix initialized to Zero
B = np.full((Nr * Nz),dpdz)
r = np.linspace(0, R, Nr)
z = np.linspace(0, L, Nz)
Reversed_U_z = np.zeros([Nr,Nz])

mu = 0.01
gamma = 0.25
def Index_Function(a,b):
    index = a + b *Nr
    return index

def Reverse_Index(index):
    Y_Index = int(index / Nr)
    X_Index = int(index % Nr)
    return X_Index,Y_Index
    
for j in range(Nz):
    for i in range(Nr):
        radius = i * delta_r
        #print (radius)
        if j == 0:  # Boundary condition at z = 0
            A[Index_Function(i,j), Index_Function(j,i)] = 1
            B[Index_Function(i,j)] = 15
        elif j == (Nz - 1):  # Boundary condition at z = L
            A[Index_Function(i,j), Index_Function(j,i)] = 1
            B[Index_Function(i,j)] = 15
        elif i == 0:  # Boundary condition at r = 0 (axis)
            A[Index_Function(i,j), Index_Function(j,i+1)] = 1
            A[Index_Function(i,j), Index_Function(j,i)] = -1
            B[Index_Function(i,j)] = 0

        elif i == (Nr - 1):  # Boundary condition at r = R (disk)
            A[Index_Function(i,j), Index_Function(j,i)] = 1
            B[Index_Function(i,j)] = 0
        else:
            A[Index_Function(i,j), Index_Function(j,i+1)] = ((mu) / (delta_r ** 2 )) + ((mu) / (radius * delta_r))
            A[Index_Function(i,j), Index_Function(j,i)] = -1 * (((2 * mu) / (delta_r ** 2)) + (mu / (radius * delta_r)) + (2 * gamma / (delta_z ** 2)))
            A[Index_Function(i,j), Index_Function(j,i-1)] = mu / (delta_r ** 2)
            A[Index_Function(i,j), Index_Function(j+1,i)] = gamma / (delta_z ** 2)
            A[Index_Function(i,j), Index_Function(j-1,i)] = gamma / (delta_z ** 2)

# Solve the linear system
U_z = np.linalg.solve(A, B)

# Reshape velocity array for plotting
for a, b in enumerate (U_z):
    Reversed_U_z_Index = Reverse_Index(a)
    Reversed_U_z[Reversed_U_z_Index[0]][Reversed_U_z_Index[1]] = b
#U_z = U_z.reshape((Nr, Nz)).T

# Create an inverted parabolic velocity profile
#for i in range(Nz):
#    U_z[i, :] = dpdz / (2 * mu) * z[i] * (L - z[i])

# Plotting
fig = plt.figure(figsize=(10, 6))
ax = fig.add_subplot(111, projection='3d')
R, Z = np.meshgrid(r, z)
ax.plot_surface(Z, R, Reversed_U_z, cmap='viridis')
ax.set_title('Parabolic Velocity profile')
ax.set_xlabel('Z-axis')
ax.set_ylabel('R-axis')
ax.set_zlabel('Velocity U_z')
plt.show()

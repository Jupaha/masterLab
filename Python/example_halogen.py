"""compute eigenvalues of linearised halogen lamp system"""
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

# parameters for the halogen lamp system
R_c = 0.37
L = 15e-6
c_w = 7e-3
T_u = 295
b = 1.5e-12
k = 0.9

# parameters for general equation form
params = np.zeros(6)
params[0] = -R_c / (L * T_u**k)
params[1] = k
params[2] = 1 / L
params[3] = -b / c_w
params[4] = T_u
params[5] = R_c / (c_w * T_u**k)

# non-linear differential system has the folling general form

# def f0(xs, u, p):
#     return p[0] * xs[1]**p[1] * xs[0] + p[2] * u

# def f1(xs, u, p):
#     return p[3] * (xs[1] - p[4])**4 + p[5] * xs[1]**p[1] * xs[0]**2

# no. of discretization elements per variable
n = 20
# current values for linearization
x0 = np.linspace(0, 100, n)
# temperature values for linearization
x1 = np.linspace(300, 10000, n)
# matrix with compinations of current and temperature (x0/x1)
x = np.zeros((n, n, 2))
# copmute x values
for i in range(n):
    for j in range(n):
        x[i][j] = [x0[i], x1[j]]

def df0_dx0(xs, p):
    """computes df0/dx0 with x values xs and parameters p"""
    return p[0] * xs[1]**p[1]

def df0_dx1(xs, p):
    """computes df0/dx1 with x values xs and parameters p"""
    return p[0] * p[1] * xs[1]**(p[1] - 1) * xs[0]

def df1_dx0(xs, p):
    """computes df1/dx0 with x values xs and parameters p"""
    return p[5] * xs[1]**p[1] * 2 * xs[0]

def df1_dx1(xs, p):
    """computes df1/dx1 with x values xs and parameters p"""
    return 4 * p[3] * (xs[1] - p[4])**3 + p[5] * p[1] * xs[1]**(p[1] - 1) * xs[0]**2

A = np.zeros((n, n, 2, 2))
eigVal = np.zeros((n, n, 2))

for i in range(n):
    for j in range(n):
        A[i, j] = np.array(
            [[df0_dx0(x[i, j], params),
              df0_dx1(x[i, j], params)],
             [df1_dx0(x[i, j], params),
              df1_dx1(x[i, j], params)]])
        eigVal[i, j], _ = np.linalg.eig(A[i, j])

print('Max. Eigenwert: {:.2E}'.format(np.amax(eigVal)))
print('Min. Eigenwert: {:.2E}'.format(np.amin(eigVal)))
print('Erforderliche Schrittweite: {:.2E}'.format(1/(2*np.amax(eigVal))))

# plot first eigenvalues dependent of x
fig1 = plt.figure()
ax1 = fig1.gca(projection='3d')
X, Y = np.meshgrid(x0, x1, indexing='ij')
Z = eigVal[:, :, 0]
surf = ax1.plot_surface(
    X, Y, Z, cmap=plt.get_cmap('viridis'), linewidth=0.1, antialiased=True)
ax1.set_xlabel('Strom in A')
ax1.set_ylabel('Temperatur in K')
ax1.set_zlabel('Eigenwert 0')

# plot first eigenvalues dependent of x
fig2 = plt.figure()
ax2 = fig2.gca(projection=Axes3D.name)
Z = eigVal[:, :, 1]
surf = ax2.plot_surface(
    X, Y, Z, cmap=plt.get_cmap('viridis'), linewidth=0.1, antialiased=True)
ax2.set_xlabel('Strom in A')
ax2.set_ylabel('Temperatur in K')
ax2.set_zlabel('Eigenwert 1')

plt.show()

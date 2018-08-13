"""compute eigenvalues of non-linear DGL system"""
import itertools as it
import numpy as np
import sympy as sp

def gen_states(N):
    """generate symbolic variables for N states"""
    return sp.symbols(' '.join(['x' + str(n) for n in range(N)]))


def gen_xe(N, n, r):
    """generate array with all combinations of x[i]"""
    # array to store x ranges
    xs = np.zeros((N, n))
    for i in range(N):
        xs[i] = np.linspace(r[i][0], r[i][1], n)
    return list(it.product(*xs))


def comp_eigenvalues(f, x, r, n):
    """compute eigenvalues of non linear DGL system"""

    # dimension of problem
    N = len(f)

    # system matrix is the jacobian matrix df_i / d_xi
    A = f.jacobian(x)

    # generate array with all combinations of x[i]
    xe = gen_xe(N, n, r)

    # solve A for all combinations of xe entrys and compute eigenvalues
    A_all = np.zeros((n**N, N, N), dtype=complex)
    eigV = np.zeros((n ** N, N), dtype=complex)
    for i in range(n ** N):
        A_temp = A
        for j in range(N):
            A_temp = A_temp.subs(x[j], xe[i][j])
        A_all[i] = np.array(A_temp).astype(np.float)
        eigV[i], _ = np.linalg.eig(A_all[i])

    print(eigV)
    # print some information

    if np.amin(eigV.real) < 0:
        print('System is stable')
    else:
        print('System is not stable')

    if np.iscomplex(eigV).any():
        print('System has complex eigenvalues:')
        eigV_comp = [ev for ev in eigV if np.iscomplex(ev).any()]
        print(eigV_comp)
        print('Max. eigenvalue: {:.2E}'.format(np.amax(eigV)))
        print('Min. eigenvalue: {:.2E}'.format(np.amin(eigV)))
    else:
        print('Max. eigenvalue: {:.2E}'.format(np.amax(eigV.real)))
        print('Min. eigenvalue: {:.2E}'.format(np.amin(eigV.real)))

    print('minimum step size: {:.2E}'.format(1 / (2 * np.amax(np.absolute(eigV)))))
    return eigV


def main():
    """main function with example how to use comp_eigenvalues()"""

    # parameters for a simplyfied halogen lamp system
    R_c = 0.37 # electrical resistance @ T_u in ohm
    L = 15e-6 # inductance in H
    c_w = 7e-3 # heat capacity of filament in J/K
    T_u = 295 # ambient temperature in K
    b = 1.5e-12 # coefficient to compute radiant power in W/K^4
    k = 0.9 # exponent to compute electrical resistance at hot temperatures

    # parameters for general equation form
    p = np.zeros(6)
    p[0] = -R_c / (L * T_u**k)
    p[1] = k
    p[2] = 1 / L
    p[3] = -b / c_w
    p[4] = T_u
    p[5] = R_c / (c_w * T_u**k)

    # create state vector
    x = gen_states(2)

    # define DGLs
    f = sp.Matrix([
        p[0] * x[1]**p[1] * x[0],
        p[3] * (x[1] - p[4])**4 + p[5] * x[1]**p[1] * x[0]**2
    ])

    # define ranges of x[i]
    r = [(0, 100), (300, 10000)]

    # define steps in ranges
    n = 20

    # compute eigenvalues
    comp_eigenvalues(f, x, r, n)


if __name__ == '__main__':
    main()

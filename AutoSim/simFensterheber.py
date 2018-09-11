import numpy as np
import scipy.io as sio
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from math import sqrt, pi, sin
mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.unicode'] = True
mpl.rcParams['font.family'] = 'Minion Pro'
mpl.rcParams['font.size'] = 12
mpl.rcParams['xtick.labelsize'] = 12
mpl.rcParams['ytick.labelsize'] = 12
mpl.rcParams['lines.linewidth'] = 0.8
mpl.rcParams['figure.facecolor'] = 'white'
mpl.rcParams['figure.figsize'] = [29.7/2.54, 21/2.54]

def sign(x):
    if x == 0:
        return 0
    if x < 0:
        return -1
    if x > 0:
        return 1
    return 2


def progress(step, max_steps):
    percents_last = int(100 * (step-1 / max_steps))
    percents = int(100 * (step / max_steps))
    if step > 1:
        if percents != percents_last:
            print('simulated {}%    '.format(percents), end='\r')
    else:
        print('simulated {}%    '.format(percents), end='\r')


start_position = 1 # 1: upper stop, 0: lower stop

U_B = 12.5  # batterie voltage
R = 0.875    # motor resistance (0.75)
L = 1500e-6 # motor inductance
k = 0.021   # motor constant
J = 1.5e-5  # moment of inertia
b = 2e-5    # viscouss friction motor (2)

# mechanical parameters of window regulator
m = 8.26   # equivalent window mass
c = 45e3    # spring constant of mech. system
d = 310    # damping constant of mech. system
eta_G = 0.395 # worm gear efficiency
r = 0.025   # radius of cable drum (0.025)
G = 73      # worm gear translation ratio

# parameter for variable friction
H = G * (k**2 * eta_G + b * R) / (k * r) # constant for easier use
I_0 = 1.2
I_up_min = 6
I_up_max = 7.5
I_down_min = 3
I_down_max = 4.5
I_deltaminmax = I_up_max - I_up_min
I_m = I_down_min - I_0
I_delta_updown = I_up_min - I_down_min
I_n = (I_up_min + I_up_max + I_down_min + I_down_max) / 4
F_r_min = (I_down_min - I_0) * H
F_r_max = I_deltaminmax * H + F_r_min

# misc.
g = 9.81      # gravity
# phi0 = -58.4  # angle of loose cable
if start_position:
    x_max = 0.51 # upper stop of window (0.545 down)
    x_start = x_max
else:
    x_max = 0.52  # upper stop of window (0.545 down)
    x_start = 0
phi_max = x_max * (G/r)
phi_loose = 61


# ripple parameter
a = 0.595                   # current influence factor
I_RS = 0.0726               # peak value current ripple at fast speed
n_x = 20                     # corner speed
alpha_start = 0.45 * 2 * pi # start angle
z = 10                      # commutator segments

iced = False
F_break = H * 0.78 * U_B / R

blocked = False
h_B = 0.1


def u_pwm(t, f, pw, amp):
    T = 1 / f
    if (t % T) / T < pw:
        return amp
    else:
        return 0
    

def u_test_pwm(t):
    if t < 0.1:
        return 0
    if t >= 0.1 and t < 5:
        return u_pwm(t, 25e3, 0.8, -12)



# # returns a voltage curve for testing the simulation
# def U_test(t):
#     if t < 0.1:
#         return 0
#     if t >= 0.1 and t < 2:
#         return -U_B
#     if t >= 2 and t < 2.5:
#         return 0
#     if t >= 2.5 and t < 4.5:
#         return -U_B
#     if t >= 4.5 and t < 5:
#         return 0
#     if t >= 5 and t < 9:
#         return U_B
#     if t >= 9 and t < 9.00005:
#         return -U_B
#     if t >= 9.00005:
#         return -U_B
#     return 0

# ### VARIABLES ###
t = 0
u = 0
di = 0
i_M = 0
i_R = 0
i = 0
dw = 0
w = 0
n = 0
phi = (x_start * (G / r))
ddx2 = 0
dx2 = 0
x2 = x_start
F_d = 0
F_c = 0
F_r = 0
M_A = 0
M_d = 0
M_c = 0
M_b = 0
M_L = 0


def set_phi():
    global phi, mech
    if x2 == x_max:
        phi = (x_start * (G / r)) + phi_loose
        mech['phi_switch'] = (x_max * (G / r))
        return
    if x2 == 0:
        phi = -phi_loose
        mech['phi_switch'] = 0
        return
    phi = (x_start * (G / r))
    return


def i_idle ():
    return u / (R + (eta_G * k ** 2) / b)


# computes the overlayed current ripple from the commutation
def i_ripple():
    return (((a * abs(i_M) + 1) * I_RS * n) / sqrt(n * n + n_x ** 2)) * sin((phi + alpha_start) * z)


mech = {
    'direction': [0, 0],
    'state': 5,
    'phi_switch': (x_max * (G / r)),
    'block_window': False,
    'block_motor': False,
    'loose': False,
    'self_lock': False
}


def signbit(num):
    if num > 0:
        return False
    else:
        return True


def mech_state():
    global mech, phi, iced
    if mech['state'] == 1:
        if abs(M_A) < abs(M_L) and signbit(M_L) != signbit(w):
            mech['state'] = 3
            return
        if x2 > x_max:
            mech['state'] = 4
            return
        if blocked and x2 > x_max - h_B:
            mech['state'] = 4
            return
        return
    if mech['state'] == 2:
        if x2 < 0:
            mech['state'] = 4
            return
        if abs(M_A) < abs(M_L) and signbit(M_A) != signbit(w):
            mech['state'] = 3
            return
        return
    if mech['state'] == 3:
        if abs(M_A) > abs(M_L) and signbit(M_A) == signbit(M_L) and i > 0:
            mech['state'] = 1
            return
        if signbit(M_A) != signbit(M_L) and M_A > 0:
            mech['state'] = 1
            return
        if abs(M_A) > abs(M_L) and signbit(M_A) == signbit(M_L) and i < 0:
            mech['state'] = 2
            return
        if abs(M_A) > abs(M_L) and signbit(M_A) != signbit(M_L) and i < 0:
            mech['state'] = 2
            return
        return
    if mech['state'] == 4:
        if abs(u) < 1 and abs(i_M) < 0.1:
            if x2 <= 0:
                phi = -phi_loose
                mech['phi_switch'] = 0
                mech['state'] = 5
                return
            elif x2 >= x_max:
                phi = (x_max * (G / r)) + phi_loose
                mech['phi_switch'] = (x_max * (G / r))
                mech['state'] = 5
                return
            elif blocked and x2 >= x_max - h_B:
                phi = ((x_max-h_B) * (G / r)) + phi_loose
                mech['phi_switch'] = ((x_max-h_B) * (G / r))
                mech['state'] = 5
                return
        elif abs(M_A) < abs(M_L) and signbit(M_A) == signbit(M_L) and M_L != 0 and u == 0 and abs(i) < 0.1:
            mech['state'] = 3
            return
        elif signbit(M_A) != signbit(M_L) and M_A > 0:
            mech['state'] = 1
            return
        elif x2 > 0 and x2 < x_max and dx2 > 0:
            mech['state'] = 1
            return
        elif x2 > 0 and x2 < x_max and dx2 < 0:
            mech['state'] = 2
            return
        return
    elif mech['state'] == 5:
        if w > 0 and phi > mech['phi_switch']:
            mech['state'] = 1
            return
        if w < 0 and phi < mech['phi_switch']:
            if iced:
                mech['state'] = 6
                return
            mech['state'] = 2
            return
        return
    elif mech['state'] == 6:
        if F_c + F_d - F_r < -F_break:
            mech['state'] = 2
            iced = False
            return
        return
    else:
        if phi > 0 and i > 0:
            mech['state'] = 1
            return
        if phi < phi_max and i < 0:
            mech['state'] = 2
            return
        if (phi < 0 and i > 0) or (phi > phi_max and i < 0):
            mech['state'] = 5
            return
        return


def spring_damper_forces():
    global F_c, F_d
    F_c = c * (phi * (r / G) - x2)
    F_d = d * (w * (r / G) - dx2)


# compute friction force
def friction_force():
    global F_r
    dx_hys = 0.01
    if ddx2 > 0:
        if dx2 > dx_hys:
            F_r = (F_r_min + x2 * (F_r_max - F_r_min) / x_max)
            return
        F_r = -(F_r_min + x2 * (F_r_max - F_r_min) / x_max)
        return
    if ddx2 < 0:
        if dx2 < -dx_hys:
            F_r = -(F_r_min + x2 * (F_r_max - F_r_min) / x_max)
            return
        F_r = (F_r_min + x2 * (F_r_max - F_r_min) / x_max)
        return
    return


# apply gravitational force if going up
def F_m():
    if dx2 > 0:
        return m * g
    return 0

if start_position:
    u_meas = sio.loadmat('u_down')['u']
else:
    u_meas = sio.loadmat('u_up')['u']

dt = 2e-6 # step size solver
# t_end = (len(u_meas)-1)*dt # simulation stop time
t_end = 0.5
length = int(t_end / dt) + 1
vars = 24


matrix = np.zeros((length, vars))

def sim():
    global t, u, di, i_M, i_R, i, dw, w, n, phi, ddx2, dx2, x2, F_d, F_c, F_r, M_A, M_d, M_c, M_b, M_L, mech, matrix
    matrix[0][0] = 0
    matrix[0][1] = 0
    idx = 0
    print("starting sim...")

    while idx < length:

        if t == 0:
            set_phi()

        if w == 0:
            mech['direction'][0] = 0
        elif signbit(w):
            mech['direction'][0] = -1
        elif not signbit(w):
            mech['direction'][0] = 1

        # u = u_meas[idx]
        u = u_pwm(t, 25e3, 0.8, -12)

        if abs(u) <= 0.1:
            i_M = 0
            i_R = 0
            di = -i / dt
            i = 0
        else:
            di = (u - R * i_M - k * w) / L
            i_M += di * dt
            i_R = i_ripple()
            i = i_M + i_R

        # di = (u - R * i_M - k * w) / L
        # i_M += di * dt
        # i_R = i_ripple()
        # i = i_M + i_R

        if mech['state'] == 1:
            M_A = k * eta_G * i_M
            M_d = F_d * (r / G)
            M_c = F_c * (r / G)
            M_b = b * w
            M_L = M_d + M_c + M_b
            dw = (M_A - M_b - M_d - M_c) / (J * eta_G)
            w += dw * dt
            n = w / (2 * pi)
            phi += w * dt
            spring_damper_forces()
            friction_force()
            ddx2 = (F_c + F_d - F_r - F_m()) / m
            dx2 += ddx2 * dt
            x2 += dx2 * dt

        if mech['state'] == 2:
            M_A = k * eta_G * i_M
            M_d = F_d * (r / G)
            M_c = F_c * (r / G)
            M_b = b * w
            M_L = M_d + M_c + M_b
            dw = (M_A - M_b - M_d - M_c) / (J * eta_G)
            w += dw * dt
            n = w / (2 * pi)
            phi += w * dt
            spring_damper_forces()
            friction_force()
            ddx2 = (F_c + F_d - F_r - F_m()) / m
            dx2 += ddx2 * dt
            x2 += dx2 * dt

        if mech['state'] == 3:
            w = 0
            dw = 0
            M_A = k * eta_G * i_M
            M_d = F_d * (r / G)
            M_c = F_c * (r / G)
            M_b = b * w
            M_L = M_d + M_c + M_b
            n = w / (2 * pi)
            phi += w * dt
            spring_damper_forces()
            # friction_force()
            ddx2 = (F_c + F_d - F_r - F_m()) / m
            dx2 += ddx2 * dt
            x2 += dx2 * dt
            if x2 > x_max:
                x2 = x_max
            if x2 < 0:
                x2 = 0

        if mech['state'] == 4:
            M_A = k * eta_G * i_M
            M_d = F_d * (r / G)
            M_c = F_c * (r / G)
            M_b = b * w
            M_L = M_d + M_c + M_b
            dw = (M_A - M_b - M_d - M_c) / (J * eta_G)
            w += dw * dt
            n = w / (2 * pi)
            phi += w * dt
            spring_damper_forces()
            friction_force()
            ddx2 = 0
            dx2 = 0
            if x2 < 0:
                x2 = 0
            if x2 > x_max:
                x2 = x_max
            if blocked and x2 > x_max - h_B:
                x2 = x_max - h_B

        if mech['state'] == 5:
            M_A = k * eta_G * i_M
            M_d = F_d * (r / G)
            M_c = F_c * (r / G)
            M_b = b * w
            M_L = M_d + M_c + M_b
            dw = (M_A - M_b - M_d - M_c) / (J * eta_G)
            w += dw * dt
            n = w / (2 * pi)
            phi += w * dt
            F_c = 0
            F_d = 0
            F_r = 0
            ddx2 = 0
            dx2 = 0

        if mech['state'] == 6:
            M_A = k * eta_G * i_M
            M_d = F_d * (r / G)
            M_c = F_c * (r / G)
            M_b = b * w
            M_L = M_d + M_c + M_b
            dw = (M_A - M_b - M_d - M_c) / (J * eta_G)
            w += dw * dt
            n = w / (2 * pi)
            phi += w * dt
            spring_damper_forces()
            friction_force()
            ddx2 = 0
            dx2 = 0

        t += dt

        if w == 0:
            mech['direction'][1] = 0
        elif signbit(w):
            mech['direction'][1] = -1
        elif not signbit(w):
            mech['direction'][1] = 1

        mech_state()

        matrix[idx][0] = t
        matrix[idx][1] = u
        matrix[idx][2] = di
        matrix[idx][3] = i_M
        matrix[idx][4] = i_R
        matrix[idx][5] = i
        matrix[idx][6] = dw
        matrix[idx][7] = w
        matrix[idx][8] = phi
        matrix[idx][9] = ddx2
        matrix[idx][10] = dx2
        matrix[idx][11] = x2
        matrix[idx][12] = n
        matrix[idx][13] = M_A
        matrix[idx][14] = M_c
        matrix[idx][15] = M_d
        matrix[idx][16] = M_b
        matrix[idx][17] = M_L
        matrix[idx][18] = mech['state']
        matrix[idx][19] = mech['direction'][0]
        matrix[idx][20] = mech['direction'][1]
        matrix[idx][21] = F_c
        matrix[idx][22] = F_d
        matrix[idx][23] = F_r

        progress(idx+1, length)

        idx += 1




if __name__ == "__main__":
    sim()

    sio.savemat('sim.mat', {'sim': matrix})

    f, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(5, sharex=True)
    ax1.plot(matrix[:, 0], matrix[:, 1])
    ax1.grid(True)
    ax2.plot(matrix[:, 0], matrix[:, 5])
    ax2.grid(True)
    ax3.plot(matrix[:, 0], matrix[:, 7])
    ax3.grid(True)
    ax4.plot(matrix[:, 0], matrix[:, 11])
    ax4.grid(True)
    ax5.plot(matrix[:, 0], matrix[:, 18])
    ax5.grid(True)

    f, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(5, sharex=True)
    ax1.plot(matrix[:, 0], matrix[:, 13])
    ax1.grid(True)
    ax2.plot(matrix[:, 0], matrix[:, 14])
    ax2.grid(True)
    ax3.plot(matrix[:, 0], matrix[:, 15])
    ax3.grid(True)
    ax4.plot(matrix[:, 0], matrix[:, 16])
    ax4.grid(True)
    ax5.plot(matrix[:, 0], matrix[:, 17])
    ax5.grid(True)

    plt.show()


    # matrix.astype('float64').tofile('test')


# bool write_data_bin(int len) {
#      clock_t t2 = std::clock()

#      FILE * f
#      f = fopen("C:\\temp\\data.bin", "wb")
#      for (int i = 0 i < len i++) {
#           for (int k = 0 k < vars k++) {
#                fwrite(&matrix[i][k], sizeof(double), 1, f)
#           }
#      }
#      fclose(f)

#      clock_t t3 = std::clock()
#      std::cout << "writing time: " << double(t3-t2) / CLOCKS_PER_SEC << " s\n"

#      return True
# }


# int main(void) {
#      std::cout.precision(6)
#      clock_t t1 = std::clock()
#      sim()
#      clock_t t2 = std::clock()
#      std::cout << "simulation time: " << double(t2-t1) / CLOCKS_PER_SEC << " s\n"

#      system("del C:\\temp\\data.bin")
#      write_data_bin(length)

#      std::cin.ignore()

#      return 0
# }

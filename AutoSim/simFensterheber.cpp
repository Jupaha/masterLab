// #include <bits/c++config>
#undef _GLIBCXX_DEPRECATED
#define _GLIBCXX_DEPRECATED 

#include <ctime>
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
#include <cstdio>
//#include <boost/tuple/tuple.hpp>
//#include <boost/foreach.hpp>

#include <windows.h>

//#include "gnuplot-iostream.h"

// #using <mscorlib.dll>
// #using <disnet.dll>

// using namespace System;

void pause_if_needed() {
#ifdef _WIN32
     // For Windows, prompt for a keystroke before the Gnuplot object goes out of scope so that
     // the gnuplot window doesn't get closed.
     std::cout << "Press enter to exit." << std::endl;
     std::cin.get();
#endif
}

#define XOR(exp1, exp2) (!((exp1) || (exp2)) || ((exp1) && (exp2)))

int sign(double x) {
    if (x == 0) {
        return 0;
    }
    if (x < 0) {
        return -1;
    }
    if (x > 0) {
        return 1;
    }
    return 2;
}

const double pi = 3.14159265358979323846264338328;

// ### PARAMETER ###
// motor parameter
static double U_B = 12.5;  // batterie voltage
static double R = 0.75;    // motor resistance
static double L = 1500e-6; // motor inductance
static double k = 0.021;   // motor constant
static double J = 1.5e-5;  // moment of inertia
static double b = 2e-5;    // viscouss friction motor

// mechanical parameters of window regulator
static double m = 8.26;   // equivalent window mass
static double c = 45e3;    // spring constant of mech. system
static double d = 310;    // damping constant of mech. system
static double eta_G = 0.395; // worm gear efficiency
static double r = 0.025;   // radius of cable drum
static double G = 73;      // worm gear translation ratio

// parameter for variable friction
static double H = G * (pow(k, 2) * eta_G + b * R) / (k * r); // constant for easier use
static double I_0 = 1.2;
static double I_up_min = 6;
static double I_up_max = 7.5;
static double I_down_min = 3;
static double I_down_max = 4.5;
static double I_deltaminmax = I_up_max - I_up_min;
static double I_m = I_down_min - I_0;
static double I_delta_updown = I_up_min - I_down_min;
static double I_n = (I_up_min + I_up_max + I_down_min + I_down_max) / 4;
static double F_r_min = (I_down_min - I_0) * H;
static double F_r_max = I_deltaminmax * H + F_r_min;

// misc.
static double g = 9.81;      // gravity
// static double phi0 = -58.4;  // angle of loose cable
static double x_max = 0.49; // upper stop of window
static double phi_max = x_max * (G/r);
static double phi_loose = 55;
static double x_start = x_max;

// ripple parameter
static double a = 0.595;                   // current influence factor
static double I_RS = 0.0726;               // peak value current ripple at fast speed
static double n_x = 20;                     // corner speed
static double alpha_start = 0.45 * 2 * pi; // start angle
static double z = 10;                      // commutator segments

static bool iced = true;
static double F_break = H * 0.78 * U_B / R;

static bool blocked = true;
static double h_B = 0.1;

// returns a voltage curve for testing the simulation
double U_test(double t) {
    if (t < 0.1) {
        return 0;
    }
    if (t >= 0.1 && t < 2) {
        return -U_B;
    }
    if (t >= 2 && t < 2.5) {
        return 0;
    }
    if (t >= 2.5 && t < 4.5) {
        return -U_B;
    }
    if (t >= 4.5 && t < 5) {
        return 0;
    }
    if (t >= 5 && t < 9) {
        return U_B;
    }
    if (t >= 9 && t < 9.00005) {
        return -U_B;
    }
    if (t >= 9.00005) {
        return -U_B;
    }
    return 0;
}

// ### VARIABLES ###
static double t = 0;
static double u = 0;
static double di = 0;
static double i_M = 0;
static double i_R = 0;
static double i = 0;
static double dw = 0;
static double w = 0;
static double n = 0;
static double phi = (x_start * (G / r));
static double ddx2 = 0;
static double dx2 = 0;
static double x2 = x_start;
static double F_d = 0;
static double F_c = 0;
static double F_r = 0;
static double M_A = 0;
static double M_d = 0;
static double M_c = 0;
static double M_b = 0;
static double M_L = 0;


void set_phi(void) {
    if (x2 == x_max) {
        phi = (x_start * (G / r)) + phi_loose;
        return;
    }
    if (x2 == 0) {
        phi = -phi_loose;
        return;
    }
    phi = (x_start * (G / r));
    return;
}


double i_idle () {
    return u / (R + (eta_G*pow(k, 2))/b);
}

// computes the overlayed current ripple from the commutation
double i_ripple() {
    return (((a * std::abs(i_M) + 1) * I_RS * n) / sqrt(n * n + n_x * n_x)) *
        sin((phi + alpha_start) * z);
}

struct mech_system {
    int direction[2] = {0, 0};
    int state = 5;
    double phi_switch = (x_max * (G / r));
    bool block_window = false;
    bool block_motor = false;
    bool loose = false;
    bool self_lock = false;
} mech;


void mech_state (void) {
    switch(mech.state) {
    case 1:
        if (std::abs(M_A) < std::abs(M_L) && std::signbit(M_L) != std::signbit(w)) {
            mech.state = 3;
            break;
        }
        if (x2 > x_max) {
            mech.state = 4;
            break;
        }
        if (blocked && x2 > x_max - h_B) {
            mech.state = 4;
            break;
        }
        break;
    case 2:
        if (x2 < 0) {
            mech.state = 4;
            break;
        }
        if (std::abs(M_A) < std::abs(M_L) && std::signbit(M_A) != std::signbit(w)) {
            mech.state = 3;
            break;
        }
        break;
    case 3:
         if (std::abs(M_A) > std::abs(M_L) && std::signbit(M_A) == std::signbit(M_L)
             && i > 0) {
             mech.state = 1;
             break;
         }
         if (std::signbit(M_A) != std::signbit(M_L) && M_A > 0) {
             mech.state = 1;
             break;
         }
         if (std::abs(M_A) > std::abs(M_L) && std::signbit(M_A) == std::signbit(M_L)
             && i < 0) {
             mech.state = 2;
             break;
         }
         if (std::abs(M_A) > std::abs(M_L) && std::signbit(M_A) != std::signbit(M_L)
             && i < 0) {
             mech.state = 2;
             break;
         }
         break;
     case 4:
         if (std::abs(u) < 1 || (std::abs(i_M) < 0.1))  {
             if (x2 <= 0) {
                 phi = -phi_loose;
                 mech.phi_switch = 0;
                 mech.state = 5;
                 break;
             }
             else if (x2 >= x_max) {
                 phi = (x_max * (G / r)) + phi_loose;
                 mech.phi_switch = (x_max * (G / r));
                 mech.state = 5;
                 break;
             }
             else if (blocked && x2 >= x_max - h_B) {
                 phi = ((x_max-h_B) * (G / r)) + phi_loose;
                 mech.phi_switch = ((x_max-h_B) * (G / r));
                 mech.state = 5;
                 break;
             }
         }
         else if (std::abs(M_A) < std::abs(M_L) && std::signbit(M_A) == std::signbit(M_L)
             && M_L != 0 && u == 0) {
             mech.state = 3;
             break;
         }

         else if (std::signbit(M_A) != std::signbit(M_L) && M_A > 0) {
             mech.state = 1;
             break;
         }

         else if (x2 > 0 && x2 < x_max && dx2 > 0)
         {
             mech.state = 1;
             break;
         }
         else if (x2 > 0 && x2 < x_max && dx2 < 0)
         {
             mech.state = 2;
             break;
         }
         break;
     case 5:
         if (w > 0 && phi > mech.phi_switch) {
             mech.state = 1;
             break;
         }
         if (w < 0 && phi < mech.phi_switch) {
             if (iced) {
                 mech.state = 6;
                 break;
             }
             mech.state = 2;
             break;
         }
         break;
     case 6:
         if (F_c + F_d - F_r < -F_break) {
             mech.state = 2;
             iced = false;
             break;
         }
         break;
     default:
         if (phi > 0 && i > 0) {
             mech.state = 1;
             break;
         }
         if (phi < phi_max && i < 0) {
             mech.state = 2;
             break;
         }
         if ((phi < 0 && i > 0) || (phi > phi_max && i < 0)) {
             mech.state = 5;
             break;
         }
         break;
     }
}

void spring_damper_forces(void) {
    F_c = c * (phi * (r / G) - x2);
    F_d = d * (w * (r / G) - dx2);
    return;
}

// compute friction force
void friction_force(void) {
    const double dx_hys = 0.01;
    if (ddx2 > 0) {
        if (dx2 > dx_hys) {
            F_r = (F_r_min + x2 * (F_r_max - F_r_min) / x_max);
            return;
        }
        F_r = -(F_r_min + x2 * (F_r_max - F_r_min) / x_max);
        return;
    }
    if (ddx2 < 0) {
        if (dx2 < -dx_hys) {
            F_r = -(F_r_min + x2 * (F_r_max - F_r_min) / x_max);
            return;
        }
        F_r = (F_r_min + x2 * (F_r_max - F_r_min) / x_max);
        return;
    }
    return;
}

// apply gravitational force if going up
double F_m(void) {
    if (dx2 > 0) {
        return m * g;
    }
    return 0;
}



const double dt = 10e-6; // step size solver
const double t_end = 10; // simulation stop time
const int length = int(t_end / dt) + 1;
const int vars = 24;

double matrix[length][vars];


void sim(void) {
    matrix[0][0] = 0;
    matrix[0][1] = 0;
    int idx = 0;
    std::cout << "starting sim...";

    while (idx < length) {

        if (t == 0) {
            set_phi();
        }

        u = U_test(t);
        if (w == 0) {
            mech.direction[0] = 0;
        } else if (std::signbit(w)) {
            mech.direction[0] = -1;
        } else if (!std::signbit(w)) {
            mech.direction[0] = 1;
        }

        if (u == 0) {
            i_M = 0;
            i_R = 0;
            di = -i / dt;
            i = 0;
        }
        else {
            di = (u - R * i_M - k * w) / L;
            i_M += di * dt;
            i_R = i_ripple();
            i = i_M + i_R;
        }

        if (mech.state == 1) {
            M_A = k * eta_G * i_M;
            M_d = F_d * (r / G);
            M_c = F_c * (r / G);
            M_b = b * w;
            M_L = M_d + M_c + M_b;
            dw = (M_A - M_b - M_d - M_c) / (J * eta_G);
            w += dw * dt;
            n = w / (2 * pi);
            phi += w * dt;
            spring_damper_forces();
            friction_force();
            ddx2 = (F_c + F_d - F_r - F_m()) / m;
            dx2 += ddx2 * dt;
            x2 += dx2 * dt;
        }

        if (mech.state == 2) {
            M_A = k * eta_G * i_M;
            M_d = F_d * (r / G);
            M_c = F_c * (r / G);
            M_b = b * w;
            M_L = M_d + M_c + M_b;
            dw = (M_A - M_b - M_d - M_c) / (J * eta_G);
            w += dw * dt;
            n = w / (2 * pi);
            phi += w * dt;
            spring_damper_forces();
            friction_force();
            ddx2 = (F_c + F_d - F_r - F_m()) / m;
            dx2 += ddx2 * dt;
            x2 += dx2 * dt;
        }

        if (mech.state == 3) {
            w = 0;
            dw = 0;
            M_A = k * eta_G * i_M;
            M_d = F_d * (r / G);
            M_c = F_c * (r / G);
            M_b = b * w;
            M_L = M_d + M_c + M_b;
            n = w / (2 * pi);
            phi += w * dt;
            spring_damper_forces();
            // friction_force();
            ddx2 = (F_c + F_d - F_r - F_m()) / m;
            dx2 += ddx2 * dt;
            x2 += dx2 * dt;
            if (x2 > x_max) {
                x2 = x_max;
            }
            if (x2 < 0) {
                x2 = 0;
            }
        }

        if (mech.state == 4) {
            M_A = k * eta_G * i_M;
            M_d = F_d * (r / G);
            M_c = F_c * (r / G);
            M_b = b * w;
            M_L = M_d + M_c + M_b;
            dw = (M_A - M_b - M_d - M_c) / (J * eta_G);
            w += dw * dt;
            n = w / (2 * pi);
            phi += w * dt;
            spring_damper_forces();
            friction_force();
            ddx2 = 0;
            dx2 = 0;
            if (x2 < 0) {
                x2 = 0;
            }
            if (x2 > x_max) {
                x2 = x_max;
            }
            if (blocked && x2 > x_max - h_B) {
                x2 = x_max - h_B;
            }
        }

        if (mech.state == 5) {
            M_A = k * eta_G * i_M;
            M_d = F_d * (r / G);
            M_c = F_c * (r / G);
            M_b = b * w;
            M_L = M_d + M_c + M_b;
            dw = (M_A - M_b - M_d - M_c) / (J * eta_G);
            w += dw * dt;
            n = w / (2 * pi);
            phi += w * dt;
            F_c = 0;
            F_d = 0;
            F_r = 0;
            ddx2 = 0;
            dx2 = 0;
        }

        if (mech.state == 6) {
            M_A = k * eta_G * i_M;
            M_d = F_d * (r / G);
            M_c = F_c * (r / G);
            M_b = b * w;
            M_L = M_d + M_c + M_b;
            dw = (M_A - M_b - M_d - M_c) / (J * eta_G);
            w += dw * dt;
            n = w / (2 * pi);
            phi += w * dt;
            spring_damper_forces();
            friction_force();
            ddx2 = 0;
            dx2 = 0;
        }

        t += dt;

        if (w == 0) {
            mech.direction[1] = 0;
        } else if (std::signbit(w)) {
            mech.direction[1] = -1;
        } else if (!std::signbit(w)) {
            mech.direction[1] = 1;
        }

        mech_state();

        matrix[idx][0] = t;
        matrix[idx][1] = u;
        matrix[idx][2] = di;
        matrix[idx][3] = i_M;
        matrix[idx][4] = i_R;
        matrix[idx][5] = i;
        matrix[idx][6] = dw;
        matrix[idx][7] = w;
        matrix[idx][8] = phi;
        matrix[idx][9] = ddx2;
        matrix[idx][10] = dx2;
        matrix[idx][11] = x2;
        matrix[idx][12] = n;
        matrix[idx][13] = M_A;
        matrix[idx][14] = M_c;
        matrix[idx][15] = M_d;
        matrix[idx][16] = M_b;
        matrix[idx][17] = M_L;
        matrix[idx][18] = mech.state;
        matrix[idx][19] = mech.direction[0];
        matrix[idx][20] = mech.direction[1];
        matrix[idx][21] = F_c;
        matrix[idx][22] = F_d;
        matrix[idx][23] = F_r;

        idx += 1;
    }

    std::cout << "time steps: " << length << "\n";
}


bool write_data_bin(int len) {
     clock_t t2 = std::clock();

     FILE * f;
     f = fopen("C:\\temp\\data.bin", "wb");
     for (int i = 0; i < len; i++) {
          for (int k = 0; k < vars; k++) {
               fwrite(&matrix[i][k], sizeof(double), 1, f);
          }
     }
     fclose(f);

     clock_t t3 = std::clock();
     std::cout << "writing time: " << double(t3-t2) / CLOCKS_PER_SEC << " s\n";

     return true;
}


int main(void) {
     std::cout.precision(6);
     clock_t t1 = std::clock();
     sim();
     clock_t t2 = std::clock();
     std::cout << "simulation time: " << double(t2-t1) / CLOCKS_PER_SEC << " s\n";

     system("del C:\\temp\\data.bin");
     write_data_bin(length);

     std::cin.ignore();

     return 0;
}

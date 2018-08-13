#include "simPMDCM_Windowlever.h"


#define XOR(exp1, exp2) (!((exp1) || (exp2)) || ((exp1) && (exp2)))

static double pi = 3.14159265358979323846264338328;

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
static double phi0 = -58.4;  // angle of loose cable
static double x_max = 0.49; // upper stop of window
static double phi_max = x_max * (G/r);
static double phi_loose = 55;

// ripple parameter
static double a = 0.595;                   // current influence factor
static double I_RS = 0.0726;               // peak value current ripple at fast speed
static double n_x = 20;                     // corner speed
static double alpha_start = 0.45 * 2 * pi; // start angle
static double z = 10;                      // commutator segments

static int get_params()
{
    std::vector<double> params;
    try
    {
        params = import_params("params.json");
    }
    catch (int e)
    {
        std::cout << "error importing parameters from json file\n";
        return 1;
    }

    if (params.size() != 21)
    {
        std::cout << "parsed parameters don't match!\n";
        std::cout << "using default parameters...\n";
        return 1;
    }

    R = params[0];
    L = params[1];
    k = params[2];
    J = params[3];
    b = params[4];
    eta_G = params[5];
    G = params[6];
    r = params[7];
    m = params[8];
    c = params[9];
    d = params[10];
    I_0 = params[11];
    I_up_max = params[12];
    I_up_max = params[13];
    I_down_min = params[14];
    I_down_max = params[15];
    x_max = params[16];
    a = params[17];
    I_RS = params[18];
    n_x = params[19];
    z = params[20];

    return 0;
}

// returns a voltage curve for testing the simulation
double U_test(double t) {
     if (t < 2) {
          return U_B;
     }
     if (t >= 2 && t < 2.5) {
          return 0;
     }
     if (t >= 2.5 && t < 5) {
          return U_B;
     }
     if (t >= 5 && t < 5.5) {
          return 0;
     }
     if (t >= 5.5 && t < 9.5) {
          return -U_B;
     }
     if (t >= 9.5 && t < 10) {
          return 0;
     }
     if (t >= 10) {
          return U_B;
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
// double phi = (x_max * (G / r)) + 50
static double phi = -phi_loose;
static double ddx2 = 0;
static double dx2 = 0;
// double x2 = x_max;
static double x2 = 0;
static double F_d = 0;
static double F_c = 0;
static double F_r = 0;
static double M_A = 0;
static double M_d = 0;
static double M_c = 0;
static double M_b = 0;
static double M_L = 0;


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
     double phi_switch = 0;
     bool block_window = false;
     bool block_motor = false;
     bool loose = false;
     bool self_lock = false;
} mech;


// void spring_damper_forces(double phi, double x2, double w, double dx2) {
//      if (XOR(phi * (r / G) < x_max, phi * (r / G) > 0)) {
//           F_c = c * (phi * (r / G) - x2);
//           F_d = d * (w * (r / G) - dx2);
//           return;
//      }
//      F_c = 0;
//      F_d = 0;
//      return;
// }

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
          if (mech.direction[0] != mech.direction[1]) {
               mech.state = 5;
               break;
          }
          break;
     case 2:
          if (std::abs(M_A) < std::abs(M_L) && std::signbit(M_A) != std::signbit(w)) {
               mech.state = 3;
               break;
          }
          if (x2 < 0) {
               mech.state = 4;
               break;
          }
          if ((mech.direction[0] < 0) != (mech.direction[1] < 0)) {
               mech.state = 5;
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
          if (std::abs(M_A) > std::abs(M_L) && std::signbit(M_A) != std::signbit(M_L)
              && i < 0) {
               mech.state = 2;
               break;
          }
          break;
     case 4:
          if ((mech.direction[0] < 0) != (mech.direction[1] < 0)) {
               mech.state = 5;
               break;
          }
          if (std::abs(M_A) < std::abs(M_L) && std::signbit(M_A) == std::signbit(M_L)
              && M_L != 0 && u == 0) {
               mech.state = 3;
               break;
          }

          if (std::signbit(M_A) != std::signbit(M_L) && M_A > 0) {
               mech.state = 1;
               break;
          }

          if (x2 > 0 && x2 < x_max && dx2 > 0) {
               mech.state = 1;
               break;
          }
          if (x2 > 0 && x2 < x_max && dx2 < 0) {
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
               mech.state = 2;
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

void spring_damper_forces(double phi, double x2, double w, double) {
     F_c = c * (phi * (r / G) - x2);
     F_d = d * (w * (r / G) - dx2);
     return;
}

// compute friction force
void friction_force(double x2, double dx2) {
     if (dx2 >= -0.01 && dx2 <= 0.01) {
          F_r = 0;
          return;
     }
     if (std::signbit(dx2) == true) {
          F_r = -(F_r_min + x2 * (F_r_max - F_r_min) / x_max);
          return;
     }
     F_r = (F_r_min + x2 * (F_r_max - F_r_min) / x_max);
     return;
}

// apply gravitational force if going up
double F_m(double dx2) {
     if (dx2 > 0) {
          return m * g;
     }
     return 0;
}

std::vector<double> load_window_lever()
{
    std::vector<double> ret;
    get_params();
    ret.push_back(R);
    ret.push_back(L);
    return ret;
}

double solve_window_lever(double Us, double dt) 
{
	//if (t > 3)
	//{
	//	u = 12;
	//}
	//else
	//{
	//	u = 0;
	//}
	u = Us;
	if (std::signbit(w) != 0) {
		mech.direction[0] = std::signbit(w);
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
		spring_damper_forces(phi, x2, w, dx2);
		friction_force(x2, dx2);
		ddx2 = (F_c + F_d - F_r - F_m(dx2)) / m;
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
		spring_damper_forces(phi, x2, w, dx2);
		friction_force(x2, dx2);
		ddx2 = (F_c + F_d - F_r - F_m(dx2)) / m;
		dx2 += ddx2 * dt;
		x2 += dx2 * dt;
	}

	if (mech.state == 3) {
		w = 0;
		//dw = -w / dt;
		dw = 0;
		M_A = k * eta_G * i_M;
		M_d = F_d * (r / G);
		M_c = F_c * (r / G);
		M_b = b * w;
		M_L = M_d + M_c + M_b;
		n = w / (2 * pi);
		phi += w * dt;
		spring_damper_forces(phi, x2, w, dx2);
		friction_force(x2, dx2);
		ddx2 = (F_c + F_d - F_r - F_m(dx2)) / m;
		dx2 += ddx2 * dt;
		x2 += dx2 * dt;
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
		spring_damper_forces(phi, x2, w, dx2);
		friction_force(x2, dx2);
		//ddx2 = (F_c + F_d - F_r - F_m(dx2)) / m;
		ddx2 = 0;
		dx2 = 0;
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
		// friction_force(x2, dx2);
		F_r = 0;
		ddx2 = 0;
		dx2 = 0;
	}

	t += dt;

	if (std::signbit(w) != 0) {
		mech.direction[1] = std::signbit(w);
	}

	mech_state();

	return i;
}
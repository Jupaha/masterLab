#include "simPMDCM_StaticLoad.h"

// Motorparameter
static double R = 0.63; // Ohm
static double L = 950E-6; 
static double J = 1e-5; // kgm^2
static double b = 3.5E-5; // Nms
static double eta_G = 0.7; // -
static double G = 74; // -
static double k = 0.021;
static double k_M = k; //0.03*0.5; // Nm/A
static double k_w = k; //0.025*0.6; // Vs/rad

// Rippleparameter
static double a = 0.6;
static double I_RS = 0.2;
static double n_x = 27;
static double z = 10; // -
static double phi_start = 0;

// Load Torques for different spinning directions
static double M_up = 0.17;
static double M_down = -0.12;

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

    if (params.size() != 14)
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
    a = params[7];
    I_RS = params[8];
    n_x = params[9];
    z = params[10];
    phi_start = params[11];
    M_up = params[12];
    M_down = params[13];

    return 0;
}

static double t = 0;
static double u = 0;
static double di = 0;
static double i_M = 0;
static double i_R = 0;
static double i = 0;
static double dw = 0;
static double w = 0;
static double phi = phi_start;
static double M = 0;


double I_ripple(double n, double i_M, double phi)
{
	return (((a*i_M+1)*I_RS*n)/sqrt(n*n+n_x*n_x))*sin((phi+phi_start)*z);
}

double U_test_window(void)
{
	if (t < 0)
	{
		return 0;
	}
	if (t >= 0 && t < 4.5)
	{
		return -12;
	}
	if (t >= 4.5 && t < 5)
	{
		return 0;
	}
	if (t >= 5 && t < 8)
	{
		return 12;
	}
	if (t >= 8)
	{
		return -12;
	}
	return 0;
}

std::vector<double> load_PMDCM_static_load()
{
    std::vector<double> ret;
    get_params();
    ret.push_back(R);
    ret.push_back(L);
    return ret;
}

double solve_PMDCM_StaticLoad(float Us, float dt)
{
	u = Us;
	t = t + dt;
	di = (u - k_w * w - R * i) / L;
	i_M = i_M + di * dt;
	i_R = I_ripple(w/(2*M_PI), i_M, phi);
	i = i_M + i_R;
    if (w >= 0)
    {
        M = M_up;
    }
    if (w < 0)
    {
        M = M_down;
    }
	dw = (i_M * k_M - w * b - M) / J;
    w = w + dw*dt;
	phi = phi + w * dt;

	return i;
}
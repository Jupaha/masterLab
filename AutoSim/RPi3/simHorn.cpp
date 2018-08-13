#include "simHorn.h"

#define n 20000

// 15 parameters
// parameters have to be ordered like following in the json params file!
static double R = 1.8; // ohm - coil resistance
static double L = 800e-6; // H - coil inductance
static double R_C = 12; // ohm - series resistance of the capacitance
static double C = 2e-6; // F - parallel capacitance to switch

// parameters of transistor pulse generator
static double f = 429; // frequency
static double PW = 0.63; // pulse width

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

    if (params.size() != 6)
    {
        std::cout << "parsed parameters don't match!\n";
        std::cout << "using default parameters...\n";
        return 1;
    }

    R = params[0];
    L = params[1];
    R_C = params[2];
    C = params[3];
    f = params[4];
    PW = params[5];

    return 0;
}

static double t = 0;
static double i = 0;
static double di = 0;
static double u_C = 0;
static double u_D = 0;
static int state = 1;

bool Son()
{
	double T = 1/f;
	double TM = fmod(t, T);
	if ((TM / T) < PW)
	{
		return true;
	}
	return false;
}


// the state discribes the state of the transistor (switch) and
// the diode of the equivalent circuit
// 1: switch closed, diode blocking
// 1 -> 2: Ton false
// 2: switch open, diode blocking
// 2 -> 3: U_D < 0
// 2 -> 1: Ton true
// 3: switch open, diode conducting
// 3 -> 2: i >= 0
// 3 -> 1: Ton true
int getState()
{
    if (state == 1)
    {
        if (Son() == false)
        {
            return 2;
        }
        return 1;
    }
    if (state == 2)
    {
        if (u_D < 0)
        {
            return 3;
        }
        if (Son() == true)
        {
            return 1;
        }
        return 2;
    }
    if (state == 3)
    {
        if (i >= 0)
        {
            return 2;
        }
        if (Son() == true)
        {
            return 1;
        }
        return 3;
    }
    return 1;
}


std::vector<double> load_horn()
{
    std::vector<double> ret;
    get_params();
    ret.push_back(R);
    ret.push_back(L);
    return ret;
}

double solve_horn(double u, double dt)
{
	t += dt;
    if (state == 1 || state == 3)
    {
        i = (L*i + dt*u) / (L  + dt*R);
        u_C = 0;
        u_D = 0;
    }
    if (state == 2)
    {
        i = (L*i + dt*(u-u_C)) / (L + dt*(R+R_C));
        u_C += dt*(i/C);
        u_D = u_C + R_C*i;
    }
    state = getState();

	return i;
}
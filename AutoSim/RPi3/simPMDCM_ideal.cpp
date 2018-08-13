#include "simPMDCM_ideal.h"


// Motorparameter
static double J = 5e-4; // kgm^2
static double b = 1e-3; // Nms
static double k = 0.05;
static double R = 1;
static double L = 1e-3;

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

    if (params.size() != 5)
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

    return 0;
}

static double t = 0;
static double U = 0;
static double dI = 0;
static double I = 0;
static double dw = 0;
static double w = 0;

std::vector<double> load_PMDCM_ideal()
{
    std::vector<double> ret;
    get_params();
    ret.push_back(R);
    ret.push_back(L);
    return ret;
}

double solve_PMDCM_Ideal(float Us, float dt)
{
	U = Us;
	t = t + dt;
	dI = (U - k * w - R * I) / L;
	I = I + dI * dt;
	dw = (I * k - w * b) / J;
	w = w + dt * dw;
	return I;
}
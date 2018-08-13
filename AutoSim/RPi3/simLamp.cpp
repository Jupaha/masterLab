#include "simLamp.h"

// Parameter
static double R_c = 0.37; // Ohm
static double L = 15e-6; // H
static double T_U = 295; // K
static double C_W = 7e-3; // J/K
static double b = 1.5e-12; // W/K^4
static double k = 0.9; // -

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

    R_c = params[0];
    L = params[1];
    T_U = params[2];
    C_W = params[3];
    b = params[4];
    k = params[5];

    return 0;
}

// variables
static double t = 0;
static double di = 0;
static double i = 0;
static double R_W = R_c;
static double P_e = 0;
static double P_t = 0;
static double dT_W = 0;
static double T_W = T_U;

std::vector<double> load_lamp()
{
    std::vector<double> ret;
    get_params();
    ret.push_back(5*R_c);
    ret.push_back(L);
    return ret;
}

double solve_lamp(double u, double dt)
{
    P_e = u*u/R_W;
    P_t = b*std::pow(T_W-T_U, 4);
    dT_W = (P_e-P_t)/C_W;
    T_W = T_W + dT_W*dt;
    R_W = R_c*std::pow(T_W/T_U, k);
    di = (u-i*R_W)/L;
    i = i + di*dt;
    return i;
}
#include "simBLDCM.h"

// ### PARAMETER ###
static double U_An = 12.0; // V - nominal voltage
static double I_An = 1; // A - nomial current
static double w_n = 100; // rad/s - nominal angular speed
static double p = 2; // pole pairs
static double _R_A = 0.5; // normalized armature resistance
static double T_L = 1; // s - mechanical time constant
static double a = -0.0933e-3;
static double b = 304e-3;
static double T_k = 0.2e-3; // s - commutation time constant
static double x_k = 0.85; // normalized commutation time
static double k_k = 1.4; // commutation factor
static double _k_L = 1; // normalized mechanical load


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

    if (params.size() != 12)
    {
        std::cout << "parsed parameters don't match!\n";
        std::cout << "using default parameters...\n";
        return 1;
    }

    U_An = params[0];
    I_An = params[1];
    w_n = params[2];
    p = params[3];
    _R_A = params[4];
    T_L = params[5];
    a = params[6];
    b = params[7];
    T_k = params[8];
    x_k = params[9];
    k_k = params[10];
    _k_L = params[11];

    return 0;
}

// ### VARIABLES ###
static double t = 0; // s - time
static double _u = 0; // V - normalized voltage
static double i = 0; // A - current
static double _i = 0; // A - normalized current
static double _i_A = 0; // A - normalized steady current
static double _i_R = 0; // A - normalized polynomial current
static double _i_k = 0; // A - normalized commutation current
static double _dw = 0; // rad/s^2 - normalized angular acceleration
static double _w = 0; // rad/s - normalized angular speed
static double _phi = 0; // rad - normalized angule
static double x = 0; // - normalized time between commutations
static double t_k = 0; // s - time immediately prior commutation
static double _i_VK = 0; // A - current immediately prior commutation
static bool state = false; // commutation state, true while commutating

std::vector<double> load_BLDCM()
{
    std::vector<double> ret;
    get_params();
    ret.push_back(U_An*_R_A/I_An);
    ret.push_back(T_k*U_An*_R_A / I_An);
    return ret;
}


double solve_BLDCM(double u, double dt)
{
    if (std::abs(u) < 0.5)
    {
        i = 0;
    }
    else
    {
        _u = u / U_An;

        _dw = (_i - _k_L*_w*std::abs(_w));
        _w = _w + _dw*dt;
        _phi = _phi + _w*dt;

        _i_A = (_u - (1-_R_A)*_w) / _R_A;

        x = 2*p*w_n*_phi/(2*M_PI) - std::floor(2*p*w_n*_phi/(2*M_PI));
        _i_R = (a*_w + b)
            * (-(256/3)*std::pow(x-0.5, 8)+(16/3)*std::pow(x-0.5, 2)-(11/27));

        if (state == false)
        {
            if (x > x_k)
            {
                t_k = t;
                state = true;
            }
        }
        else
        {
            if (x < x_k)
            {
                state = false;
            }
        }

        if (state == false)
        {
            _i_k = 0;
        }
        else
        {
            _i_k = -k_k*_i_VK*std::exp(-(t-t_k)/T_k);
        }

        _i = _i_A + _i_R + _i_k;
        i = _i * I_An;
    }

    t += dt;
    
	return i;
}

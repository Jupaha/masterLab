#include "simRL.h"

static double R = 2; // Ohm
static double L = 1e-3; // H

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

    if (params.size() != 2)
    {
        std::cout << "parsed parameters don't match!\n";
        std::cout << "using default parameters...\n";
        return 1;
    }

    R = params[0];
    L = params[1];

    return 0;
}

static double t = 0;
static double di = 0;
static double i = 0;

std::vector<double> load_RL()
{
    std::vector<double> ret;
    get_params();
    ret.push_back(R);
    ret.push_back(L);
    return ret;
}

double solve_RL(double U, double dt)
{
	t += dt;
	di = (U - R*i) / L;
	i += di*dt;
	return i;
}

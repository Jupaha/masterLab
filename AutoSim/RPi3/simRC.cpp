#include "simRC.h"

static double R = 1; // Ohm
static double C = 0.0005; // F

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
    C = params[1];

    return 0;
}

static double t = 0;
static double i = 0;
static double idt = 0;

std::vector<double> load_RC()
{
    std::vector<double> ret;
    get_params();
    ret.push_back(R);
    ret.push_back(0);
    return ret;
}

double solve_RC(double u, double dt)
{
	t += dt;
	i = (u - idt / C) / R;
	idt += i*dt;
	return i;
}

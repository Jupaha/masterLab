#include "testFunctions.h"

static const double pi = 3.14159265358979323846;
static double dt_test = 10e-6;
static double t = 0;
static double i = 0;
static double i_DC = 0;
static double step = 0.001;
static double max = 32.767;
static double amp = 32.767;
static double offset = 0;
static double freq = 0;
static double T = 0;
static double x = 0;
static double sign = 1;
static double PW = 0.5;
static double R = 1;
static double L = 20e-6;

double get_dt()
{
    return dt_test;
}

double load_DC()
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

    if (params.size() != 1)
    {
        std::cout << "parsed parameters don't match!\n";
        std::cout << "using default parameters...\n";
        return 1;
    }

    i_DC = params[0];
    dt_test = 10e-6;

    return 0;
}
double comp_DC()
{
    return i_DC;
}


double load_ramp()
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

    if (params.size() != 4)
    {
        std::cout << "parsed parameters don't match!\n";
        std::cout << "using default parameters...\n";
        return 1;
    }

    step = params[0];
    max = params[1];
    sign = params[2];
    dt_test = params[3];

    return 0;
}
double comp_ramp()
{
    i += step;
    if (i > max)
    {
        i = 0;
    }

    return i*sign;
}


double load_sine()
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

    if (params.size() != 3)
    {
        std::cout << "parsed parameters don't match!\n";
        std::cout << "using default parameters...\n";
        return 1;
    }

    amp = params[0];
    offset = params[1];
    freq = params[2];
    dt_test = 10e-6;

    return 0;
}
double comp_sine()
{
    i = amp*std::sin(2*pi*freq*t) + offset;
    t += dt_test;

    return i;
}


double load_triangle()
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

    if (params.size() != 3)
    {
        std::cout << "parsed parameters don't match!\n";
        std::cout << "using default parameters...\n";
        return 1;
    }

    amp = params[0];
    offset = params[1];
    freq = params[2];
    dt_test = 10e-6;

    T = 1/freq;

    return 0;
}
double comp_triangle()
{
    x = fmod(t, T) / T;
    if (x <= 0.5)
    {
        i = (2*x - 0.5)*amp + offset;
    }
    else
    {
        i = (1.5 - 2*x)*amp + offset;
    }
    t += dt_test;

    return i;
}


double load_pulse()
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

    if (params.size() != 4)
    {
        std::cout << "parsed parameters don't match!\n";
        std::cout << "using default parameters...\n";
        return 1;
    }

    amp = params[0];
    offset = params[1];
    freq = params[2];
    PW = params[3];
    dt_test = 10e-6;

    T = 1/freq;

    return 0;
}
double comp_pulse()
{
    x = fmod(t, T) / T;
    if (x <= PW)
    {
        i = amp + offset;
    }
    else
    {
        i = 0;
    }
    t += dt_test;

    return i;
}

std::vector<double> load_RL_HP()
{
    std::vector<double> ret;

    std::vector<double> params;
    try
    {
        params = import_params("params.json");
    }
    catch (int e)
    {
        std::cout << "error importing parameters from json file\n";
        ret.push_back(999);
        return ret;
    }

    if (params.size() != 3)
    {
        std::cout << "parsed parameters don't match!\n";
        std::cout << "using default parameters...\n";
        ret.push_back(999);
        return ret;
    }

    R = params[0];
    L = params[1];
    offset = params[2];

    ret.push_back(R);
    ret.push_back(L);
    ret.push_back(offset);

    return ret;
}
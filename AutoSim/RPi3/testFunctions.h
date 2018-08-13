#pragma once

#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include "misc.h"
#include "peri.h"

double get_dt();

double load_DC();
double comp_DC();

double load_ramp();
double comp_ramp();

double load_sine();
double comp_sine();

double load_triangle();
double comp_triangle();

double load_pulse();
double comp_pulse();

std::vector<double> load_RL_HP();
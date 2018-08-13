#pragma once

#include <cmath>
#include <iostream>
#include <vector>
#include "misc.h"

/// load parameters and return tau for R/C
std::vector<double> load_horn();

/// solve equations
double solve_horn(double u, double dt);
#pragma once

#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include "misc.h"

std::vector<double> load_BLDCM();

double solve_BLDCM(double u, double dt);
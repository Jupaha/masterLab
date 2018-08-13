#pragma once

#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include "misc.h"

std::vector<double> load_horn();

double solve_RC(double u, double dt);
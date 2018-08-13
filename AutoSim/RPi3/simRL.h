#pragma once

#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include "misc.h"

std::vector<double> load_RL();

double solve_RL(double U, double dt);
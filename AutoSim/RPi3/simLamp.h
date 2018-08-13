#pragma once

#include <ctime>
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
#include <cstdio>
#include "misc.h"

std::vector<double> load_lamp();

double solve_lamp(double u, double dt);
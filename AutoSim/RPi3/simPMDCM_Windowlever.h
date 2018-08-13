#pragma once

#include <ctime>
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
#include <cstdio>
#include "misc.h"

std::vector<double> load_window_lever();

double solve_window_lever(double Us, double dt);
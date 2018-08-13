#pragma once

#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <unistd.h>
#include <chrono>
#include <vector>
#include "peri.h"
#include "misc.h"

std::vector<double> load_PMDCM_ideal();

double solve_PMDCM_Ideal(float Us, float dt);
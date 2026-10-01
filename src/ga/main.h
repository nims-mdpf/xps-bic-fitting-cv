#pragma once

#define _USE_MATH_DEFINES
#include <cstdio>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <random>
#include <algorithm>
#include <numeric>
#include <cmath>

#include "util.h"
#include "parameterStructure.h"
#include "genetic_algorithm.h"

using namespace std;

int readcsv(const string &fname, Spectrum &data);
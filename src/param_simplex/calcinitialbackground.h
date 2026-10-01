#pragma once
#include <cstdio>
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include "settings.h"
#include "minimumpointsparameter.h"
#include "optimizeparameters.h" // To use make_background_Linear()

using namespace std;

int calcInitialBackground(const int ndata, const vector<double> &bindingEnergy,
	const vector<double> &spectrum, const vector<double> &smoothedSpectrum,
	const vector<MinimumPoints*> &minimumPointsList,
	vector<double> &background, vector<double> &diffSpectrum);

int calcInitialBackgroundShirley(const int ndata, const vector<double> &bindingEnergy,
	const vector<double> &spectrum, const vector<MinimumPoints*> &minimumPointsList,
	vector<double> &background,
	vector<double> &diffSpectrum, const double initBGLower);

int calcInitialBackgroundLinear(const int ndata,
	const vector<double> &spectrum, const vector<double> &smoothedSpectrum,
	vector<double> &background, vector<double> &diffSpectrum);

int calcInitialBackgroundLinear2(const int ndata, const vector<double> &bindingEnergy,
	const vector<double> &spectrum, const vector<double> &smoothedSpectrum,
	const vector<MinimumPoints*> &minimumPointsList,
	vector<double> &background, vector<double> &diffSpectrum);

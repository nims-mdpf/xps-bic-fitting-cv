#pragma once
#include <vector>
#include "settings.h"


class PeakParameter
{
public:
	double height;       // Height
	double position;     // Position (eV)
	double sigma;        // sigma (eV)
	double gamma;        // gamma (eV)
	PeakParameter(double _h, double _p, double _s, double _g) :
		height(_h), position(_p), sigma(_s), gamma(_g) {}
	~PeakParameter() {}
};


class FittingParameters
{
public:
	int numpeak;
	double bgHigher;
	double bgLower;
	std::vector<PeakParameter*> peakList;
	//Settings* settings;
	vector<double> x;
	vector<double> s;
	vector<double> f;
	vector<double> qpq;
	vector<double> bg;
	vector<vector<double> > peaks;
	// void make_spectrum_and_peak();
	// void make_background();
	int update();
	// vector<double> signal;
    // vector<vector<double> > peaks;

	// int outputCommonParameters(const string ofname);
	// int makePeaks(vector<double> &x, DataObj *d);
	// int resizer();
	~FittingParameters() {}
};
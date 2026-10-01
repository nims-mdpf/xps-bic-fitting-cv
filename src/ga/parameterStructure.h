#pragma once

#include <cstdio>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include "tblfunc4.h"
#include <gsl/gsl_sf_gamma.h>

using namespace std;

typedef struct _data_obj {
    FuncTbl3 *Vg;
} DataObj;

struct InfoGenetic{
	int MaxOfLoop;
    int generation;
    int individual;
    int split;
    int seed;
	int num_island;
    double p_upper;
    double p_lower;
	double gw_alpha;
	double gw_beta;
	double lw_alpha;
	double lw_beta;
    double prob_cross;
	double prob_exchange;
    double prob_mutation;
	double beta_log_upper;
	double beta_log_lower;
    double beta;
    double lambda;
	double bgh_alpha;
	double bgh_beta;
	double bgl_alpha;
	double bgl_beta;
    // double bg_stddev;
	// double shift_stddev;
    double intensity_alpha;
	double intensity_beta;
	double intensity_max;
	string noiseType;

	DataObj d;
};

class PeakParameter
{
public:
	double area;
	double position;
	double width;
	double lorentzRatio;
	PeakParameter(double _a, double _p, double _w, double _mu) :
		area(_a), position(_p), width(_w), lorentzRatio(_mu) {}
	PeakParameter() {}
	~PeakParameter() {}
};

// Conv Voigt Parameters
class CVP
{
public:
	double area;
	double position;
	double gauss_width;
	double lorentz_width;
	CVP(double _a, double _p, double _w, double _mu) :
		area(_a), position(_p), gauss_width(_w), lorentz_width(_mu) {}
	CVP() {}
	~CVP() {}
};

class Spectrum
{
public:
	int size;
	vector<double> x;
	vector<double> y;
	~Spectrum() {}
};

class FittingParameters
{
public:
	int peaks_size;
	vector<CVP> parameters;
	vector<double> signal;
    vector<vector<double> > peaks;

	int outputCommonParameters(const string ofname);
	int makePeaks(vector<double> &x, DataObj *d);
	int resizer();
	~FittingParameters() {}
};

class Node
{
public:
	vector<int> pair;
	Node(int _a, int _b){
		pair.resize(2);
		pair[0] = _a;
		pair[1] = _b;
	}
	~Node() {}
};

class Series
{
public:
	int index;
	string path;
	FittingParameters fitting_parameters;
	vector <double> background;
	vector <double> fitting;
	Spectrum data;
	double low;
	double high;

	int makeSumSpectrum(DataObj *d);
	int makeLocalBackground();
	int outputFitting(const string kwd, DataObj *d);

	~Series() {}
};

class SolutionStruct
{
public:
	Series series;
	double error;
	double energy;
	double log_prior_density;

	int calc_log_likelihood(InfoGenetic &info);
	int calc_log_prior_density(InfoGenetic &info);
	int calc_energy(InfoGenetic &info);

	~SolutionStruct() {}
};

class InputClass
{
public:
	vector<SolutionStruct> geneticParames;
	InfoGenetic info;
	vector<string> labels;
	double best_energy;
	vector<double> log;
	SolutionStruct bestGeneticParame;
	~InputClass() {}
};
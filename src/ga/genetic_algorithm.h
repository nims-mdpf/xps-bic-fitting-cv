// #pragma once

#define _USE_MATH_DEFINES
#include <cstdio>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <random>
#include <algorithm>
#include <iterator>
#include <numeric>
#include <cmath>
#include <thread>
#include <mutex>
//#include <experimental/filesystem>
//#include <filesystem>
#include <gsl/gsl_sf_gamma.h>
#include <iomanip>

#include "util.h"
#include "parameterStructure.h"
#include "tblfunc4.h"

#define DEBUGTBLFUNC 0

using namespace std;

int SetupDataObj3();
int execDistributedGeneticAlgorithm(int &random_seed, SolutionStruct &solution);
int execGeneticAlgorithm(InputClass& input);
int initializeParameters(vector<SolutionStruct> &geneticParames, InfoGenetic &info);
int SamplingFromPriorDistributions(SolutionStruct &geneticParame, InfoGenetic &info);
int sortPositionParameter(vector<CVP> &parameters);
int calcFitness(vector<SolutionStruct> &geneticParames, InfoGenetic info, vector<double> &energy_list);
vector<SolutionStruct> selectIndividual(vector<SolutionStruct> &geneticParames, InfoGenetic &info);
vector<double> fitness2prob(vector<double> fitness);
vector<int> prob2indexList(const InfoGenetic &info, vector<double> prob);
int calcSelectionPressure(InfoGenetic &info, vector<double> &fitness_list);
int calcPorbList(vector<SolutionStruct> &geneticParames,  InfoGenetic &info,
vector<int> &index_list, vector<vector<int> > &local_index_list);
vector<SolutionStruct> mutateGene(vector<SolutionStruct> &geneticParames, InfoGenetic &info);
vector<SolutionStruct> crossOverGene(vector<SolutionStruct> &geneticParames, const InfoGenetic &info);
vector<SolutionStruct> calcCrossover(SolutionStruct &struct_dad, SolutionStruct &struct_mam, const InfoGenetic &info);
double generateChild(double &dad, double &mam, const InfoGenetic &info);
double random_mt(const int &random_seed, double upper, double lower);
double random_gamma_mt(const int &random_seed, double alpha, double beta);
double random_gauss_mt(const int &random_seed, double mean, double stddev);
int outputEnergylog(vector<InputClass> &inputs, const string &ofname);
vector<int> makeRandomIndexList(const int &size, const int &seed);
int calcMinimumEnergyResult(int &num_island, SolutionStruct &result, vector<InputClass> &inputs);
int saveFiles(const string ss, SolutionStruct &solution, vector<InputClass> &inputs, DataObj d);
int saveAllIsland(vector<InputClass> &inputs);

int get_numThreadLimit();
int get_omp_num_threads();

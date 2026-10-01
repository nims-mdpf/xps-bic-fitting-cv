#pragma once
#define _USE_MATH_DEFINES
#include <cmath>
#include <algorithm>
#include <vector>
#include <random>
#include <iostream>
#include <fstream>
#include <thread>
#include "peakparameter.h"
#include "tblfunc4.h"

typedef struct _data_obj {
    FuncTbl3 *Vg;
} DataObj;
// DataObj SetupDataObj();

extern DataObj d;

using namespace std;

int SetupDataObj3();
// int SetupDataObj2(char* argv0);

double lorentzFunc(const double x, const double c, const double w);
double gaussFunc(const double x, const double c, const double w);
double voigtFunc(const double x, const PeakParameter* pp);
double voigtArea(const PeakParameter* pp);
int make_spectrum_and_peak(const vector<double> &xx, const vector<PeakParameter*> &peakList, vector<double> &sumSpectrum, vector<vector<double> > &peaks);
int make_spectrum(const vector<double> &xx, const vector<PeakParameter*> &peakList, vector<double> &sumSpectrum);
double ran0();
int poissonNoise(double lambda);
double gaussNoise(double sigma, double lambda);
int savitzkyGolayGeneral(const vector<double> &y, const int numSmooth, const int &order, vector<double> &s);
int existFile(const char* path);
double calcRMSE(const vector<double> &v1, const vector<double> &v2, const int errType=1);

/*
* ---------------------------------------------------------------------
* Release memory of a vector variable like vector<T*>
* T is an arbitrary class, e.g. int, double, MyClass, etc.
* ---------------------------------------------------------------------
* // Example
* vector<MyClass*> vec;
* vec.push_back(new MyClass(...));
* ...
* // When the vector is no longer needed
* releaseMemory(vec);
* ---------------------------------------------------------------------
*/
template<typename T>
void releaseMemoryT(vector<T*> &vList)
{
	for (T* v : vList) { delete v; }
	vector<T*>().swap(vList);
}

/*
* ---------------------------------------------------------------------
* Make a copy of a vector variable like vector<T*>
* T is an arbitrary class, e.g. int, double, MyClass, etc.
* ---------------------------------------------------------------------
*/
template<typename T>
void copyVectorT(const vector<T*> &vin, vector<T*> &vout)
{
	releaseMemoryT(vout);
	for (T* v : vin) { vout.push_back(new T(*v)); }
}

/*
* "in" function like python
*/
template<typename Container>
bool exist_in(const Container& c, const typename Container::value_type& v) {
	return (c.end() != std::find(c.begin(), c.end(), v));
}


//
// http://yuta1402.hatenablog.jp/entry/2017/12/15/233615
//
// template<class Iterator, class Function>
// inline Function parallel_for_each(Iterator begin, Iterator end, Function f)
// {
//     if(begin == end) {
//         return std::move(f);
//     }

//     std::size_t num_threads = std::thread::hardware_concurrency();
//     std::size_t step = std::max<std::size_t>(1, std::distance(begin, end) / num_threads);

//     std::vector<std::thread> threads;

//     for(; begin < end - step; begin += step) {
//         threads.emplace_back([=, &f](){
//             std::for_each(begin, begin + step, f);
//         });
//     }

//     threads.emplace_back([=, &f](){
//         std::for_each(begin, end, f);
//     });

//     for(auto&& t : threads) {
//         t.join();
//     }

//     return std::move(f);
// }

// template<class Container, class Function>
// inline Function parallel_for_each(Container&& c, Function f)
// {
//     return parallel_for_each(std::begin(c), std::end(c), f);
// }


//
// https://rextester.com/JRMW42020
//
// template <typename InputIter, typename UnaryFunction>
// void parallel_for_each(InputIter first, InputIter last, UnaryFunction f)
// {
//     auto range = std::distance(first, last);
//     size_t approxNumThreads = std::thread::hardware_concurrency();

//     // number of elements for one thread
//     size_t jobsForThread = range / approxNumThreads;
//     // number of elements for one thread + remainder
//     size_t jobsForMainThread = range % approxNumThreads + jobsForThread;

//     // minus main thread
//     std::vector<std::thread> threads(approxNumThreads - 1);

//     InputIter block_start = first;
//     InputIter block_end = first + jobsForThread;

//     for(size_t i = 0; i < approxNumThreads - 1; ++i)
//     {
//         threads[i] = std::thread(std::for_each<InputIter, UnaryFunction>,
//                                  block_start, block_end, f);

//         block_start = block_end;
//         block_end += jobsForThread;
//     }

//     std::for_each(block_start, block_start + jobsForMainThread, f);

//     for(auto& t : threads)
//     {
//         if(t.joinable())
//         {
//             t.join();
//         }
//     }
// }

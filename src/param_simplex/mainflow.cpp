#include "mainflow.h"
#include <algorithm>
#include <omp.h>

/*
* ---------------------------------------------------------------------
* Main flow of this program
* There are two major loops for the round-robin smoothing.
* 1st loop:
*   - Savitzky-Golay smoothing
*   - Detect minimum points
*   - Remove unnecessary minimum points
*   - Get a boundary of energy range
*   - Make initial background spectrum by the iterative Shirley method
*   - Detect initial peaks
*
* After 1st loop, get a global boundary of energy range
*
* 2nd loop:
*   - Optimize the peaks and background 
* ---------------------------------------------------------------------
* INPUT
*    originalBindingEnergy: Original binding energies (scaled)
*    originalSpectrum:      Original spectrum (scaled)
*    settings:              Settings
* ---------------------------------------------------------------------
* H. Shinotsuka 2019/08/08
* ---------------------------------------------------------------------
*/
void mainFlow(const vector<double> &originalBindingEnergy,
	const vector<double> &originalSpectrum, 
	vector<FittingParameters> &fps, MyTimer &mytimer)
{
	cout << "mainFlow()" << endl;
	cout << endl;
	cout << "#=========================================================#" << endl;
	cout << "#  List of FittingParameters is given.                    #" << endl;
	cout << "#  Each FittingParameters will be optimized.              #" << endl;
	cout << "#=========================================================#" << endl;
	cout << endl;

	// Number of data points
	const int originalNdata = (int)originalBindingEnergy.size();

	int idHigherGlobal, idLowerGlobal;
	// getGlobalEnergyRange(optModelList, idHigherGlobal, idLowerGlobal);
	// cout << "Debug. idHigherGlobal and idLowerGlobal are set to all data";
	idHigherGlobal = 0;
	idLowerGlobal = originalNdata - 1;
	cout << "idHigherGlobal = " << idHigherGlobal << endl;
	cout << "idLowerGlobal  = " << idLowerGlobal << endl;

	// Slice between [idHigherGlobal, idLowerGlobal]
	int ndata = idLowerGlobal - idHigherGlobal + 1;
	// int m_SG = om->mSG;
	vector<double> bindingEnergySlice = cutoffEnergyRange(originalBindingEnergy, idHigherGlobal, idLowerGlobal);
	vector<double> spectrumSlice = cutoffEnergyRange(originalSpectrum, idHigherGlobal, idLowerGlobal);
	// vector<double> backgroundSlice = cutoffEnergyRange(om->background, idHigherGlobal, idLowerGlobal);

	// Set global spect
	spect.setdata(bindingEnergySlice, spectrumSlice);

	// Start optimization using thread
	int numModel = fps.size();
	printf("Number of models = %d\n", numModel);

    // unsigned int NumOfCores=thread::hardware_concurrency();

	// // int numThreadLimit = NumOfCores/2;
	// int numThreadLimit = NumOfCores - 1;
	// numThreadLimit = max(1, numThreadLimit);
	// printf("Total number of cores = %d\n", NumOfCores);
	// printf("Upper limit of number of cores = %d\n", numThreadLimit);

	int numThreadLimit = get_numThreadLimit();
	omp_set_num_threads(numThreadLimit);

	// Parallelization by openMP
	#pragma omp parallel for schedule(static, 1)
		for(int ll=0; ll<numModel; ll++){ 
			optimizeFlow_wrapper(ll, ref(fps[ll]), originalBindingEnergy, originalSpectrum, idHigherGlobal, idLowerGlobal);
		}

	// // Partially parallel loop. Reduce CPU load
	// int ll = 0;
	// do {
	// 	vector<thread> thread_list;
	// 	for(int i=0; i < numThreadLimit; i++){
	// 		printf("dbg. parallel loop for ll = %d, core ID = %d\n", ll, i);
	// 		thread_list.emplace_back(thread( optimizeFlow_wrapper, ll, ref(fps[ll]), originalBindingEnergy, originalSpectrum, idHigherGlobal, idLowerGlobal ));
	// 		ll++;
	// 		if (ll >= numModel){break;}
	// 	}
	// 	for(thread &th : thread_list){ th.join(); }
	// } while(ll < numModel - 1);

	// Parallel loop --- Too high load to CPU, at most 155 threads
	// vector<thread> thread_list(numModel);
	// for(int ll=0; ll<numModel; ll++){ 
	// 	thread_list[ll] = thread( optimizeFlow_wrapper, ll, ref(fps[ll]), originalBindingEnergy, originalSpectrum, idHigherGlobal, idLowerGlobal); 
	// }
	// for(thread &th : thread_list){ th.join(); }

	// Serial loop
	// for(int ll=0; ll<numModel; ll++){ 
	// 	optimizeFlow_wrapper(ll, fps[ll], originalBindingEnergy, originalSpectrum, idHigherGlobal, idLowerGlobal); 
	// }


	// // 非効率．後半にピーク本数の多いものがたくさんいるため，一部のスレッドだけ計算時間が長くなる．
	// // https://www.slideshare.net/yohhoy/boostjp14-parfor-20140301-31781446
	// //並列化 powered by C++11 
	// int N = numModel;
	// // int m = max(thread::hardware_concurrency(), 1u); 
	// int m = NumOfCoresLimit;
	// vector<thread> worker; 
	// for (int i = 0; i < m; ++i) { 
	// 	worker.emplace_back([&](int id) { 
	// 		int r0 = N/m * id + min(N%m, id); 
	// 		int r1 = N/m * (id+1) + min(N%m, id+1); 
	// 		for (int j = r0; j < r1; ++j) { 
	// 			// process(data[j]); 
	// 			// optimizeFlow_wrapper(optModelList[j]);
	// 			optimizeFlow_wrapper(j, fps[j], originalBindingEnergy, originalSpectrum, idHigherGlobal, idLowerGlobal); 
	// 		}
	// 	}, i);
	// } for (auto& t : worker) t.join();	


	mytimer.show("after optimization");
	cout << "mainFlow() finished" << endl;
	return;
}

int optimizeFlow_wrapper(int ll, FittingParameters& fp, const vector<double> &originalBindingEnergy,
	const vector<double> &originalSpectrum, int idHigherGlobal, int idLowerGlobal)
{
	// FittingParameters* fp = &fps[ll];
	cout << "Optimization for ll = " << ll << endl;


	fp.x = spect.x;
	fp.update();

	bool isConverged = false;
	int numTotalOptimizeLoop = 0;

	// Optimize paramters between [idHigherGlobal, idLowerGlobal]
	if (optimizeFlow(ll, fp, spect.x, spect.y, isConverged, numTotalOptimizeLoop)) {
		cout << "Warning at peakSeparationMain()." << endl;
		exit(1);
	}


	//Output all figures for m_SG
	if (1) {
		ostringstream oss;
		oss << setw(4) << setfill('0') << ll;
		output_result("peak_search_ll" + oss.str(), spect.x, spect.y, fp.bg, fp.peakList, originalBindingEnergy, originalSpectrum, idHigherGlobal, idLowerGlobal);
	}
	return 0;
}

int get_numThreadLimit(){
	// Start optimization using thread
    int NumOfCores = (int)thread::hardware_concurrency();
	printf("Total number of cores = %d\n", NumOfCores);
	int omp_num_threads = get_omp_num_threads();
	printf("OMP_NUM_THREADS = %d\n", omp_num_threads);

	int numThreadLimit = min(max(1, omp_num_threads), max(1, NumOfCores - 1));
	printf("Upper limit of number of cores = %d\n", numThreadLimit);
	return numThreadLimit;
}

int get_omp_num_threads() {
	int num_threads = 0;
	const char *tmp = getenv("OMP_NUM_THREADS");
	string env_var(tmp ? tmp : "");
	if (env_var.empty()) {
		cerr << "[ERROR] No such variable found!" << endl;
		//exit(EXIT_FAILURE);
	} else {
		num_threads = stod(env_var);
	}
	return num_threads;
}

/*
* ---------------------------------------------------------------------
* Make a list for round-robin smoothing.
* Each element means a windows size of a Savitzky-Golay smoothing.
* ---------------------------------------------------------------------
*/
vector<int> make_smoothingList()
{
	// Make a list of {5, 5, 5, ..., 7, 7, 7, ..., ..., 15, 15, 15}
	vector<int> smList;
	vector<vector<int>> smoothingWindowList = { { 5,480 },{ 7,480 },{ 9,480 },{ 11,480 },{ 13,480 },{ 15,480 } };
	for (int i = 0; i < smoothingWindowList.size(); i++) {
		smoothingWindowList[i][0], smoothingWindowList[i][1];
		for (int j = 0; j < smoothingWindowList[i][1]; j++) {
			smList.push_back(smoothingWindowList[i][0]);
		}
	}

	cout << "Number of smoothing = " << smList.size() << endl;
	for (auto &sg_window : smList) { cout << sg_window << " "; }
	cout << endl;

	return smList;
}

/*
* ---------------------------------------------------------------------
* Make a pick up list, which contains {0, 1}.
* The size of the list is the same as the number of round-robin smoothing.
* When "0", smoothing is done, but optimization is skipped.
* When "1", both smoothing and optimization is done.
* ---------------------------------------------------------------------
* e.g. {0, 1, 0, 1, ..., 0, 0, 0, 1, 0, 0, 0, 1, ...}
* ---------------------------------------------------------------------
*/
vector<int> make_pickupList(const int n)
{
	vector<int> pickupList;
	for (int i = 1; i <= n; i++) {
		int pickup = 1;
		if (i <= 40) {
			if (i % 2 != 0) { pickup = 0; }
		}
		else if (i <= 80) {
			if (i % 4 != 0) { pickup = 0; }
		}
		else if (i <= 240) {
			if (i % 8 != 0) { pickup = 0; }
		}
		else if (i <= 480) {
			if (i % 16 != 0) { pickup = 0; }
		}
		else if (i <= 960) {
			if (i % 16 != 0) { pickup = 0; }
		}
		else {
			if (i % 32 != 0) { pickup = 0; }
		}
		pickupList.push_back(pickup);
	}

	int numPickupList = (int)pickupList.size();
	cout << "Number of pickupList = " << numPickupList << endl;
	int numDoOptimize = accumulate(pickupList.begin(), pickupList.end(), 0);
	int numSkipOptimize = n - numDoOptimize;
	cout << "Number of DoOptimize = " << numDoOptimize << endl;
	cout << "Number of SkipOptimize = " << numSkipOptimize << endl;

	if (n != numPickupList) {
		cout << "Error: Number of pickupList is inconsitent with number of smoothing." << endl;
		exit(1);
	}

	return pickupList;
}

/*
* ---------------------------------------------------------------------
* Cut off energy range
* INPUT
*    x: vector data
*    idS: ID of a starting point
*    idE: ID of a end point
* RETURN
*    A sliced data
* ---------------------------------------------------------------------
*/
vector<double> cutoffEnergyRange(const vector<double> &x, const int idS, const int idE) {
	return vector<double>(x.begin() + idS, x.begin() + idE + 1);
}

/*
* Get two IDs of global energy range which will be used at the parameter optimization
*/
void getGlobalEnergyRange(const vector<OptimalModel*> &optModelList, int &idHigherGlobal, int &idLowerGlobal)
{
	idHigherGlobal = optModelList[0]->idHigher;
	idLowerGlobal = optModelList[0]->idLower;
	//cout << "--------------------------------------------" << endl;
	for (auto om : optModelList) {
		//cout << om->idHigher << " " << om->idLower << endl;
		idHigherGlobal = min(idHigherGlobal, om->idHigher);
		idLowerGlobal = max(idLowerGlobal, om->idLower);
	}
	//cout << "--------------------------------------------" << endl;
	//cout << idHigherGlobal << " " << idLowerGlobal << endl;
	//cout << "--------------------------------------------" << endl;
}


// Get a minimum and maximum peak width from a peakList_Temp
void summarize_peakList(const vector<PeakParameter*> &peakList_Temp,
	double &min_width, double &max_width)
{
	min_width = 0.0;
	max_width = 0.0;

	if (peakList_Temp.size() == 0) {
		return;
	}

	double width = 0;
	width = peakList_Temp[0]->sigma + peakList_Temp[0]->gamma;
	min_width = width;
	max_width = min_width;

	for (int j = 0; j < peakList_Temp.size(); j++) {
		width = peakList_Temp[j]->sigma + peakList_Temp[j]->gamma;
		min_width = min(min_width, width);
		max_width = max(max_width, width);
	}
}

void getBoundary(vector<MinimumPoints*> &minimumPointsList, int &idHigher, int &idLower,
	double &bgHigher, double &bgLower)
{
	if (minimumPointsList.size() < 2) {
		cout << "Error: minimumPointList does not have two or more elements." << endl;
		exit(1);
	}

	double mostMinimumPoint_height = minimumPointsList.front()->height;
	idHigher = minimumPointsList.front()->id;
	bgHigher = minimumPointsList.front()->height;
	idLower = idHigher;
	for (int i = 1; i < minimumPointsList.size(); i++) {
		if (minimumPointsList[i]->height < mostMinimumPoint_height) {
			mostMinimumPoint_height = minimumPointsList[i]->height;
			idLower = minimumPointsList[i]->id;
			bgLower = minimumPointsList[i]->height;
		}
	}
	//if (idHigher == idLower) {
	//	cout << "Warning: minimumPointList does not have a range." << endl;
	//}
}

int hasNegative(const vector<double>& v)
{
	for (int i = 0; i < v.size(); i++) {
		if (v[i] < 0.0) {
			return 1;
		}
	}
	return 0;
}

void pickupModel(vector<OptimalModel*> &optModelList, const double nrmseThreshold) {
	// Remove optModels whose NRMSE are greater than a threshold
	int cnt = 0;
	const int numModel = (int)optModelList.size();
	for (int i = numModel - 1; i >= 0; i--) {
		if (optModelList[i]->nrmseSmoothedSpectrum > nrmseThreshold) {
			optModelList.erase(optModelList.begin() + i);
			cnt++;
		}
	}
	cout << cnt << " models are removed because of a large NRMSE values" << endl;
	return;
}
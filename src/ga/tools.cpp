#define _USE_MATH_DEFINES
#define SQR(x) ((x)*(x))

#include <cstdio>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <limits>
#include <iomanip>

#include "parameterStructure.h"
#include "util.h"
#include "tblfunc4.h"

// ガウス分布に基づく事前確率の、対数の負符号
double neg_log_prior_Gauss(const double x, const double nu, const double xi)
{
    return 1.0 / (2.0 * xi) * SQR(x - nu) + 0.5 * log(2 * M_PI * xi);
}

// ガンマ分布に基づく事前確率の、対数の負符号
double neg_log_prior_Gamma(const double x, const double eta, const double theta)
{
    return -(eta - 1) * log(x) + x / theta + gsl_sf_lngamma(eta) + eta * log(theta);
}

int Series::makeSumSpectrum(DataObj *d){

    fitting_parameters.makePeaks(data.x, d);
    fitting.resize(data.size);
    for(auto& f : fitting){ f = 0.0; }

    int num_peaks = fitting_parameters.peaks_size;
    for (int n = 0, size = data.size; n < size; n++){
        for (int k = 0; k < num_peaks; k++){
            fitting[n] += ( fitting_parameters.peaks[k][n]);
        }
    }

    return 0;
}

int Series::makeLocalBackground(){
    int ndata = fitting.size();
    double all_area = 0;
    vector<double> areaP(ndata);
    background.resize(ndata);
    areaP[0] = 0.0;
    for (int i = 1; i < ndata; i++) {
        areaP[i] = areaP[i - 1] + (fitting[i - 1] + fitting[i]) * 0.5;
    }
    all_area = areaP[ndata - 1];

    for (int i = 0; i < ndata; i++) {
        background[i] = (high - low) * (all_area - areaP[i]) / all_area + low;
    }

    return 0;
}

int FittingParameters::makePeaks(vector<double> &x, DataObj *d){
    double a, p, g, l;

    vector<CVP> param(peaks_size);
    for (int k = 0; k < peaks_size; k++){
        p = parameters[k].position;
        a = 1.0;
        l = parameters[k].lorentz_width;
        g = parameters[k].gauss_width;
        CVP pp(a, p, g, l);
        param[k] = pp;
    }

    makeConvVoigt(d, x, param, peaks, signal);

    // double vmax, vmin;
    // for (int k = 0; k < peaks_size; k++) {
    //     vmax = *max_element(signal.begin(), signal.end());
    //     for(auto& s : peaks[k]){ s /= vmax; }
    // }

    // for (int k = 0; k < peaks_size; k++) {
    //     for(auto& s : peaks[k]){ s *= parameters[k].area; }
    // }

    double vmax;
    for (int k = 0; k < peaks_size; k++) {
        vmax = *max_element(peaks[k].begin(), peaks[k].end());
        vmax = max(vmax, 1e-8);
        for(auto& s : peaks[k]){ s *= (parameters[k].area / vmax); }
    }


    return 0;
}

// int FittingParameters::resizer(){
//     lorentz_width.resize(peaks_size);
// 	position.resize(peaks_size);
// 	sigma.resize(peaks_size);

//     return 0;
// }

// int SolutionStruct::calc_log_prior_density(InfoGenetic &info){

//     log_prior_density = 0.0;
//     double tau, ii, jj, si;
//     int num_peaks = common_parameters.peaks_size;

//     // log prior density ; shift
//     for (auto& s : series_list) {
//         log_prior_density += neg_log_prior_Gauss(s.shift, 0.0, info.shift_stddev);
//     }

//     for (int k = 0; k < num_peaks; k++){
//         tau = ( 1. / common_parameters.sigma[k] );
//         for ( auto& ni : node ){
//             ii = series_list[ni.pair[0]].intensity[k];
//             jj = series_list[ni.pair[1]].intensity[k];
//             si = ( tau * tau / 2.0 ) * SQR( ii - jj );
//             log_prior_density += si;
//         }
//     }

//     return 0;
// }

int SolutionStruct::calc_log_likelihood(InfoGenetic &info){
    series.makeSumSpectrum(&info.d);
    series.makeLocalBackground();

    double ei, yi, fi;
    int ndata = series.data.size;
    error = 0.0;
    for ( int i = 0; i < ndata; i++ ) {
        double yi = series.data.y[i];
        double fi = ( series.fitting[i] + series.background[i] );
        
        // ei = ( fi - yi * log(fi) + gsl_sf_lngamma( (yi + 1.0) ) );
        // error += ei;

        if(info.noiseType == "poisson"){
            error += ( fi - yi * log(fi) + gsl_sf_lngamma( (yi + 1.0) ) );
        } else if (info.noiseType == "gauss"){
            error += (fi - yi) * (fi - yi);
        } else {
		    cout << "Error. Undefined noise type " << info.noiseType << endl;
		    exit(1);
        }
	}
    return 0;
}




int SolutionStruct::calc_energy(InfoGenetic &info){

    calc_log_likelihood(info);
    // calc_log_prior_density(info);
    energy = error;

    return 0;
}

int Series::outputFitting(const string ofname, DataObj *d){
    // string ofname;
	ofstream ofs;
    stringstream ss;

    makeSumSpectrum(d);

    double step = abs( data.x[0] - data.x[1] );
    int num_peaks = fitting_parameters.peaks_size;

    ss << "Binding_Energy(eV),spectrum,fitting,background";
	for (int k = 0; k < num_peaks; k++){
        ss << ",peaks[" << k << "]";
    }

	ss << endl;
	for (int i = 0; i < data.size; i++) {
		ss << data.x[i]
			<< "," << data.y[i]
			<< "," << ( fitting[i] + background[i] )
            << "," << background[i];
		for (int k = 0; k < num_peaks; k++){
			ss << "," << fitting_parameters.peaks[k][i] + background[i];
		}
		ss << endl;
	}

	// ofname = kwd;
	ofs.open(ofname);
	if (ofs) {
		// cout << "Output to " + ofname << endl;
		ofs << ss.str();
		ofs.close();
	}
	else {
		cout << "\n!!! Failed to create a file: " << ofname << " !!!\n" << endl;
	}
	ss.str("");

	return 0;
}

int FittingParameters::outputCommonParameters(const string ofname){
    // // string ofname;
	// ofstream ofs;
    // stringstream ss;

    // // ss << "compound,shift,gauss_width,lorentz_width,sigma";
    // ss << "compound,height,position,gauss_width,lorentz_width";
    // ss << endl;

    // for (int k=0; k < peaks_size; k++){
    //     ss << k << ","
    //         << parameters[k].area << ","
    //         << parameters[k].position << ","
    //         << parameters[k].gauss_width << ","
    //         << parameters[k].lorentz_width;
    //     ss << endl;
    // }


	// // ofname = kwd;
	// ofs.open(ofname);
	// if (ofs) {
	// 	// cout << "Output to " + ofname << endl;
	// 	ofs << ss.str();
	// 	ofs.close();
	// }
	// else {
	// 	cout << "\n!!! Failed to create a file: " << ofname << " !!!\n" << endl;
	// }
	// ss.str("");

    FILE *fp;
	if( (fp = fopen(ofname.c_str(), "w") ) == NULL ) {
		printf("Error. Cannot open %s\n", ofname.c_str());
        exit(1);
	}
	fprintf(fp, "compound,area,position,gauss_width,lorentz_width\n");
    for (int k=0; k < peaks_size; k++){
        fprintf(fp, "%d,%.5e,%.3f,%.3f,%.3f\n", k, parameters[k].area, parameters[k].position, parameters[k].gauss_width, parameters[k].lorentz_width);
    }
	fclose(fp);

	return 0;
}

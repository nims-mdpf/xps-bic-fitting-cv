#define _USE_MATH_DEFINES
#define _CRT_SECURE_NO_WARNINGS
#include "optimizeparameters.h"
#include <complex.h>
#include "cerf.h"

#include <gsl/gsl_multimin.h>
#include <gsl/gsl_sf_gamma.h>

Spectrum spect;

void Spectrum::setdata(const vector<double> &xi, const vector<double> &yi)
{
	ndata = xi.size();
	x = xi;
	y = yi;
	xmax = *max_element(x.begin(), x.end());
	xmin = *min_element(x.begin(), x.end());
	ymax = *max_element(y.begin(), y.end());
	ymin = *min_element(y.begin(), y.end());
}

void restrictPeaks(vector<PeakParameter*> &peakList)
{
	for (int j = 0; j < peakList.size(); j++) {
		//shinotsuka20180209
		//peakList[j]->height = max(peakList[j]->height, 0.0);
		//peakList[j]->width = max(peakList[j]->width, 0.0);
		peakList[j]->height = max(peakList[j]->height, 1e-15);
		peakList[j]->sigma = max(peakList[j]->sigma, 1e-15);
		peakList[j]->gamma = max(peakList[j]->gamma, 1e-15);
		//peakList[j]->width = min(peakList[j]->width, 5.0);   //shinotsuka20180219, test. not used
		// peakList[j]->lorentzRatio = max(peakList[j]->lorentzRatio, settings.lorentzRatioCutoffs[0]);
		// peakList[j]->lorentzRatio = min(peakList[j]->lorentzRatio, settings.lorentzRatioCutoffs[1]);
	}
}

// int removeInfinitesimalPeaks(vector<PeakParameter*> &peakList)
// {
// 	int numRemovedPeaks = 0;
// 	const double threshold = 1e-14;  // Threshold of the peak area to be removed

// 	if (peakList.size() == 0) {
// 		return 0;
// 	}

// 	for (int j = (int)peakList.size() - 1; j >= 0; j--) {
// 		if (peakList[j]->width <= threshold || peakList[j]->height <= threshold) {
// 			cout << "Remove a peakList[" << j << "] "
// 				<< peakList[j]->height << " " << peakList[j]->position
// 				<< " " << peakList[j]->width << " " << peakList[j]->lorentzRatio
// 				<< " because its height or width is infinitesimal small." << endl;
// 			peakList.erase(peakList.begin() + j);
// 			numRemovedPeaks++;
// 		}
// 	}
// 	return numRemovedPeaks;
// }

//shinotsuka20180530
//
// Remove very small peaks relative to a whole peak area
// threshold : Threshold of the peak area to be removed
//

int removeSmallPeaksByWholeArea(vector<PeakParameter*> &peakList, const double threshold = 0.01)
{
	if(peakList.size()==1) return 0;
	int numRemovedPeaks = 0;
	// Remove a peak whose area is much less than the maximum peak area.
	vector<double> areas(peakList.size());
	double stdArea = 0.0;
	for (int j = 0; j < peakList.size(); j++) {
		// Peak area of a Voigt function
		areas[j] = voigtArea(peakList[j]);
		// areas[j] = peakList[j]->height / voigt(0.0, peakList[j]->sigma, peakList[j]->gamma);
		stdArea += areas[j];
		// cout << "PeakList[" << j << "] area = " << areas[j] << endl;
	}
	for (int j = (int)peakList.size() - 1; j >= 0; j--) {
		if (areas[j] / stdArea < threshold) {
			// cout << "Remove a peakList[" << j << "] "
			// 	<< peakList[j]->height << " " << peakList[j]->position
			// 	<< " " << peakList[j]->sigma << " " << peakList[j]->gamma
			// 	<< " because its area " << areas[j] << " / " << stdArea << " (whole area) < " << threshold << " (threshold)" << endl;
			peakList.erase(peakList.begin() + j);
			numRemovedPeaks++;
		}
		if(peakList.size()==1) break;
	}	
	return numRemovedPeaks;
}

//shinotsuka20180601
//
// Remove very narrow peaks 
//
int removeSmallPeaksByPeakWidth(vector<PeakParameter*> &peakList, const double energy_step, const double threshold_fwhm = 0.02)
{
	if(peakList.size()==1) return 0;
	int numRemovedPeaks2 = 0;
	// TODO. width of Voigt ?

	for (int i = (int)peakList.size() - 1; i >= 0; i--) {
		double width = peakList[i]->sigma + peakList[i]->gamma;
		if (width < 0.5 * max(energy_step, threshold_fwhm)) {
			// cout << "Remove a peakList[" << i << "] by HWHM "
			// 	<< peakList[i]->height << " " << peakList[i]->position
			// 	<< " " << peakList[i]->sigma << " " << peakList[i]->gamma 
			// 	<< " because HWHM < threshold width " << energy_step << " " << threshold_fwhm << endl;
			peakList.erase(peakList.begin() + i);
			numRemovedPeaks2++;
		}
		if(peakList.size()==1) break;
	}
	return numRemovedPeaks2;
}


// int removeSmallPeaks(vector<PeakParameter*> &peakList)
// {
// 	int numRemovedPeaks = 0;
// 	// TODO. area of Voigt ?
// 	// now skip 

// 	// double threshold = 0.01;  // Threshold of the peak area to be removed

// 	// // Remove a peak whose area is much less than the maximum peak area.
// 	// vector<double> areas(peakList.size());
// 	// double stdArea = 0.0;
// 	// for (int j = 0; j < peakList.size(); j++) {
// 	// 	areas[j] = peakList[j]->height * peakList[j]->width * (peakList[j]->lorentzRatio * M_PI + (1.0 - peakList[j]->lorentzRatio) * sqrt(M_PI / log(2.0)));
// 	// 	stdArea = max(stdArea, areas[j]);
// 	// }
// 	// for (int j = (int)peakList.size() - 1; j >= 0; j--) {
// 	// 	if (areas[j] / stdArea < threshold) {
// 	// 		cout << "Remove a peakList[" << j << "] "
// 	// 			<< peakList[j]->height << " " << peakList[j]->position
// 	// 			<< " " << peakList[j]->sigma << " " << peakList[j]->gamma
// 	// 			<< " because its area " << areas[j] << " / " << stdArea << " (std. area) < " << threshold << " (threshold)" << endl;
// 	// 		peakList.erase(peakList.begin() + j);
// 	// 		numRemovedPeaks++;
// 	// 	}
// 	// }
// 	return numRemovedPeaks;
// }

// int removePeaksByHWHM(vector<PeakParameter*>&peakList, const double threshold_width)
// {
// 	int numRemovedPeaks2 = 0;
// 	for (int i = (int)peakList.size() - 1; i >= 0; i--) {
// 		//cout << i << " " << areas[j] << endl;
// 		// if (peakList[i]->width < threshold_width) {
// 		if ((peakList[i]->sigma < threshold_width) && (peakList[i]->gamma < threshold_width)) {
// 			cout << "Remove a peakList[" << i << "] by HWHM"
// 				<< peakList[i]->height << " " << peakList[i]->position
// 				<< " " << peakList[i]->sigma << " " << peakList[i]->gamma << endl
// 				<< " because the HWHM is smaller than threshold width " << endl;
// 			peakList.erase(peakList.begin() + i);
// 			numRemovedPeaks2++;
// 		}
// 		if(peakList.size()==1) break;
// 	}
// 	return numRemovedPeaks2;
// }

int removePeaksByEnergyRange(const vector<double> &bindingEnergy, vector<PeakParameter*> &peakList)
{
	int numRemovedPeaks3 = 0;
	for (int i = (int)peakList.size() - 1; i >= 0; i--) {
		if (peakList[i]->position > bindingEnergy.front()
			|| peakList[i]->position < bindingEnergy.back()) {
			// cout << "Remove a peakList[" << i << "] by energy range "
			// 	<< peakList[i]->height << " " << peakList[i]->position
			// 	<< " " << peakList[i]->sigma << " " << peakList[i]->gamma << endl
			// 	<< " because its position " << peakList[i]->position << " is out of chiSquare calculation range " << endl;
			peakList.erase(peakList.begin() + i);
			numRemovedPeaks3++;
		}
		if(peakList.size()==1) break;
	}
	return numRemovedPeaks3;
}





/*
* ---------------------------------------------------------------------
* Make a background with bgHigher, bgLower and ffp[]
* INPUT
*    bgHigher:           background height at the higher binding energy side
*    bgLower:            background height at the lower binding energy side
*    ffp:                superposition of the peaks free from any background
* OUTPUT
*    qpq:                Q(x)/(P(x)+Q(x)) in the Shirley equation
*    background:         The Shirley background
* ---------------------------------------------------------------------
*/
int make_background(const double &bgHigher, const double &bgLower,
	const vector<double> &ffp, 
	vector<double> &qpq, vector<double> &background)
{
	int ret = -1;
	if (settings.backgroundType == "Shirley") {
		ret = make_background_Shirley(bgHigher, bgLower, ffp, qpq, background);
	}
	else if (settings.backgroundType == "Linear") {
		// Linear background using two points {bgHigher, bgLower}
		background.resize(ffp.size());
		ret = make_background_Linear(bgHigher, bgLower, background);
	}
	else {
		cout << "Error at make_background(). Unknown settings.backgroundType = " << settings.backgroundType << endl;
		exit(1);
	}
	return ret;
}

int make_background_Shirley(const double &bgHigher, const double &bgLower,
	const vector<double> &ffp, vector<double> &qpq, vector<double> &background)
{
	const int ndata = (int)ffp.size();
	if (ndata == 0) {
		return 1;
	}
	qpq.resize(ndata);
	background.resize(ndata);
	vector<double> areaQ(ndata);
	areaQ[ndata - 1] = 0.0;
	for (int i = ndata - 2; i >= 0; i--) {
		areaQ[i] = areaQ[i + 1] + (ffp[i] + ffp[i + 1]) * 0.5;
	}

	double sumNew = areaQ[0];
	// Background as functions of bgHigher and bgLower
	for (int i = 0; i < ndata; i++) {
		qpq[i] = areaQ[i] / sumNew;
		background[i] = (bgHigher - bgLower) * qpq[i] + bgLower;
	}
	return 0;
}

int make_background_Linear(const double &yH, const double &yL, vector<double> &background)
{
	const int ndata = (int)background.size();
	for (int i = 0; i < ndata; i++) {
		background[i] = ((double)(-i) * yL + (double(i - ndata + 1)) * yH) / (double)(-ndata + 1);
	}
	return 0;
}

/*
* ---------------------------------------------------------------------
* Calculate a chi square value.
* INPUT
*    x[]:                Binding energies (eV)
*    spectrum[]:         Raw spectrum
*    peakList[]:         The peak parameters
*    bgHigher:           background height at the higher binding energy side
*    bgLower:            background height at the lower binding energy side
* OUTPUT
*    qpq[]:              Q(x)/(P(x)+Q(x)) in the Shirley equation
*    s[]:                A residual spectrum
*    newChiSquare:       The chi square value
* ---------------------------------------------------------------------
*/
int calcChiSquare(const vector<double> &x, const vector<double> &spectrum,
	const vector<PeakParameter*> &peakList, 
	const double &bgHigher, const double &bgLower, 
	vector<double> &background,	vector<double> &qpq, vector<double> &s, 
	double &newChiSquare, double &newReducedChiSquare)
{
	const int ndata = (int)x.size();

	if (ndata != spectrum.size()) {
		cout << "Error. Invalid array size, x[] and spectrum[]. " << ndata << " " << spectrum.size() << endl;
		return 1;
	}

	// Update a superposition of the peaks free from any background
	vector<double> ffp;
	if (make_spectrum(x, peakList, ffp)) {
		cout << "Error. Failed to make a superposed spectrum (calcChiSquare)." << endl;
		return 1;
	}
	if (make_background(bgHigher, bgLower, ffp, qpq, background)) {
		cout << "Error. Failed to make a background (calcChiSquare)." << endl;
		return 1;
	}

	// Vector s[]
	s.resize(ndata);

	for (int i = 0; i < ndata; i++) {
		s[i] = (ffp[i] + background[i]) - spectrum[i];
	}

	newChiSquare = 0.0;
	//for (int i = idHigher; i <= idLower; i++) {
	for (int i = 0; i < ndata; i++) {
		newChiSquare += pow(s[i], 2.0) / (ffp[i] + background[i]);
	}

	//int ndata_for_chi2 = (idLower - idHigher + 1);
	int ndata_for_chi2 = ndata;
	int ndf_for_chi2 = ndata_for_chi2 - (4 * (int)peakList.size() + 2);
	newReducedChiSquare = newChiSquare / ndf_for_chi2;

	// This is time consuming.
	//outputChiSquare(x, spectrum, ffp, background, s);

	//if (newChiSquare < 0.0) {
	//	cout << "dbg20190809b, negative chi-square" << endl;
	//	cout << "bgHigher, bgLower = " << bgHigher << ", " << bgLower << endl;
	//	for (int i = 0; i < ndata; i++) {
	//		cout << s[i] << " " << ffp[i] << " " << background[i] << endl;
	//	}
	//}

	return 0;
}

void outputChiSquare(const vector<double> &x, const vector<double> &spectrum, const vector<double> &ffp, const vector<double> &bg, const vector<double> &s)
{
	static int call_count = 0;
	const int ndata = (int)x.size();
	string ofname;
	ofstream ofs;
	stringstream ss;

	ss << "x,spectrum,background,fitted,chiSquare" << endl;
	string kwd = "output_chiSquare";

	ofname = kwd + ".csv";

	//shinotsuka20181122
	char ofnamec[1024];
	sprintf(ofnamec, "output_chiSquare_%05d.csv", ++call_count);
	ofname = string(ofnamec);
	//cout << ofname << endl;
	//return;

	ofs.open(ofname);
	for (int i = 0; i < ndata; i++) {
		ss << x[i] 
			<< "," << spectrum[i]
			<< "," << bg[i]
			<< "," << ffp[i] + bg[i]
			<< "," << pow(s[i] / spectrum[i], 2.0) << endl;
	}
	ss << endl;

	if (ofs) {
		cout << ofname << ", call_count=" << call_count << endl;
		ofs << ss.str();
		ofs.close();
	}
	else {
		cout << "\n!!! Failed to create a file: " << ofname << " !!!\n" << endl;
	}
}

int output_history(const int m_SG, vector<double> x, vector<double> y, vector<double> background, vector<PeakParameter*> peakList)
{
	string ofname;
	ofstream ofs;
	static int call_count = 0;
	const int ndata = (int)x.size();
	const int originalNdata = (int)x.size();
	cout << "<< m_SG = " << m_SG << "  history No. " << call_count << ">>" << endl;

	// Output a list of peak parameters and their areas
	cout << "Number of peaks=" << peakList.size() << endl;

	if (peakList.size() == 0) {
		cout << "Skip to output file because there is no peaks." << endl;
		return 0;
	}

	//shinotsuka20181205a
	// Before output the result, denormalize spectra and parameters
	//denormalize_spectra_and_parameters(x, y, background, peakList,
	//	x, y, smoothedSpectrum, settings);

	stringstream ss;
	ss.str("");


	// Output fitted spectra
	vector<double> sumSpectrum;
	vector<vector<double>> peaks;
	make_spectrum_and_peak(x, peakList, sumSpectrum, peaks);

	ss << "Binding_Energy(eV),spectrum,background,fitting";
	for (int j = 0; j < peakList.size(); j++) { ss << ",peak[" << j << "]"; }
	ss << endl;

	for (int i = 0; i < x.size(); i++) {
		ss << x[i]
			<< "," << y[i]
			<< "," << background[i]
			<< "," << background[i] + sumSpectrum[i];
		for (int j = 0; j < peakList.size(); j++) {
			ss << "," << background[i] + peaks[j][i];
		}
		ss << endl;
	}


	//ofname = "result_" + kwd + ".csv";
	char ofnamec[1024];
	//sprintf(ofnamec, "optimize_history_No_%05d.csv", ++call_count);
	sprintf(ofnamec, "optimize_history_mSG_%04d_No_%05d.csv", m_SG, ++call_count);
	ofname = string(ofnamec);

	ofs.open(ofname);
	if (ofs) {
		cout << "Output to " + ofname << endl;
		ofs << ss.str();
		ofs.close();
	}
	else {
		cout << "\n!!! Failed to create a file: " << ofname << " !!!\n" << endl;
	}
	ss.str("");

	return 0;
}



/*
* ---------------------------------------------------------------------
* Copy from {peakList_in, bgHigher_in, bgLower_in} 
*        to {peakList_out, bgHigher_out, bgLower_out}
* INPUT
*    peakList_in[]:         Peak parameters
*    bgHigher_in:           Background height at the higher binding energy side
*    bgLower_in:            Background height at the lower binding energy side
* OUTPUT
*    peakList_out[]:        Peak parameters
*    bgHigher_out:          Background height at the higher binding energy side
*    bgLower_out:           Background height at the lower binding energy side
* ---------------------------------------------------------------------
*/
void copy_peakparameter(const vector<PeakParameter*> &peakList_in, const double bgHigher_in, const double bgLower_in,
	vector<PeakParameter*> &peakList_out, double &bgHigher_out, double &bgLower_out)
{
	//shinotsuka20190807
	copyVectorT(peakList_in, peakList_out);
	bgHigher_out = bgHigher_in;
	bgLower_out = bgLower_in;
}



// double my_f (const gsl_vector *v, void *params)
// {
// 	double x, y;
// 	double *p = (double *)params;

// 	x = gsl_vector_get(v, 0);
// 	y = gsl_vector_get(v, 1);

// 	return p[2] * (x - p[0]) * (x - p[0]) +
// 	p[3] * (y - p[1]) * (y - p[1]) + p[4];
// }

// vector<double> xxx;
// vector<double> yyy;
// Settings settingsCopy;

void convert_gsl_vector_to_peakparameter(const gsl_vector *v, vector<PeakParameter*> &peakList, double &bgHigher, double &bgLower)
{
	int numpeak = (v->size - 2) / 4;
	vector<PeakParameter*>().swap(peakList);
	// printf("calc_poisson_error. numpara = %d, numpeak = %d\n", numpara, numpeak);
	for (int j = 0; j < numpeak; j++) {
		double height = gsl_vector_get(v, 0+j*4);
		double position = gsl_vector_get(v, 1+j*4);
		double sigma = gsl_vector_get(v, 2+j*4);
		double gamma = gsl_vector_get(v, 3+j*4);
		//cout << height << "," << position << "," << width << "," << lorentzRatio << endl;
		peakList.push_back(new PeakParameter(height, position, sigma, gamma));
	}

	bgHigher = gsl_vector_get(v, numpeak*4);
	bgLower = gsl_vector_get(v, 1+numpeak*4);
	//cout << bgHigher << "," << bgLower << endl; getchar();
}

double calc_log_likelihood_Poisson(const vector<double> &s, const vector<double> &y){
	// s is fitted data, float
	// y is measured data, int
    double error = 0.0;
    int ndata = (int)s.size();
    for ( int i = 0; i < ndata; i++ ) {
        error += ( s[i] - y[i] * log(s[i]) + gsl_sf_lngamma( (y[i] + 1.0) ) );
    }
    return error;
}
double calc_log_likelihood_Gauss(const vector<double> &s, const vector<double> &y){
	// s is fitted data, float
	// y is measured data, float
    double error = 0.0;
    int ndata = (int)s.size();
    for ( int i = 0; i < ndata; i++ ) {
        error += ( s[i] - y[i] ) * ( s[i] - y[i] );
    }
	// error = sqrt(error/ndata);
    return error;
}

double calc_poisson_error(const gsl_vector *v, void *params)
{
	// cout << v->size << endl;
	// const int ndata = xxx.size();
	const int ndata = spect.ndata;
	vector<double> ffp;
	vector<double> qpq;
	vector<double> background;
	vector<double> s;
	ffp.resize(ndata);

	vector<PeakParameter*> peakList;
	double bgHigher;
	double bgLower;
	convert_gsl_vector_to_peakparameter(v, peakList, bgHigher, bgLower);

	// if (make_spectrum(xxx, peakList, ffp)) {
	if (make_spectrum(spect.x, peakList, ffp)) {
		cout << "Error. Failed to make a superposed spectrum (calc_poisson_error)." << endl;
		return 1;
	}
	background.resize(ndata);
	if (make_background(bgHigher, bgLower, ffp, qpq, background)) {
		cout << "Error. Failed to make a background (calc_poisson_error)." << endl;
		return 1;
	}

	// Vector s[]
	s.resize(ndata);
	for (int i = 0; i < ndata; i++) {
		s[i] = ffp[i] + background[i];
		// Protect not to be < 0
		s[i] = max(s[i], 1e-15);
	}

	double errorfunc = 0.0;
	if(settings.noiseType == "Poisson"){
		// errorfunc = calc_log_likelihood_Poisson(s, yyy);
		errorfunc = calc_log_likelihood_Poisson(s, spect.y);
	} else if(settings.noiseType == "Gauss"){
		//errorfunc = calc_log_likelihood_Gauss(s, yyy);
		errorfunc = calc_log_likelihood_Gauss(s, spect.y);
	} else {
		cout << "Error. Undefined noise type " << settings.noiseType << endl;
		exit(1);
	}
	// errorfunc = calc_log_likelihood(s, yyy);

	//printf("errorfunc = %f\n", errorfunc);
	return errorfunc;
}

/*
* ---------------------------------------------------------------------
* Calculate convolution Voigt function V(x,sigma,gamma) and its derivatives
* INPUT
*	x:		position
*	sigma:	STD of Gauss function
*	gamma:	half width of Lorentz function
* OUTPUT
*	v:		V(x,sigma,gamma)
*	dvdx:	dV(x,sigma,gamma)/dx
*	dvds:	dV(x,sigma,gamma)/d(sigma)
*	dvdg:	dV(x,sigma,gamma)/d(gamma)
* ---------------------------------------------------------------------
*/
void voigt_derivs(double x, double s, double g, double &v, double &dvdx, double &dvds, double &dvdg){
    double value = 0.0;
    double _Complex z, wofz, ztemp;
	s = max(s, 1e-8);
	g = max(g, 1e-8);
    z = (x + I * g) / (sqrt(2.0) * s);
    wofz = w_of_z(z);
	double ssquare = s * s;
    double scube = s * ssquare;
	double sq2pi = sqrt(2.0 * M_PI);

    v = creal(wofz) / (sq2pi * s);

    ztemp = wofz * (x - s + I * g) * (x + s + I * g);
    dvds = g / (M_PI * scube) + creal(ztemp) / (sq2pi * s * scube);

    ztemp = wofz * (g - I * x);
    dvdg = -1.0 / (M_PI * ssquare) + creal(ztemp) / (sq2pi * scube);

    ztemp = wofz * (x + I * g);
    dvdx = - creal(ztemp) / (sq2pi * scube);

}

/* The gradient of f, df = (df/dx, df/dy). */
void
my_df (const gsl_vector *v, void *params,
       gsl_vector *df)
{
	const int ndata = spect.ndata;
	vector<double> ffp;
	vector<double> qpq;
	vector<double> background;
	vector<double> func;
	ffp.resize(ndata);

	vector<PeakParameter*> peakList;
	double bgHigher;
	double bgLower;
	convert_gsl_vector_to_peakparameter(v, peakList, bgHigher, bgLower);

	if (make_spectrum(spect.x, peakList, ffp)) {
		cout << "Error. Failed to make a superposed spectrum (my_df)." << endl;
		exit(1);
	}
	background.resize(ndata);
	if (make_background(bgHigher, bgLower, ffp, qpq, background)) {
		cout << "Error. Failed to make a background (my_df)." << endl;
		exit(1);
	}

	// Vector func[]
	func.resize(ndata);
	for (int i = 0; i < ndata; i++) {
		func[i] = ffp[i] + background[i];
		// Protect not to be < 0
		func[i] = max(func[i], 1e-15);
	}

	gsl_vector_set_zero(df);
	int cnt = 0;
	for (int j = 0; j < peakList.size(); j++) {
		double mu, h, s, g;
		mu = peakList[j]->position;
		h = peakList[j]->height;
		s = peakList[j]->sigma;
		g = peakList[j]->gamma;

		// Calculate V(0,s,g), dVdx(0,s,g), dVds(0,s,g), dVdg(0,s,g)
		double v0, dvdx0, dvds0, dvdg0;
		voigt_derivs(0.0, s, g, v0, dvdx0, dvds0, dvdg0);

		double dedm = 0.0;
		double deda = 0.0;
		double deds = 0.0;
		double dedg = 0.0;
		for(int i=0; i<ndata; i++){
			double x = spect.x[i];
			double y = spect.y[i];
			// double v = a * voigt(x - mu, s, g);
			// double v = voigt(x - mu, s, g);

			// Calculate V(x-mu,s,g), dVdx(x-mu,s,g), dVds(x-mu,s,g), dVdg(x-mu,s,g)
			double v, dvdx1, dvds1, dvdg1;
			voigt_derivs(x - mu, s, g, v, dvdx1, dvds1, dvdg1);

			///////////////////////////////////////////
			double f = func[i];
			// dedm += (2.0) * (y - f) * h / v0 * dvdx(x - mu, s, g);
			// deda += (-2.0) * (y - f) * v / v0;
			// deds += (-2.0) * (y - f) * h * (dvds(x - mu, s, g) * v0 - v * dvds0) / (v0 * v0);
			// dedg += (-2.0) * (y - f) * h * (dvdg(x - mu, s, g) * v0 - v * dvdg0) / (v0 * v0);
			dedm += (2.0) * (y - f) * h / v0 * dvdx1;
			deda += (-2.0) * (y - f) * v / v0;
			deds += (-2.0) * (y - f) * h * (dvds1 * v0 - v * dvds0) / (v0 * v0);
			dedg += (-2.0) * (y - f) * h * (dvdg1 * v0 - v * dvdg0) / (v0 * v0);
		}
		gsl_vector_set(df, cnt++, dedm);
		gsl_vector_set(df, cnt++, deda);
		gsl_vector_set(df, cnt++, deds);
		gsl_vector_set(df, cnt++, dedg);

	}

	double deda = 0.0;
	double dedb = 0.0;
	for(int i=0; i<ndata; i++){
		double y = spect.y[i];
		double f = func[i];
		deda += (-2.0) * (y - f) * qpq[i];
		dedb += (-2.0) * (y - f) * (1.0 - qpq[i]);
	}
	gsl_vector_set(df, cnt++, deda);
	gsl_vector_set(df, cnt++, dedb);

	// printf("debug. cnt = %d, numPeak = %d\n", cnt, peakList.size());
}



/* Compute both f and df together. */
void
my_fdf (const gsl_vector *x, void *params,
        double *f, gsl_vector *df)
{
  *f = calc_poisson_error(x, params);
  my_df(x, params, df);
}


void restrict_parameter(gsl_vector *x)
{

	int numpeak = (x->size - 2) / 4;
	for (int j = 0; j < numpeak; j++) {
		// Restrict height
		gsl_vector_set (x, 0+j*4, max(1e-15, gsl_vector_get(x, 0+j*4)));

		// Restrict position
		double p = gsl_vector_get(x, 1+j*4);
		// p = min(max(xxx.back(), p), xxx.front());
		p = min(max(spect.xmin, p), spect.xmax);
		gsl_vector_set (x, 1+j*4, p);

		// Restrict sigma
		gsl_vector_set (x, 2+j*4, max(1e-8, gsl_vector_get(x, 2+j*4)));
		// gsl_vector_set (x, 2+j*4, max(0.0, gsl_vector_get(x, 2+j*4)));
		// Restrict gamma
		gsl_vector_set (x, 3+j*4, max(1e-8, gsl_vector_get(x, 3+j*4)));
		// gsl_vector_set (x, 3+j*4, max(0.0, gsl_vector_get(x, 3+j*4)));

		// if((gsl_vector_get(x, 2+j*4) == 0.0) && (gsl_vector_get(x, 3+j*4) == 0.0)){
		// 	// If both widths are zero, set pure lorentzian with small gamma
		// 	gsl_vector_set (x, 2+j*4, 0.0);
		// 	gsl_vector_set (x, 3+j*4, 1e-8);
		// }

	}

	// Restrict background parameters
	gsl_vector_set (x, x->size-2, max(gsl_vector_get(x, x->size-2), 1e-15));
	gsl_vector_set (x, x->size-1, max(gsl_vector_get(x, x->size-1), 1e-15));

	if (settings.backgroundType == "Shirley") {
		// background at higer BE side > lower side
		gsl_vector_set (x, x->size-2, max(gsl_vector_get(x, x->size-2), gsl_vector_get(x, x->size-1)));
	}



	// cout << "After restriction" << endl;
	// //gsl_vector_fprintf(stdout, x, "%f, "); getchar();
	// for(int i=0; i<x->size; i++){cout << gsl_vector_get(x, i) << ", ";}
	// cout << endl; //getchar();
}


int optimizeParametersGSL(const int ndata, const int m_SG, const vector<double> &x, const vector<double> &spectrum,
	const bool isPeakOptimize, const bool isBgOptimize, const bool isFixedLorentzRatio, const double fixedLorentzRatio,
	vector<double> &background, vector<PeakParameter*> &peakList,
	bool isNumberPeakFixed, int loopMaxOptimisze,
	double &optimizedChiSquare, double &optimizedReducedChiSquare,
	double &bgHigher, double &bgLower, int &loop)
{
	// cout << "optimizeParametersGSL() Initial number of peaks = " << peakList.size() << endl;
	// for(auto pp : peakList){ cout << pp->height << ", " << pp->position << ", " << pp->sigma << ", " << pp->gamma << endl; }
	// cout << bgHigher << ", " << bgLower << endl;

	int num_parameter = 0;
	if (isPeakOptimize) { num_parameter += (int)peakList.size() * 4; }
	if (isBgOptimize) { num_parameter += 2; }

	if (num_parameter == 0) {
		cout << "Warning: There is no paraemter to be optimized." << endl;
		return 0;
	}

	// set global variables
	// xxx = x;
	// yyy = spectrum;
	// settingsCopy = settings;

	double *par;

	const gsl_multimin_fminimizer_type *T =	gsl_multimin_fminimizer_nmsimplex2;
	// const gsl_multimin_fminimizer_type *T =	gsl_multimin_fminimizer_nmsimplex;
	// const gsl_multimin_fminimizer_type *T =	gsl_multimin_fminimizer_nmsimplex2rand;
	
	gsl_multimin_fminimizer *s = NULL;
	gsl_vector *ss, *gslx;
	gsl_multimin_function minex_func;

	size_t iter = 0;
	size_t iterMax = 1000;
	int numNotMove = 0;
	// int numNotMoveLimit = 200;
	int numNotMoveLimit = 100;
	double fvalold = 1e301;
	double fvalnew = 1e300;

	int status;
	double size;

	/* Starting point */
	gslx = gsl_vector_alloc (num_parameter);
	int cnt = -1;
	for (int j = 0; j < peakList.size(); j++) {
		gsl_vector_set (gslx, ++cnt, peakList[j]->height);
		gsl_vector_set (gslx, ++cnt, peakList[j]->position);
		gsl_vector_set (gslx, ++cnt, peakList[j]->sigma);
		gsl_vector_set (gslx, ++cnt, peakList[j]->gamma);
	}
	gsl_vector_set (gslx, ++cnt, bgHigher);
	gsl_vector_set (gslx, ++cnt, bgLower);
	

	/* Set initial step sizes */
	/* 探索範囲を初期値の1.0倍とする，すなわち，0倍から2倍の範囲を探索する */
	ss = gsl_vector_alloc (num_parameter);
	gsl_vector_memcpy( ss, gslx );
	gsl_vector_scale( ss, 1.0 );
	// gsl_vector_scale( ss, 0.5 );

	// 各パラメータの探索範囲の初期設定
	// ピーク高さ：スペクトルのyの範囲の0.1倍
	// ピーク位置：スペクトルのxの範囲の0.1倍
	// ガウス幅とローレンツ幅：エネルギー間隔*1.0または現在の値*1.0のうち大きいほう
	// BG端点高さ：yの範囲の0.05倍または現在の値の0.5倍のうち小さいほう
	double width_min_step = 1.0 * abs(x[0] - x[1]);
	double position_step = abs(x.front() - x.back()) * 0.1;
	double height_step = settings.spectrum_range * 0.1;
	// double bgH_step = gsl_vector_get(ss, );
	for (int j = 0; j < peakList.size(); j++) {
		gsl_vector_set (ss, 0+j*4, height_step);
		gsl_vector_set (ss, 1+j*4, position_step);
		gsl_vector_set (ss, 2+j*4, max(width_min_step, gsl_vector_get(ss, 2+j*4)));
		gsl_vector_set (ss, 3+j*4, max(width_min_step, gsl_vector_get(ss, 3+j*4)));
	}
	double bgH_step = min(gsl_vector_get(gslx, gslx->size-2)*0.5, height_step*0.5);
	double bgL_step = min(gsl_vector_get(gslx, gslx->size-1)*0.5, height_step*0.5);
	gsl_vector_set (ss, ss->size-2, bgH_step);
	gsl_vector_set (ss, ss->size-1, bgL_step);


	//gsl_vector_fprintf(stdout, ss, "%f, "); getchar();


	/* Initialize method and iterate */
	minex_func.n = num_parameter;
	minex_func.f = calc_poisson_error;
	minex_func.params = par;

	s = gsl_multimin_fminimizer_alloc (T, num_parameter);
	gsl_multimin_fminimizer_set (s, &minex_func, gslx, ss);

	do
	{
		iter++;

		// cout << "Check" << endl; for(int i=0; i<s->x->size; i++){cout << gsl_vector_get(s->x, i) << ", ";}; cout << endl; 

		status = gsl_multimin_fminimizer_iterate(s);

		if (status)
			break;

		size = gsl_multimin_fminimizer_size (s);
		status = gsl_multimin_test_size (size, 1e-2);

		// printf("NM. %5d, %12.3f, %8.3f", (int)iter, s->fval, size); for(int i=0; i<s->x->size; i++){printf(", %8.3f", gsl_vector_get(s->x, i));}; cout << endl;

		// if (status == GSL_SUCCESS) {
		// 	printf ("converged to minimum at\n");
		// }
		// printf("%5d f() = %7.3f size = %.3f\n", (int)iter, s->fval, size);

		// fvalnew = s->fval;

		/* 境界条件のせいで収束判定が機能しない．そこで境界条件の制約を課してから関数値を計算 */
		restrict_parameter(s->x);
		fvalnew = calc_poisson_error(s->x, par);
		double relerr = abs((fvalnew - fvalold) / fvalold);
		// printf("%5d f() = %7.3f size = %.3f, relerr = %.5f\n", (int)iter, fvalnew, size, relerr);

		/* 指定回数動かなければ収束したとみなす */
		// if(fvalnew == fvalold){
		if(relerr < 1e-5){
			numNotMove++;
		} else {
			fvalold = fvalnew;
			numNotMove = 0;
		}
		if(numNotMove >= numNotMoveLimit){
			cout << "Take it converged because it does not move " << numNotMoveLimit << " times" << endl;
			break;
		}

	}
	while (status == GSL_CONTINUE && iter < iterMax);

	if (status == GSL_SUCCESS){
		printf ("GSLNM. Converged to minimum at ");
	} else {
		printf ("GSLNM. Not converged to minimum at ");
	}
	printf("%5d f() = %7.3f size = %.3f\n", (int)iter, fvalnew, size);

	convert_gsl_vector_to_peakparameter(s->x, peakList, bgHigher, bgLower);
	
	// for(auto& p : peakList){
	// 	cout << p->height << ", " << p->position << ", " << p->sigma << ", " << p->gamma << ", ";
	// }
	// cout << bgHigher << ", " << bgLower << endl;

	// cout << "step size " << endl;
	// for(int i=0; i<ss->size; i++){cout << gsl_vector_get(ss, i) << ", ";}
	// cout << endl;


	/* 最後のrestrict実行により，変数および関数値が変化しうることに注意 */
	// cout << "optimizedChiSquare (before restriction) = " << fvalnew << endl;
	optimizedChiSquare = calc_poisson_error(s->x, par);
	// cout << "optimizedChiSquare (after  restriction) = " << optimizedChiSquare << endl; //getchar();

	gsl_vector_free(gslx);
	gsl_vector_free(ss);
	gsl_multimin_fminimizer_free (s);



	//return status;

	// if(par){ delete[] par; }
	
	return 0;
}


int optimizeParametersCG(const int ndata, const int m_SG, const vector<double> &x, const vector<double> &spectrum,
	const bool isPeakOptimize, const bool isBgOptimize, const bool isFixedLorentzRatio, const double fixedLorentzRatio,
	vector<double> &background, vector<PeakParameter*> &peakList,
	bool isNumberPeakFixed, int loopMaxOptimisze,
	double &optimizedChiSquare, double &optimizedReducedChiSquare,
	double &bgHigher, double &bgLower, int &loop)
{
	// cout << "optimizeParametersCG() Initial number of peaks = " << peakList.size() << endl;
	// for(auto pp : peakList){ cout << pp->height << ", " << pp->position << ", " << pp->sigma << ", " << pp->gamma << endl; }
	// cout << bgHigher << ", " << bgLower << endl;

	int num_parameter = 0;
	if (isPeakOptimize) { num_parameter += (int)peakList.size() * 4; }
	if (isBgOptimize) { num_parameter += 2; }

	if (num_parameter == 0) {
		cout << "Warning: There is no paraemter to be optimized." << endl;
		return 0;
	}

	double *par;

    const gsl_multimin_fdfminimizer_type *T;
    gsl_multimin_fdfminimizer *s;

	gsl_vector *gslx; //*ss, 
	gsl_multimin_function_fdf minex_func;

	size_t iter = 0;
	size_t iterMax = 1000;
	int numNotMove = 0;
	// int numNotMoveLimit = 200;
	int numNotMoveLimit = 20;
	double fvalold = 1e301;
	double fvalnew = 1e300;

	int status;
	// double size;

	/* Starting point */
	gslx = gsl_vector_alloc (num_parameter);
	int cnt = -1;
	for (int j = 0; j < peakList.size(); j++) {
		gsl_vector_set (gslx, ++cnt, peakList[j]->height);
		gsl_vector_set (gslx, ++cnt, peakList[j]->position);
		gsl_vector_set (gslx, ++cnt, peakList[j]->sigma);
		gsl_vector_set (gslx, ++cnt, peakList[j]->gamma);
	}
	gsl_vector_set (gslx, ++cnt, bgHigher);
	gsl_vector_set (gslx, ++cnt, bgLower);
	

	/* Initialize method and iterate */
	minex_func.n = num_parameter;
	minex_func.f = calc_poisson_error;
	minex_func.df = my_df;
	minex_func.fdf = my_fdf;
	minex_func.params = par;

	T = gsl_multimin_fdfminimizer_conjugate_fr;
	// T = gsl_multimin_fdfminimizer_conjugate_pr;
	// T = gsl_multimin_fdfminimizer_vector_bfgs2;
	// T = gsl_multimin_fdfminimizer_vector_bfgs;
	// T = gsl_multimin_fdfminimizer_steepest_descent;

	s = gsl_multimin_fdfminimizer_alloc (T, num_parameter);
	// gsl_multimin_fdfminimizer_set (s, &minex_func, gslx, 0.01, 1e-4);
	// gsl_multimin_fdfminimizer_set (s, &minex_func, gslx, 0.01, 1e-3);
	gsl_multimin_fdfminimizer_set (s, &minex_func, gslx, 10.0, 1e-4);


	do
	{
		iter++;

		// cout << "Check" << endl; for(int i=0; i<s->x->size; i++){cout << gsl_vector_get(s->x, i) << ", ";}; cout << endl; 

		status = gsl_multimin_fdfminimizer_iterate (s);

		// if (status)
		// 	break;
		if (status && iter >= 5)
			break;

		status = gsl_multimin_test_gradient (s->gradient, 1e-3);

		// For debug
		printf("CG. %5d, %12.3f, %5d", (int)iter, s->f, status); for(int i=0; i<s->x->size; i++){printf(", %8.3f", gsl_vector_get(s->x, i));}; cout << endl;

		// fvalnew = s->fval;

		/* 境界条件のせいで収束判定が機能しない．そこで境界条件の制約を課してから関数値を計算 */
		restrict_parameter(s->x);
		fvalnew = calc_poisson_error(s->x, par);
		double relerr = abs((fvalnew - fvalold) / fvalold);
		// printf("%5d f() = %7.3f size = %.3f, relerr = %.5f\n", (int)iter, fvalnew, size, relerr);

		/* 指定回数動かなければ収束したとみなす */
		// if(fvalnew == fvalold){
		if(relerr < 1e-5){
			numNotMove++;
		} else {
			fvalold = fvalnew;
			numNotMove = 0;
		}
		if(numNotMove >= numNotMoveLimit){
			cout << "Take it converged because it does not move " << numNotMoveLimit << " times" << endl;
			break;
		}

	}
	while ((status == GSL_CONTINUE && iter < iterMax) || iter < 5);
	// while (status == GSL_CONTINUE && iter < iterMax);

	if (status == GSL_SUCCESS){
		printf ("GSLCG. Converged to minimum at ");
	} else {
		printf ("GSLCG. Not converged to minimum at ");
	}
	printf("%5d f() = %7.3f \n", (int)iter, fvalnew);

	convert_gsl_vector_to_peakparameter(s->x, peakList, bgHigher, bgLower);
	
	// for(auto& p : peakList){
	// 	cout << p->height << ", " << p->position << ", " << p->sigma << ", " << p->gamma << ", ";
	// }
	// cout << bgHigher << ", " << bgLower << endl;

	// cout << "step size " << endl;
	// for(int i=0; i<s->dx->size; i++){cout << gsl_vector_get(s->dx, i) << ", ";}
	// cout << endl;


	/* 最後のrestrict実行により，変数および関数値が変化しうることに注意 */
	// cout << "optimizedChiSquare (before restriction) = " << fvalnew << endl;
	optimizedChiSquare = calc_poisson_error(s->x, par);
	// cout << "optimizedChiSquare (after  restriction) = " << optimizedChiSquare << endl; //getchar();

	gsl_multimin_fdfminimizer_free (s);
	gsl_vector_free(gslx);

	return 0;
}

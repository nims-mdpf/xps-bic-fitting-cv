#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
//#include <gsl/gsl_rng.h>
#include "tblfunc4.h"

// #include "cerf.h"

#ifdef _MSC_VER
#include <stdarg.h>
void err(int eval, const char *fmt, ...)
{
	printf("Error: Status = %d\n", eval);
	va_list ap;
	va_start(ap, fmt);
	vprintf(fmt, ap);
	va_end(ap);
	printf("\n");
	exit(eval);
}
#else
#include <err.h>
#endif

// #define DEBUGTBLFUNC 1

void getIDandPortion(const double x, const double xmin, const double xmax, const double dx, const int nx, int *id, double *p)
{
	if (x <= xmin) {
		*id = 0;
		*p = 0.0;
	}
	else if (x >= xmax) {
		*id = nx - 2;
		*p = 1.0;
	}
	else {
		*id = (int)((x - xmin) / dx);
		*p = (x - (xmin + (*id) * dx)) / dx;
	}
}

void getIDandPortionLog(const double x, const double xmin, const double xmax, const double dx, const int nx, int *id, double *p)
{
	if (x < xmin) {
		*id = 0;
		*p = 0.0;
	}
	else if (x > xmax) {
		*id = nx - 2;
		*p = 1.0;
	}
	else {
		double tmp = log(x / xmin);
		// *id = (int)(tmp / dx);
		// *p = (tmp - (*id) * dx) / dx;
		*id = (int)(tmp / dx);
		*p = (tmp - (*id) * dx) / dx;
		*id = (*id) + 1;  // xmin and rmin start from vx[1]. vx[0]=0 is reserved.
	}
}

void getIDandPortionNew(const double x, const FuncTbl3 *t, int *id, double *p)
{
	double xp = fabs(x);

	if (xp <= t->xmin) {
		*id = 0;
		*p = 0.0;
	}
	else if (xp >= t->xmax) {
		*id = t->nx - 2;
		*p = 1.0;
		cnt_ftbl_outside++;
	}
	else {
		//*id = (int)(xp*xp / (1.0 + xp*xp) * (t->nx - 1));
		double xsq = xp * xp;
		*id = (int)(xsq / (1.0 + xsq) * (t->nx - 1));
		// between VX[*id], VX[*id+1]
		*p = (xp - t->xval[*id]) / (t->xval[*id + 1] - t->xval[*id]);
	}
}

#ifndef SQR
#define SQR(x) ((x)*(x))
#endif

// 0 <= px <= 1
// 0 <= py <= 1
// f00 = value at (px, py) = (0, 0)
// f10 = value at (px, py) = (1, 0)
// f01 = value at (px, py) = (0, 1)
// f11 = value at (px, py) = (1, 1)
double bilinear_interpolation(double px, double py, double f00, double f10, double f01, double f11)
{
	// printf("%e, %e, %e, %e\n", f00, f01, f10, f11);
	// printf("%e\n", (1.0 - px) * (1.0 - py) * f00
	// 	+ px * (1.0 - py) * f10
	// 	+ (1.0 - px) * py * f01
	// 	+ px * py * f11);
	return (1.0 - px) * (1.0 - py) * f00
		+ px * (1.0 - py) * f10
		+ (1.0 - px) * py * f01
		+ px * py * f11;
}

// Initialize counter
long cnt_ftbl_total = 0;
long cnt_ftbl_outside = 0;

// Output counter
void output_cnt_ftbl()
{
	printf("cnt_ftbl_outide / cnt_ftbl_total = %ld / %ld = %lf\n", cnt_ftbl_outside, cnt_ftbl_total, (double)cnt_ftbl_outside / (double)cnt_ftbl_total);
}

double lorentz(double x, double r)
{
	return (r / M_PI) / (x * x + r * r);
}

double lorentzCDF(double x, double r)
{
	return 0.5 + (1.0 / M_PI) * atan2(x, r);
}

// // Voigt function by libcerf.h
// double ftbl(double x, double c, double s, double r)
// {
// 	double v = voigt(x-c, s, r);
// 	return v;
// }

// Voigt function by table
double ftbl(double x, double c, double s, double r, FuncTbl3 *t)
{
	s = fabs(s);
	r = fabs(r);
	// return ftbl(x, c, s, r);

	// printf("x, c, s, g = %f, %f, %f, %f\n", x, c, s, r);
	cnt_ftbl_total++;

	if(s == 0.0){
		// Pure lorentzian
		return lorentz(x - c, r);
	}

	double xt = (x - c) / s;
	double rt = r / s;
	// printf("xt, rt = %f, %f\n", xt, rt);
	// printf("t->rmax = %f\n", t->rmax);

	if(rt > t->rmax){
		// When rt > rmax, it is almost pure lorentzian with HWHM rt
		// Gaussian with sigma=1 behaves as a delta-function
		// printf("here, %f, %f\n", lorentz(xt, rt), (rt / M_PI) / (xt * xt + rt * rt));
		// return lorentz(xt, rt); // wrong
		return lorentz(x - c, r);
	}

	if (fabs(xt) > t->xmax) { cnt_ftbl_outside++; return 0.0; }

    int idx, idr;
    double px, pr;
	//getIDandPortionNew(x, t, &idx, &px);
	// getIDandPortion(x, t->xmin, t->xmax, t->dx, t->nx, &idx, &px);
	// getIDandPortionLog(x, t->xmin, t->xmax, t->dx, t->nx, &idx, &px);
	getIDandPortionLog(fabs(xt), t->xmin, t->xmax, t->dx, t->nx, &idx, &px);
	getIDandPortionLog(rt, t->rmin, t->rmax, t->dr, t->nr, &idr, &pr);
	
	
#ifdef DEBUGTBLFUNC
	printf("xt = %f, idx/nx = %d/%d, px = %f\n", xt, idx, t->nx, px);
	printf("rt = %f, idr/nr = %d/%d, pr = %f\n", rt, idr, t->nr, pr);
	printf("xt = %f, x[%d] = %f, x[%d] = %f\n", xt, idx, t->xval[idx], idx+1, t->xval[idx+1]);
	printf("rt = %f, r[%d] = %f, r[%d] = %f\n", rt, idr, t->rval[idr], idr+1, t->rval[idr+1]);
	printf("fval= %e, %e, %e, %e\n", t->fval[idr][idx], t->fval[idr + 1][idx], t->fval[idr][idx + 1], t->fval[idr + 1][idx + 1]);
#endif

	// return bilinear_interpolation(px, pr, t->fval[idr][idx], t->fval[idr][idx + 1], t->fval[idr + 1][idx], t->fval[idr + 1][idx + 1]);
	// return (1.0/s) * bilinear_interpolation(px, pr, t->fval[idr][idx], t->fval[idr][idx + 1], t->fval[idr + 1][idx], t->fval[idr + 1][idx + 1]);
	double v = (1.0/s) * bilinear_interpolation(px, pr, t->fval[idr][idx], t->fval[idr][idx + 1], t->fval[idr + 1][idx], t->fval[idr + 1][idx + 1]);

	// Return Voigt function whose peak area = 1
	return v;

	// // Return Voigt function normalized at x=0
	// // V(0, 1, g/s)
	// px = 0.0; idx = 0;
	// double v0 = (1.0/s) * bilinear_interpolation(px, pr, t->fval[idr][idx], t->fval[idr][idx + 1], t->fval[idr + 1][idx], t->fval[idr + 1][idx + 1]);
	// return v/v0;
}

double finttbl(double x, double c, double s, double r, FuncTbl3 *t)
{
	cnt_ftbl_total++;

	if(s == 0.0){
		// Pure lorentzian
		return lorentzCDF(x, r);
	}

	double xt = (x - c) / s;
	double rt = r / s;

	if(rt > t->rmax){
		// When rt > rmax, it is almost pure lorentzian with HWHM rt
		// Gaussian with sigma=1 behaves as a delta-function
		return lorentzCDF(xt, rt);
	}

	if(xt > t->xmax){
		return 1.0;
	} else if(xt < -(t->xmax)){
		return 0.0;
	}

	int idx, idr;
	double px, pr;
	// getIDandPortionNew(x, t, &idx, &px);
	// getIDandPortion(r, t->rmin, t->rmax, t->dr, t->nr, &idr, &pr);
	// getIDandPortionLog(x, t->xmin, t->xmax, t->dx, t->nx, &idx, &px);
	getIDandPortionLog(fabs(xt), t->xmin, t->xmax, t->dx, t->nx, &idx, &px);
	getIDandPortionLog(rt, t->rmin, t->rmax, t->dr, t->nr, &idr, &pr);

#ifdef DEBUGTBLFUNC
	printf("xt = %f, idx/nx = %d/%d, px = %f\n", xt, idx, t->nx, px);
	printf("rt = %f, idr/nr = %d/%d, pr = %f\n", rt, idr, t->nr, pr);
	printf("xt = %f, x[%d] = %f, x[%d] = %f\n", xt, idx, t->xval[idx], idx+1, t->xval[idx+1]);
	printf("rt = %f, r[%d] = %f, r[%d] = %f\n", rt, idr, t->rval[idr], idr+1, t->rval[idr+1]);
	printf("fval= %e, %e, %e, %e\n", t->fval[idr][idx], t->fval[idr + 1][idx], t->fval[idr][idx + 1], t->fval[idr + 1][idx + 1]);
#endif

	// F(|x/s|, 1, g/s)
	// double ftmp = bilinear_interpolation(px, pr, t->fval[idr][idx], t->fval[idr][idx + 1], t->fval[idr + 1][idx], t->fval[idr + 1][idx + 1]);
	double ftmp = bilinear_interpolation(px, pr, t->fint[idr][idx], t->fint[idr][idx + 1], t->fint[idr + 1][idx], t->fint[idr + 1][idx + 1]);

	if (xt < 0.0) {
		// When xt < 0
		// return 1.0 - ftmp;
		ftmp =  1.0 - ftmp;
	}
	else {
		// When xt >= 0
		// return ftmp;
	}

	// Return CDF of Voigt function whose peak area = 1
	return ftmp;

	// // Return CDF of Voigt function normalized at x=0
	// px = 0.0; idx = 0;
	// double v0 = (1.0/s) * bilinear_interpolation(px, pr, t->fval[idr][idx], t->fval[idr][idx + 1], t->fval[idr + 1][idx], t->fval[idr + 1][idx + 1]);
	// return ftmp v0;
}

// Initialize FuncTbl3* 
FuncTbl3* newFuncTbl()
{
	FuncTbl3 *v;
	v = (FuncTbl3*)malloc(sizeof(FuncTbl3));
	memset(v, 0, sizeof(FuncTbl3));
	return v;
}





#define TBLFUNC_TEST

#ifdef TBLFUNC_TEST

// #define NDATMAX 512          // スペクトルの最大データ数
// #define NXMAX 201         // 2次元関数テーブルの横軸x方向の最大データ数
// #define NRMAX 102         // 2次元関数テーブルのローレンツHWHM方向の最大データ数
#define NXMAX 1001         // 2次元関数テーブルの横軸x方向の最大データ数
#define NRMAX 501         // 2次元関数テーブルのローレンツHWHM方向の最大データ数
#define VOIGT_REF "../make_table/VoigtDB.csv"
//#define VOIGT_REF "./dat/psVoigtSum.csv"
//#define VOIGT_REF "./dat/psVoigtSum_nx10001_nr101.csv"
//#define VOIGT_REF "./dat/psVoigtSum_nx101_nr6.csv"
#define NBUF  1024


#define MAX(a, b) (((a) > (b)) ? (a) : (b))
#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#define SQR(x) ((x)*(x))

// 観測データと積分テーブル (関数内におくと自動領域に取られるので外部に置いとく)
// double X[NDATMAX], Y[NDATMAX];
double VX[NXMAX];
double VR[NRMAX];
//double VGT[NRMAX][NXMAX], IVGT[NRMAX][NXMAX];
double **VGT, **IVGT;



double Gaussian(double x);
double psVoigtSum(double x, double r);
double Voigt(double x, double s, double g);
// void evaluateErrorByGridSearch(FuncTbl3 *Vg);
// void evaluateErrorByLogGridSearch(FuncTbl3 *Vg);
// void evaluateErrorByRandomSearch(FuncTbl3 *Vg);
// void output_for_python(FuncTbl3 *Vg);

// pseudo-Voigt関数のテーブルを読み込む
// void LoadFtbl(const char *reffile, FuncTbl3 *Vg, FuncTbl3 *VgInt)
void LoadFtbl(const char *reffile, FuncTbl3 *Vg)
{
	int j;
	double x, r, v, vint;
	printf("Load from %s ...\n", reffile);

	//動的に確保
	// VGT = malloc(sizeof(double*)*NRMAX);
	// for (j = 0; j < NRMAX; j++) {
	// 	VGT[j] = malloc(sizeof(double)*NXMAX);
	// }
	VGT = new double*[NRMAX];
	for (j = 0; j < NRMAX; j++) {
		VGT[j] = new double[NXMAX];
	}

	// IVGT = malloc(sizeof(double*)*NRMAX);
	// for (j = 0; j < NRMAX; j++) {
	// 	IVGT[j] = malloc(sizeof(double)*NXMAX);
	// }
	IVGT = new double*[NRMAX];
	for (j = 0; j < NRMAX; j++) {
		IVGT[j] = new double[NXMAX];
	}

	//char *tokenPtr;
	char buf[NBUF];
	FILE *fp;

	// VX, VGT, IVGT に Gaussian-Lorenz データ読み込み
	if ((fp = fopen(reffile, "r")) == NULL)    err(EXIT_FAILURE, "Cannot open %s", reffile);
	// 最初の１行はヘッダ
	if (fgets(buf, NBUF, fp) == NULL) err(EXIT_FAILURE, "No data in %s", reffile);

	int nx = -1;
	int nr = -1;
	double r_old = -1.0;

	while (!feof(fp)) {
		if (fgets(buf, NBUF, fp) == NULL) break;

		if (sscanf(buf, "%lf,%lf,%lf,%lf", &x, &r, &v, &vint) == 0) err(EXIT_FAILURE, "Error to read a line:\n%s\n", buf);

		if (r_old != r) {
			nx = -1;
			nr++;
			r_old = r;
		}
		nx++;

		VX[nx] = x;
		VR[nr] = r;
		VGT[nr][nx] = v;
		IVGT[nr][nx] = vint;

		//fprintf(stderr, "%d, %d: %lf, %lf, %lf, %lf\n", nx, nr, VX[nx], VR[nr], VGT[nr][nx], IVGT[nr][nx]);
		//fflush(stderr);
	}
	fclose(fp);

	nx++;
	nr++;

	if (nx > NXMAX) {
		err(EXIT_FAILURE, "Insufficient memory size: nx=%d (max %d)\n", nx, NXMAX);
	}
	if (nr > NRMAX) {
		err(EXIT_FAILURE, "Insufficient memory size: nr=%d (max %d)\n", nr, NRMAX);
	}

	// 取得した値を代入
	Vg->nx = nx;
	// Vg->xmin = VX[0];
	Vg->xmin = VX[1]; // Because VX[0] = 0.0
	Vg->xmax = VX[nx - 1];
	// Vg->dx = (VX[nx - 1] - VX[0]) / (nx - 1); //undefined
	Vg->dx = log(Vg->xmax / Vg->xmin) / (nx - 2);  // log-scale
	Vg->nr = nr;
	// Vg->rmin = VR[0];
	Vg->rmin = VR[1]; // Because VR[0] = 0.0
	// Vg->rmax = VR[nr - 2]; // Because VR[nr-1] = 1.0 is only for pure lorentzian
	Vg->rmax = VR[nr - 1];
	//Vg->dr = (VR[nr - 1] - VR[0]) / (nr - 1);
	// Vg->dr = log(VR[nr - 1] / VR[0]) / (nr - 1);  // log-scale
	// Vg->dr = log(VR[nr - 2] / VR[1]) / (nr - 3);  // log-scale
	Vg->dr = log(Vg->rmax / Vg->rmin) / (nr - 2);  // log-scale
	Vg->fval = VGT;

	//Vg->c = 0.0002;
	Vg->c = 1.0/(double)(Vg->nx - 1); // not used
	Vg->xval = VX;
	Vg->rval = VR;
	//printf("NXMAX = %d\n", NXMAX);
	//printf("NRMAX = %d\n", NRMAX);

	// 関数テーブルのコピーを作ってから積分テーブルを代入
	// *VgInt = *Vg;
	// VgInt->fval = IVGT;

	// Table of CDF
	Vg->fint = IVGT;

	printf("nx = %d, xmin = %f, xmax = %f, dx = %f\n", Vg->nx, Vg->xmin, Vg->xmax, Vg->dx);
	printf("nr = %d, rmin = %f, rmax = %f, dr = %f\n", Vg->nr, Vg->rmin, Vg->rmax, Vg->dr);
	fflush(stderr);
}

#ifdef DEBUGTBLFUNC
int main(void)
{
	//const int cntmax = 100000000;
	// const int cntmax = 10;
	//repeat_calc_by_function(cntmax);
	//exit(1);

	int i, j;
	//double x, r;

	//FuncTbl3 *tbl;
	//tbl = SetupFtbl(psVoigtSum, XMIN, XMAX, NDIVX, RMIN, RMAX, NDIVR);

	FuncTbl3 *Vg = newFuncTbl();     // 疑似Voigt関数のテーブル
	// FuncTbl3 *VgInt = newFuncTbl();  // 疑似Voigt関数の積分テーブル
	
	//LoadTable(Vg, VgInt);
	// LoadFtbl(VOIGT_REF, Vg, VgInt);
	LoadFtbl(VOIGT_REF, Vg);

	// printf("Vg->c = %f\n", Vg->c);
	// printf("nx = %d, xmin = %f, xmax = %f, dx = %f\n", Vg->nx, Vg->xmin, Vg->xmax, Vg->dx);
	// printf("nr = %d, rmin = %f, rmax = %f, dr = %f\n", Vg->nr, Vg->rmin, Vg->rmax, Vg->dr);
	fflush(stderr);

	// Simple view
	for (j = 0; j < Vg->nr; j+=(Vg->nr/5)) {
		for (i = 0; i < Vg->nx; i+=(Vg->nx/10)) {
			// printf("Vg->fval[%d][%d] = %f, VgInt->fval[%d][%d] = %f\n", j, i, Vg->fval[j][i], j, i, VgInt->fval[j][i]);
			printf("Vg->fval[%d][%d] = %f, Vg->fint[%d][%d] = %f\n", j, i, Vg->fval[j][i], j, i, Vg->fint[j][i]);
		}
	}

     double x = -100.0; // x-value
	 double c = 3.0;    // peak position
     double r = 0.001;  // HWHM of lorentzian
	 double s = 2.0;    // STD of gaussian

	//  printf("%g,%g,%g,%g,%g,%g\n", x, c, s, r, ftbl(x, c, s, r, Vg), finttbl(x, c, s, r, VgInt));
	 printf("%g,%g,%g,%g,%g,%g\n", x, c, s, r, ftbl(x, c, s, r, Vg), finttbl(x, c, s, r, Vg));

	//// グリッドサーチで誤差を評価する
	//evaluateErrorByGridSearch(Vg);
	//evaluateErrorByLogGridSearch(Vg);

	////乱数を使って誤差を評価する
	////evaluateErrorByRandomSearch(Vg);

	//// For Python plot
	////output_for_python(Vg);

	//// Time
	//repeat_calc_by_table(Vg, cntmax);
	////repeat_calc_by_function(Vg, cntmax);

	//output_cnt_ftbl();

	//対話型でチェックする
	do {
		int returnscanf = scanf("%lf %lf %lf %lf", &x, &c, &s, &r);
		// printf("%g,%g,%g,%g,%g,%g\n", x, c, s, r, ftbl(x, c, s, r, Vg), finttbl(x, c, s, r, VgInt));
		printf("%g,%g,%g,%g,%g,%g\n", x, c, s, r, ftbl(x, c, s, r, Vg), finttbl(x, c, s, r, Vg));
	} while(1);

	printf("Finished\n");
    return 0;
}
#endif


double Gaussian(double x)
{
    return exp(-0.5*SQR(x-1.0));
}

double psVoigtSum(double x, double r)
{
    // return r / (1.0 + SQR(x-1.0)) + (1.0-r) * pow(2.0, -SQR(x-1.0));
    return r / (1.0 + SQR(x)) + (1.0-r) * pow(2.0, -SQR(x));
}

// double Voigt(double x, double s, double g)
// {

// 	return 0.0;
// }

#endif // TBLFUNC_TEST

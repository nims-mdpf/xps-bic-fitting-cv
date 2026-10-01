#ifndef __TBLFUNC_H__
#define __TBLFUNC_H__ 1

typedef struct _ftbl3 {
    int nx;                 // number of x data
    double xmin, xmax, dx;  // start, end and step size of x data
    int nr;                 // number of r data
    double rmin, rmax, dr;  // start, end and step size of r data
	double c;               // additional parameter
	double *xval;
	double *rval;
    double **fval;          // 2D table of the function
    double **fint;          // 2D table of the CDF
} FuncTbl3;

extern long cnt_ftbl_total;
extern long cnt_ftbl_outside;

void getIDandPortion(const double x, const double xmin, const double xmax, const double dx, const int nx, int *id, double *p);
void getIDandPortionLog(const double x, const double xmin, const double xmax, const double dx, const int nx, int *id, double *p);
void getIDandPortionNew(const double x, const FuncTbl3 *t, int *id, double *p);
double bilinear_interpolation(double px, double py, double f00, double f10, double f01, double f11);
void output_cnt_ftbl();
// double ftbl(double x, double r, FuncTbl3 *t);
// double finttbl(double x, double r, FuncTbl3 *t);

double ftbl(double x, double c, double s, double r);
double ftbl(double x, double c, double s, double r, FuncTbl3 *t);
// double ftbltbl(double x, double c, double s, double r, FuncTbl3 *t);
double finttbl(double x, double c, double s, double r, FuncTbl3 *t);
FuncTbl3* newFuncTbl();
//FuncTbl3* SetupFtbl(double(*f)(double x, double r), double xmin, double xmax, int nx, double rmin, double rmax, int nr);
void LoadFtbl(const char *reffile, FuncTbl3 *Vg);
#endif

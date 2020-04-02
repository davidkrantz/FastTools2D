#ifndef __BESSELLIB_H_
#define __BESSELLIB_H_

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <immintrin.h>
#include <gsl/gsl_sf_expint.h>

#if (defined __AVX__ || defined __SSE4_2__)
#include "math_x86.h"
#endif

typedef double complex[2];

#define pi     3.14159265358979323846264338327950288
#define MAX(a, b) ((a > b) ? a : b)
#define MIN(a, b) ((a > b) ? b : a)
#define EPS    1.110223024625157e-16
#define EULER  0.577215664901532860606512090082402431042

double K0x(double);
double K1x(double);
double E0(double);
double E1(double);
double E2(double);
double K0xy(double, double, int, const int);
double K1xy(double, double, int, const int);
double Km1xy(double, double, int, const int);
int choose_algorithm(double, double, int*);
double Cjnu(double, int, double);
double cjnk(int, int, int);
double Pochhammer(int, int);
double factorial(int);
double Uj1x(int, double);
double Uj0x(int, double);
double Uj2x(int, double);


#endif

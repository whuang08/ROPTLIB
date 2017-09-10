#ifndef WAVELET_H
#define WAVELET_H

#include <cmath>
#include <iostream>
#ifndef MATLAB_MEX_FILE
#include <f2c.h>
#else
#include "blas.h"
#include "lapack.h"
#define integer ptrdiff_t
#endif

void haarFWT_1d(int n, doublecomplex *v);
void haarFWT_1d_inverse(int n, doublecomplex *v);

void haarFWT_2d(int n1, int n2, doublecomplex *v);
void haarFWT_2d_inverse(int n, int n2, doublecomplex *v);

#endif/* WAVELET_H */

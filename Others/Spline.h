/*
This file defines function for computing cubic spline curve for various functions.
TODO: Define the class to be a spline curve, not just a class contains many static functions.

-----WH
*/

#ifndef SPLINE_H
#define SPLINE_H

#include <iostream>
#include <limits>
#include "Others/def.h"

/*Define the namespace*/
namespace ROPTLIB{

	class Spline{
	public:
		static int SplineUniformPeriodic(const realdp *Y, int n, realdp h, realdp *coefs); /* periodic boundary condition: for closed curves */
		static int SplinePeriodic(const realdp *X, const realdp *Y, int n, realdp *coefs);

		static int SplineUniformSlopes(const realdp *Y, int n, realdp h, realdp *coefs); /* Hermite boundary conditions: for open curves (derivative-free from) */
		static int SplineSlopes(const realdp *X, const realdp *Y, int n, realdp *coefs);

		static int SolveTridiagonalSystem(realdp *d, realdp *ud, realdp *ld, realdp *vec, realdp *s, int n);
		static int SolvePeriodicSystem(realdp *d, realdp *ud, realdp *ld, realdp *vec, realdp *s, int nn);

		static realdp ValSpline(const realdp *coefs, const realdp *breaks, int N, realdp t);
		static realdp ValSplineUniform(const realdp *coefs, int N, realdp h, realdp t);
		static void FirstDeri(const realdp *coefs, int N, realdp *dericoefs);
		static realdp ValFirstDeriUniform(const realdp *dericoefs, int N, realdp h, realdp t);
		static realdp ValFirstDeri(const realdp *dericoefs, const realdp *breaks, int N, realdp t);
		static void SecondDeri(const realdp *coefs, int N, realdp *dericoefs);
		static realdp ValSecondDeriUniform(const realdp *dericoefs, int N, realdp h, realdp t);
		static realdp ValSecondDeri(const realdp *dericoefs, const realdp *breaks, int N, realdp t);
	};
}; /*end of ROPTLIB namespace*/
#endif /* SPLINE_H */

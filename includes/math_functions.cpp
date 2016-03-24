/*
 * math_functions.cpp
 *
 *  Created on: Mar 22, 2016
 *      Author: kiliakis
 */

#include "math_functions.h"

namespace mymath {
// Function to implement integration of f(x) over the interval
// [a,b] using the trapezoid rule with nsub subdivisions.
ftype *cum_trapezoid(ftype *f, ftype deltaX, int nsub) {
	// initialize the partial sum to be f(a)+f(b) and
	// deltaX to be the step size using nsub subdivisions
	ftype *psum = new ftype[nsub];
	psum[0] = 0;

	// increment the partial sum
	for (int i = 1; i < nsub; i++)
		psum[i] = psum[i - 1] + (f[i] + f[i - 1]) * (deltaX / 2.0);

	return psum;

}

ftype trapezoid(ftype *f, ftype deltaX, int nsub) {
	// initialize the partial sum to be f(a)+f(b) and
	// deltaX to be the step size using nsub subdivisions
	ftype psum = f[0] + f[nsub - 1]; //f(a)+f(b);
	//ftype deltaX = (b-a)/nsub;

	// increment the partial sum
	for (int index = 1; index < nsub; index++) {
		psum = psum + 2.0 * f[index];
	}

	// multiply the sum by the constant deltaX/2.0
	psum = (deltaX / 2.0) * psum;

	// return approximation
	return psum;

}

int min(ftype * a, int size, int step) {
	int p = 0;
	ftype min = a[0];
	for (int i = 0; i < size; i += step) {
		if (a[i] < min) {
			min = a[i];
			p = i;
		}
	}
	return p;

}

int max(ftype * a, int size, int step) {
	int p = 0;
	ftype max = a[0];
	for (int i = 0; i < size; i += step) {
		if (a[i] > max) {
			max = a[i];
			p = i;
		}
	}
	return p;

}

void linspace(ftype* a, const ftype start, const ftype end, const int n) {
	ftype step = (end - start) / (n - 1);
	ftype value = start;
	for (int i = 0; i < n; ++i) {
		a[i] = value;
		value += step;
	}
}

ftype mean(const ftype data[], const int n) {
	ftype m = 0;
	for (int i = 0; i < n; ++i) {
		m += data[i];
	}
	return m / n;
}

ftype standard_deviation(const ftype data[], const int n, const ftype mean) {
	ftype sum_deviation = 0.0;
	int i;
	for (i = 0; i < n; ++i)
		sum_deviation += (data[i] - mean) * (data[i] - mean);
	return sqrt(sum_deviation / n);
}

}


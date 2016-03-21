/*
 * math_functions.h
 *
 *  Created on: Mar 21, 2016
 *      Author: kiliakis
 */

#ifndef INCLUDES_MATH_FUNCTIONS_H_
#define INCLUDES_MATH_FUNCTIONS_H_

#include <cmath>
#include "configuration.h"


namespace mymath {
// Function to implement integration of f(x) over the interval
// [a,b] using the trapezoid rule with nsub subdivisions.
ftype *trapezoid(ftype *f, ftype deltaX, int nsub) {
	// initialize the partial sum to be f(a)+f(b) and
	// deltaX to be the step size using nsub subdivisions
	ftype *psum = new ftype[nsub];
	//for (int i = 0; i < nsub; ++i)
	psum[0] = 0;

	//psum[0] = f[0] + f[nsub-1];
	//deltaX = 		
	//psum[0] = f(a) + f(b);
	//ftype deltaX = (b - a) / nsub;

	// increment the partial sum
	for (int i = 1; i < nsub; i++) 
		psum[i] = psum[i-1] + (f[i]+f[i-1])*(deltaX / 2.0);

	// multiply the sum by the constant deltaX/2.0
	/*
	for (int i = 0; i < nsub; ++i)
	{
		psum[i] = (deltaX / 2.0) * psum[i];	
	}
	*/
	// return approximation
	return psum;

}

template<typename T>
inline int min(T * a, int size, int step=1)
{
	int p = 0;
	T min = a[0];
	for (int i = 0; i < size; i+=step)
	{
		if(a[i] < min){
			min = a[i];
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


inline ftype mean(const ftype data[], const int n) {
	ftype m = 0;
	for (int i = 0; i < n; ++i) {
		m += data[i];
	}
	return m / n;
}

inline ftype standard_deviation(const ftype data[], const int n, const ftype mean) {
	ftype sum_deviation = 0.0;
	int i;
	for (i = 0; i < n; ++i)
		sum_deviation += (data[i] - mean) * (data[i] - mean);
	return sqrt(sum_deviation / n);
}


}



#endif /* INCLUDES_MATH_FUNCTIONS_H_ */

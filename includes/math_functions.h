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

ftype *cum_trapezoid(ftype *f, ftype deltaX, int nsub);

ftype trapezoid(ftype *f, ftype deltaX, int nsub);

int min(ftype * a, int size, int step = 1);

int max(ftype * a, int size, int step = 1);

void linspace(ftype* a, const ftype start, const ftype end, const int n);

ftype mean(const ftype data[], const int n);

ftype standard_deviation(const ftype data[], const int n, const ftype mean);

}

#endif /* INCLUDES_MATH_FUNCTIONS_H_ */

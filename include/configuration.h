/*
 * configuration.h
 *
 *  Created on: Mar 8, 2016
 *      Author: kiliakis
 */

#ifndef INCLUDES_CONFIGURATION_H_
#define INCLUDES_CONFIGURATION_H_

#include <complex>

typedef double ftype;
typedef unsigned int uint;
typedef std::complex<ftype> complex_t;
const unsigned PRECISION = 8;
// Compile time switches
// Needs re-compile every time a change is made

//switch on/off fixed beam distribution
#define FIXED_PARTICLES

// switch on/off timing
#define TIMING

// switch on/off result printing
#define PRINT_RESULTS

#endif /* INCLUDES_CONFIGURATION_H_ */

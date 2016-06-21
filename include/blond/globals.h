/*
 * globals.h
 *
 *  Created on: Mar 10, 2016
 *      Author: kiliakis
 */

#ifndef INCLUDES_GLOBALS_H_
#define INCLUDES_GLOBALS_H_

#include <blond/beams/Beams.h>
#include <blond/beams/Slices.h>
#include <blond/input_parameters/GeneralParameters.h>
#include <blond/input_parameters/RfParameters.h>

extern GeneralParameters* GP;
extern Beams* Beam;
extern RfParameters* RfP;
extern Slices* Slice;

// TODO num of threads is global
// should it?
extern int n_threads;

#endif /* INCLUDES_GLOBALS_H_ */

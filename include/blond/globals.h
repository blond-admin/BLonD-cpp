/*
 * globals.h
 *
 *  Created on: Mar 10, 2016
 *      Author: kiliakis
 */

#ifndef INCLUDES_GLOBALS_H_
#define INCLUDES_GLOBALS_H_
#include <blond/utilities.h>
#include <blond/input_parameters/GeneralParameters.h>
#include <blond/beams/Beams.h>
#include <blond/input_parameters/RfParameters.h>
#include <blond/beams/Slices.h>

class Beams;

API extern GeneralParameters *GP;
API extern Beams *Beam;
API extern RfParameters *RfP;
API extern Slices *Slice;

// TODO num of threads is global
// should it?
API extern int n_threads;

#endif /* INCLUDES_GLOBALS_H_ */

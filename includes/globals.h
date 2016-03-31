/*
 * globals.h
 *
 *  Created on: Mar 10, 2016
 *      Author: kiliakis
 */

#ifndef INCLUDES_GLOBALS_H_
#define INCLUDES_GLOBALS_H_

#include "../input_parameters/GeneralParameters.h"
#include "../beams/Beams.h"
#include "../input_parameters/RfParameters.h"
#include "../beams/Slices.h"

extern GeneralParameters *GP;
extern Beams *Beam;
extern RfParameters *RfP;
extern Slices *Slice;

// TODO num of threads is global
// should it?
extern int n_threads;

#endif /* INCLUDES_GLOBALS_H_ */

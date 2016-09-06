/*
 * utilities.h
 *
 *  Created on: Sep 5, 2016
 *      Author: kiliakis
 */

#ifndef TRACKERS_UTILITIES_H_
#define TRACKERS_UTILITIES_H_

#include <blond/configuration.h>
#include <blond/utilities.h>
#include <blond/input_parameters/GeneralParameters.h>
#include <blond/input_parameters/RfParameters.h>

ftype phase_modulo_below_transition(const ftype phi);


ftype phase_modulo_above_transition(const ftype phi);


std::vector<bool> is_in_separatrix(GeneralParameters *GP, RfParameters *RfP,
                                   Beams *Beam, f_vector_t &dt, f_vector_t &dE,
                                   f_vector_t total_voltage = {});

f_vector_t hamiltonian(GeneralParameters *GP, RfParameters *RfP, Beams *Beam,
                  f_vector_t &dt, f_vector_t &dE,
                  f_vector_t total_voltage = {});

ftype hamiltonian(GeneralParameters *GP, RfParameters *RfP, Beams *Beam,
                  ftype dt, ftype dE,
                  f_vector_t total_voltage = {});

#endif /* TRACKERS_UTILITIES_H_ */

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
#include <blond/constants.h>


static inline double phase_modulo_above_transition(const double phi)
{
    return phi - 2.0 * constant::pi * std::floor(phi / (2.0 * constant::pi));
}


static inline double phase_modulo_below_transition(const double phi)
{
    return phi - 2.0 * constant::pi *
           (std::floor(phi / (2.0 * constant::pi) + 0.5));
}

std::vector<int> is_in_separatrix(const GeneralParameters *GP,
                                   const RfParameters *RfP,
                                   const Beams *Beam,
                                   const f_vector_t &dt,
                                   const f_vector_t &dE,
                                   const f_vector_t total_voltage = {});

f_vector_t hamiltonian(const GeneralParameters *GP,
                       const RfParameters *RfP,
                       const Beams *Beam,
                       const double *__restrict dt,
                       const double *__restrict dE,
                       const int size,
                       const f_vector_t total_voltage = {});

double hamiltonian(const GeneralParameters *GP,
                  const RfParameters *RfP,
                  const Beams *Beam,
                  const double dt,
                  const double dE,
                  const f_vector_t total_voltage = {});

#endif /* TRACKERS_UTILITIES_H_ */

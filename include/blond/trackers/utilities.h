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

void minmax_location(const f_vector_t &x,const f_vector_t &f,
                     f_vector_t &min_x_position, f_vector_t &max_x_position,
                     f_vector_t &min_values, f_vector_t &max_values);

void potential_well_cut(const f_vector_t &theta_coord_array,
                        const f_vector_t &potential_array,
                        f_vector_t &theta_coord_sep,
                        f_vector_t &potential_well_sep);

#endif /* TRACKERS_UTILITIES_H_ */

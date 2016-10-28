/*
 * Distributions.h
 *
 *  Created on: Mar 16, 2016
 *      Author: kiliakis
 */

#ifndef BEAMS_DISTRIBUTIONS_H_
#define BEAMS_DISTRIBUTIONS_H_

#include <blond/configuration.h>
// #include <blond/constants.h>
// #include <blond/utilities.h>
#include <blond/python.h>
// #include <blond/globals.h>
// #include <blond/math_functions.h>
#include <blond/beams/Beams.h>
#include <blond/trackers/Tracker.h>
#include <blond/impedances/InducedVoltage.h>
#include <map>



void matched_from_line_density(Beams *beam,
                               FullRingAndRf *full_ring,
                               std::map<std::string, std::string> line_density_opt,
                               std::string main_harmonic = "lowest_freq",
                               TotalInducedVoltage *totVolt = nullptr,
                               std::string plot = "",
                               std::string figdir = "fig",
                               std::string half_option = "first",
                               std::map<std::string, std::string> extraVoltageDict =
                                   std::map<std::string, std::string>(),
                               int n_iterations_input = 100,
                               int seed = 0);


void matched_from_distribution_density(FullRingAndRf *full_ring,
                                       std::map<std::string, std::string> distribution_opt,
                                       std::string main_harmonic = "lowest_freq",
                                       int n_iterations_input = 1,
                                       std::map<std::string, std::string> extraVoltageDict =
                                               std::map<std::string, std::string>(),
                                       int seed = 0);


void longitudinal_bigaussian(GeneralParameters *GP, RfParameters *RfP,
                             Beams *Beam, double sigma_dt, double sigma_dE = 0,
                             int seed = 0, bool reinsertion = false);


f_vector_t distribution_density_function(const f_vector_t &action_array,
        const std::string &dist_type, double length, double exponent = 0.);

f_vector_t line_density_function(const f_vector_t &coord_array,
                                 const std::string &dist_type,
                                 const double bunch_length,
                                 const double bunch_position = 0.,
                                 double exponent = 0.);



#endif /* BEAMS_DISTRIBUTIONS_H_ */

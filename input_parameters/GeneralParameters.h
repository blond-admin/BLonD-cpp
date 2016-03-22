/*
 * general_parameters.h
 *
 *  Created on: Mar 8, 2016
 *      Author: kiliakis
 */

#ifndef INPUT_PARAMETERS_GENERALPARAMETERS_H_
#define INPUT_PARAMETERS_GENERALPARAMETERS_H_

//#include "../includes/globals.h"
#include <vector>
#include <cmath>
#include <numeric>
#include <cstring>
#include "../includes/constants.h"
#include "../includes/configuration.h"
#include "../includes/utilities.h"

enum particle_type {
	proton, electron, user_input, none
};

class GeneralParameters {

private:
	void eta_generation();
	void _eta0();
	void _eta1();
	void _eta2();
public:
	int n_sections;
	particle_type particle, particle_2;
	int n_turns;
	ftype mass, mass2;
	ftype charge, charge2;
	ftype cumulative_times;
	ftype *alpha;
	ftype *momentum;
	int alpha_order;
	ftype* ring_length;
	ftype ring_circumference;
	ftype ring_radius;
	ftype *beta;
	ftype *gamma;
	ftype *energy;
	ftype *kin_energy;
	ftype* cycle_time;
	ftype* f_rev, *omega_rev;
	std::vector<ftype> t_rev;
	ftype *eta_0, *eta_1, *eta_2;

	GeneralParameters(const int n_turns, ftype* ring_length, ftype *alpha,
			const int alpha_order, ftype *momentum,
			const particle_type particle, ftype user_mass = 0,
			ftype user_charge = 0, particle_type particle2 = none,
			ftype user_mass_2 = 0, ftype user_charge_2 = 0,
			int number_of_sectrions = 1);

	~GeneralParameters();

};

#endif /* INPUT_PARAMETERS_GENERALPARAMETERS_H_ */

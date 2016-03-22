/*
 * RfParameters.h
 *
 *  Created on: Mar 9, 2016
 *      Author: kiliakis
 */

#ifndef INPUT_PARAMETERS_RFPARAMETERS_H_
#define INPUT_PARAMETERS_RFPARAMETERS_H_

//#include "../includes/globals.h"
#include "GeneralParameters.h"
#include "../beams/Beams.h"
//#include "../includes/utilities.h"
#include "../includes/math_functions.h"
#include "../trackers/sin.h"
#include <algorithm>    // std::cops
#include <iterator>

//#include "../includes/globals.h"

enum accelerating_systems_type {
	as_single, all, first
};

class RfParameters {
public:
	RfParameters(GeneralParameters *gp, Beams *beam, int _n_rf,
			ftype *_harmonic, ftype *_voltage, ftype *_phi_offset,
			ftype* _phi_noise = NULL, ftype * _omega_rf = NULL,
			int _section_index = 1,
			accelerating_systems_type accelerating_systems = as_single);

	ftype *E_increment;
	ftype *phi_s;
	ftype *Qs;
	ftype *omega_s0;
	ftype *omega_RF_d;
	ftype *phi_RF;
	ftype *dphi_RF;
	ftype *dphi_RF_steering;
	ftype *t_RF;
	ftype *omega_RF;

	ftype eta_tracking(const Beams *beam, const int counter, const ftype dE);
	ftype eta_0(const int i);
	ftype eta_1(const int i);
	ftype eta_2(const int i);
	ftype beta(const int i);
	ftype gamma(const int i);
	ftype energy(const int i);
	ftype momentum(const int i);
	int sign_eta_0(const int i);
	// TODO why we pass rf_params??
	void calc_phi_s(const accelerating_systems_type acc_sys = as_single);

	// TODO assume input_value is an array
	// that is why we don't have any input_check function
	int counter;
	GeneralParameters *gp;
	int n_rf;
	//accelerating_systems_type accelerating_systems;
	//int n_turns;
	ftype *harmonic;
	ftype *voltage;
	ftype *phi_offset;
	ftype *phi_noise;
	int section_index;
	ftype length_ratio;
	ftype section_length;

	~RfParameters();

private:
};

#endif /* INPUT_PARAMETERS_RFPARAMETERS_H_ */

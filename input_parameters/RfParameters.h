/*
 * RfParameters.h
 *
 *  Created on: Mar 9, 2016
 *      Author: kiliakis
 */

#include "GeneralParameters.h"
#include "../beams/Beams.h"
#include "../includes/utilities.h"
#include <algorithm>    // std::cops
#include <iterator>

#ifndef INPUT_PARAMETERS_RFPARAMETERS_H_
#define INPUT_PARAMETERS_RFPARAMETERS_H_

//#include "../includes/globals.h"

class RfParameters {
public:
	RfParameters(GeneralParameters *gp, Beams *beam, int _n_rf,
			ftype *_harmonic, ftype *_voltage, ftype *_phi_offset,
			ftype* _phi_noise = NULL, ftype * _omega_rf = NULL,
			int _section_index = 1);
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

	inline ftype eta_tracking(const Beams *beam, const int counter,
			const ftype dE);
	inline ftype eta_0(const int i);
	inline ftype eta_1(const int i);
	inline ftype eta_2(const int i);
	inline ftype beta(const int i);
	inline ftype gamma(const int i);
	inline ftype energy(const int i);
	inline ftype momentum(const int i);
	int sign_eta_0(const int i);
	// TODO assume input_value is an array
	// that is why we don't have any input_check function
	int counter;
	GeneralParameters *gp;
	int n_rf;
	//int n_turns;
	ftype *harmonic;
	ftype *voltage;
	ftype *phi_offset;
	ftype *phi_noise;
	int section_index;
	ftype length_ratio;
	ftype section_length;

private:
};

/*
 :How to use RF programs:

 - For 1 RF system and constant values of V, h or phi, just input the single value
 - For 1 RF system and varying values of V, h or phi, input an array of n_turns values
 - For several RF systems and constant values of V, h or phi, input lists of single values
 - For several RF systems and varying values of V, h or phi, input lists of arrays of n_turns values
 */
// RfParameters == RfSectionParameters
// completely removed accelerating_systems
RfParameters::RfParameters(GeneralParameters *_gp, Beams *beam, int _n_rf,
		ftype *_harmonic, ftype *_voltage, ftype *_phi_offset,
		ftype* _phi_noise, ftype * _omega_RF, int _section_index) {
	this->counter = 0;
	this->gp = _gp;
	this->section_index = _section_index - 1;
	this->n_rf = _n_rf;
	//this->n_turns = gp->n_turns;
	this->harmonic = _harmonic;
	this->voltage = _voltage;
	this->phi_offset = _phi_offset;
	//this->phi_noise = _phi_noise;
	//this->omega_RF = _omega_RF;
	// TODO how to initialize this phi_s
	this->phi_s = new ftype[gp->n_turns];
	for (int i = 0; i < gp->n_turns; ++i)
		phi_s[i] = 0;

	this->section_length = gp->ring_length[section_index];
	this->length_ratio = section_length / gp->ring_circumference;

	// wiped out all the imports
	this->E_increment = new ftype[n_rf * (gp->n_turns)];
	// Don't have n_turns +1 cause of the way np.diff works
	int k = 0;
	for (int i = 0; i < n_rf; ++i) {
		for (int j = 0; j < gp->n_turns; ++j) {
			E_increment[k] = gp->energy[k + 1] - gp->energy[k];
			k++;
		}
		k++;
	}

	this->Qs = new ftype[(gp->n_turns + 1)];
	for (int i = 0; i < gp->n_turns + 1; ++i) {
		Qs[i] = sqrt(
				harmonic[i] * gp->charge * voltage[i]
						* abs(eta_0(i) * cos(phi_s[i]))
						/ (2 * pi * gp->beta[i] * gp->beta[i] * gp->energy[i]));
	}

	this->omega_s0 = new ftype[(gp->n_turns + 1)];
	for (int i = 0; i < (gp->n_turns + 1); ++i) {
		this->omega_s0[i] = Qs[i] * gp->omega_rev[i];
	}

	this->omega_RF_d = new ftype[n_rf * (gp->n_turns + 1)];
	for (int i = 0; i < n_rf * (gp->n_turns + 1); ++i) {
		this->omega_RF_d[i] = 2 * pi * gp->beta[i] * c * harmonic[i]
				/ gp->ring_circumference;
	}
	/*
	 this->omega_RF = new ftype[n_rf * (gp->n_turns + 1)];
	 if (omega_RF == NULL) {
	 dprintf("It was null!\n");
	 for (int i = 0; i < n_rf * (gp->n_turns + 1); ++i) {
	 omega_RF[i] = omega_RF_d[i];
	 }
	 }
	 */
	if (_omega_RF == NULL) {
		this->omega_RF = new ftype[n_rf * (gp->n_turns + 1)];
		std::copy(&omega_RF_d[0], &omega_RF_d[n_rf * (gp->n_turns + 1)],
				omega_RF);
	} else {
		this->omega_RF = _omega_RF;
	}

	this->phi_RF = new ftype[n_rf * (gp->n_turns + 1)];
	//this->phi_RF = (ftype *) aligned_malloc(
	//		sizeof(ftype) * n_rf * (gp->n_turns + 1));
	std::copy(&phi_offset[0], &phi_offset[n_rf * (gp->n_turns + 1)], phi_RF);

	this->dphi_RF = new ftype[n_rf];
	std::fill_n(dphi_RF, n_rf, 0);

	this->dphi_RF_steering = new ftype[n_rf];
	std::fill_n(dphi_RF_steering, n_rf, 0);

	this->t_RF = new ftype[gp->n_turns + 1];

	for (int i = 0; i < gp->n_turns + 1; ++i) {
		t_RF[i] = 2 * pi / omega_RF[i];
	}
}

inline ftype RfParameters::eta_tracking(const Beams *beam, const int counter,
		const ftype dE) {
	ftype eta = 0;
	if (gp->alpha_order == 1)
		eta = eta_0(counter);
	else {
		ftype delta = dE
				/ ((beam->gp->beta[0]) * (beam->gp->beta[0])
						* beam->gp->energy[0]);
		eta += eta_0(counter) * 1;
		if (gp->alpha_order > 0)
			eta += eta_1(counter) * delta;
		if (gp->alpha_order > 1)
			eta += eta_2(counter) * delta * delta;
		if (gp->alpha_order > 2)
			dprintf(
					"WARNING: Momentum compaction factor is implemented only up to 2nd order");
	}
	return eta;

}

inline ftype RfParameters::eta_0(const int i) {
	return gp->eta_0[section_index * (gp->n_turns + 1) + i];
}

inline ftype RfParameters::eta_1(const int i) {
	return gp->eta_1[section_index * (gp->n_turns + 1) + i];
}

inline ftype RfParameters::eta_2(const int i) {
	return gp->eta_2[section_index * (gp->n_turns + 1) + i];
}

inline ftype RfParameters::beta(const int i) {
	return gp->beta[section_index * (gp->n_turns + 1) + i];
}

inline ftype RfParameters::gamma(const int i) {
	return gp->gamma[section_index * (gp->n_turns + 1) + i];

}

inline ftype RfParameters::energy(const int i) {
	return gp->energy[section_index * (gp->n_turns + 1) + i];

}

inline ftype RfParameters::momentum(const int i) {
	return gp->momentum[section_index * (gp->n_turns + 1) + i];

}

inline int RfParameters::sign_eta_0(const int i) {
	if (gp->eta_0[section_index * (gp->n_turns + 1) + i] > 0)
		return 1;
	else if (gp->eta_0[section_index * (gp->n_turns + 1) + i] == 0)
		return 0;
	else
		return -1;
}

#endif /* INPUT_PARAMETERS_RFPARAMETERS_H_ */

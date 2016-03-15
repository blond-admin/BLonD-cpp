/*
 * Tracker.h
 *
 *  Created on: Mar 10, 2016
 *      Author: kiliakis
 */

#include "../input_parameters/GeneralParameters.h"
#include "../beams/Beams.h"
#include "../includes/utilities.h"
#include "sin.h"

#ifndef TRACKERS_TRACKER_H_
#define TRACKERS_TRACKER_H_

enum solver_type {
	simple, full
};

// TODO
// !!!!WARNING!!!!
// we only use beam->dt, dE
// bool arrays are used every time to update only the right values!!

class RingAndRfSection {

private:

public:
	bool *indices_right_outside;
	bool *indices_inside_frame;
	bool *indices_left_outside;
	void set_periodicity();
	inline void kick(const bool *update, const int index);
	inline void kick(const ftype * __restrict__ beam_dt,
			ftype * __restrict__ beam_dE, const int n_rf,
			const ftype * __restrict__ voltage,
			const ftype * __restrict__ omega_RF,
			const ftype * __restrict__ phi_RF, const int n_macroparticles,
			const ftype acc_kick, const bool * __restrict__ update);
	inline void kick(const int index);
	inline void kick(const ftype * __restrict__ beam_dt,
			ftype * __restrict__ beam_dE, const int n_rf,
			const ftype * __restrict__ voltage,
			const ftype * __restrict__ omega_RF,
			const ftype * __restrict__ phi_RF, const int n_macroparticles,
			const ftype acc_kick);
	inline void drift(const bool *update, const int index);
	inline void drift(ftype * __restrict__ beam_dt,
			const ftype * __restrict__ beam_dE, const solver_type solver,
			const ftype T0, const ftype length_ratio, const int alpha_order,
			const ftype eta_zero, const ftype eta_one, const ftype eta_two,
			const ftype beta, const ftype energy, const int n_macroparticles,
			const bool * __restrict__ update);
	inline void drift(const int index);
	inline void drift(ftype * __restrict__ beam_dt,
			const ftype * __restrict__ beam_dE, const solver_type solver,
			const ftype T0, const ftype length_ratio, const int alpha_order,
			const ftype eta_zero, const ftype eta_one, const ftype eta_two,
			const ftype beta, const ftype energy, const int n_macroparticles);
	void track();
	void horizontal_cut();
	RingAndRfSection(GeneralParameters *gp, RfParameters *rf_params,
			Beams *beam, solver_type solver = simple, ftype *PhaseLoop = NULL,
			ftype * NoiseFB = NULL, bool periodicity = false, ftype dE_max = 0,
			bool rf_kick_interp = false, ftype* Slices = NULL,
			ftype * TotalInducedVoltage = NULL, int n_threads = 1);
	RfParameters *rf_params;
	GeneralParameters *gp;
	Beams *beam;
	solver_type solver;
	ftype *PhaseLoop;
	ftype * NoiseFB;
	bool periodicity;
	ftype dE_max;
	bool rf_kick_interp;
	ftype* Slices;
	ftype * TotalInducedVoltage;
	int n_threads;
	ftype* acceleration_kick;
};

// Two versions of kick, drift one with periodicity and another without periodiciy
inline void RingAndRfSection::kick(const int index) {
	kick(beam->dt, beam->dE, rf_params->n_rf, &rf_params->voltage[index],
			&rf_params->omega_RF[index], &rf_params->phi_RF[index],
			beam->n_macroparticles, acceleration_kick[index]);
}

inline void RingAndRfSection::kick(const ftype * __restrict__ beam_dt,
		ftype * __restrict__ beam_dE, const int n_rf,
		const ftype * __restrict__ voltage, const ftype * __restrict__ omega_RF,
		const ftype * __restrict__ phi_RF, const int n_macroparticles,
		const ftype acc_kick) {

// KICK
	// TODO try to remove this if here
	int k = 0;

	for (int j = 0; j < n_rf; j++) {
#pragma omp parallel for
		for (int i = 0; i < n_macroparticles; i++)
			beam_dE[i] += voltage[k]
					* vdt::fast_sin(omega_RF[k] * beam_dt[i] + phi_RF[k]);
		k += n_macroparticles;
	}

// SYNCHRONOUS ENERGY CHANGE
#pragma omp parallel for
	for (int i = 0; i < n_macroparticles; i++)
		beam_dE[i] += acc_kick;

}

inline void RingAndRfSection::kick(const bool* update, const int index) {
	kick(beam->dt, beam->dE, rf_params->n_rf, &rf_params->voltage[index],
			&rf_params->omega_RF[index], &rf_params->phi_RF[index],
			beam->n_macroparticles, acceleration_kick[index], update);
}

inline void RingAndRfSection::kick(const ftype * __restrict__ beam_dt,
		ftype * __restrict__ beam_dE, const int n_rf,
		const ftype * __restrict__ voltage, const ftype * __restrict__ omega_RF,
		const ftype * __restrict__ phi_RF, const int n_macroparticles,
		const ftype acc_kick, const bool * __restrict__ update) {

// KICK
	// TODO try to remove this if here
	int k = 0;
	for (int j = 0; j < n_rf; j++) {
		for (int i = 0; i < n_macroparticles; i++)
			beam_dE[i] +=
					update[i] ?
							voltage[k]
									* vdt::fast_sin(
											omega_RF[k] * beam_dt[i]
													+ phi_RF[k]) :
							0;
		k += n_macroparticles;
	}

// SYNCHRONOUS ENERGY CHANGE
	for (int i = 0; i < n_macroparticles; i++)
		beam_dE[i] += update[i] ? acc_kick : 0;

}

inline void RingAndRfSection::drift(ftype * __restrict__ beam_dt,
		const ftype * __restrict__ beam_dE, const solver_type solver,
		const ftype T0, const ftype length_ratio, const int alpha_order,
		const ftype eta_zero, const ftype eta_one, const ftype eta_two,
		const ftype beta, const ftype energy, const int n_macroparticles,
		const bool * __restrict__ update) {

	int i;
	ftype T = T0 * length_ratio;

	if (solver == simple) {
		ftype coeff = eta_zero / (beta * beta * energy);

		for (i = 0; i < n_macroparticles; i++)
			beam_dt[i] += update[i] ? T * coeff * beam_dE[i] : 0;
	}

	else {
		const ftype coeff = 1. / (beta * beta * energy);
		const ftype eta0 = eta_zero * coeff;
		const ftype eta1 = eta_one * coeff * coeff;
		const ftype eta2 = eta_two * coeff * coeff * coeff;

		if (alpha_order == 1)
			for (i = 0; i < n_macroparticles; i++)
				beam_dt[i] +=
						update[i] ?
								T * (1. / (1. - eta0 * beam_dE[i]) - 1.) : 0;
		else if (alpha_order == 2)
			for (i = 0; i < n_macroparticles; i++)
				beam_dt[i] +=
						update[i] ?
								T
										* (1.
												/ (1. - eta0 * beam_dE[i]
														- eta1 * beam_dE[i]
																* beam_dE[i])
												- 1.) :
								0;
		else
			for (i = 0; i < n_macroparticles; i++)
				beam_dt[i] +=
						update[i] ?
								T
										* (1.
												/ (1. - eta0 * beam_dE[i]
														- eta1 * beam_dE[i]
																* beam_dE[i]
														- eta2 * beam_dE[i]
																* beam_dE[i]
																* beam_dE[i])
												- 1.) :
								0;
	}

}

inline void RingAndRfSection::drift(const bool *update, const int index) {
	drift(beam->dt, beam->dE, solver, gp->t_rev[index], rf_params->length_ratio,
			gp->alpha_order, rf_params->eta_0(index), rf_params->eta_1(index),
			rf_params->eta_2(index), rf_params->beta(index),
			rf_params->energy(index), beam->n_macroparticles, update);
}

inline void RingAndRfSection::drift(ftype * __restrict__ beam_dt,
		const ftype * __restrict__ beam_dE, const solver_type solver,
		const ftype T0, const ftype length_ratio, const int alpha_order,
		const ftype eta_zero, const ftype eta_one, const ftype eta_two,
		const ftype beta, const ftype energy, const int n_macroparticles) {

	int i;
	ftype T = T0 * length_ratio;

	if (solver == simple) {
		ftype coeff = eta_zero / (beta * beta * energy);
#pragma omp parallel for
		for (i = 0; i < n_macroparticles; i++)
			beam_dt[i] += T * coeff * beam_dE[i];
	}

	else {
		const ftype coeff = 1. / (beta * beta * energy);
		const ftype eta0 = eta_zero * coeff;
		const ftype eta1 = eta_one * coeff * coeff;
		const ftype eta2 = eta_two * coeff * coeff * coeff;

		if (alpha_order == 1)
			for (i = 0; i < n_macroparticles; i++)
				beam_dt[i] += T * (1. / (1. - eta0 * beam_dE[i]) - 1.);
		else if (alpha_order == 2)
			for (i = 0; i < n_macroparticles; i++)
				beam_dt[i] += T
						* (1.
								/ (1. - eta0 * beam_dE[i]
										- eta1 * beam_dE[i] * beam_dE[i]) - 1.);
		else
			for (i = 0; i < n_macroparticles; i++)
				beam_dt[i] += T
						* (1.
								/ (1. - eta0 * beam_dE[i]
										- eta1 * beam_dE[i] * beam_dE[i]
										- eta2 * beam_dE[i] * beam_dE[i]
												* beam_dE[i]) - 1.);
	}

}

inline void RingAndRfSection::drift(const int index) {
	drift(beam->dt, beam->dE, solver, gp->t_rev[index], rf_params->length_ratio,
			gp->alpha_order, rf_params->eta_0(index), rf_params->eta_1(index),
			rf_params->eta_2(index), rf_params->beta(index),
			rf_params->energy(index), beam->n_macroparticles);
}

inline void RingAndRfSection::track() {
	//omp_set_num_threads(n_threads);
	if (periodicity) {
		// Change reference of all the particles on the right of the current
		// frame; these particles skip one kick and drift
		for (int i = 0; i < beam->n_macroparticles; ++i) {
			beam->dt[i] -= indices_right_outside[i]
					* gp->t_rev[rf_params->counter + 1];
		}
		// Synchronize the bunch with the particles that are on the right of
		// the current frame applying kick and drift to the bunch; after that
		// all the particle are in the new updated frame

		kick(indices_inside_frame, rf_params->counter);
		drift(indices_inside_frame, rf_params->counter + 1);

		// find left outside particles and kick, drift them one more time
		int a = 0;
		for (int i = 0; i < beam->n_macroparticles; ++i) {
			if (beam->dt[i] < 0) {
				indices_left_outside[i] = beam->id[i] > 0;
				a++;
			} else {
				indices_left_outside[i] = false;
			}
		}
		if (a > 0) {
			// This will update only the indices_left_outside values
			//  need to test this
			for (int i = 0; i < beam->n_macroparticles; ++i) {
				beam->dt[i] += gp->t_rev[rf_params->counter + 1]
						* indices_left_outside[i];
			}
			kick(indices_left_outside, rf_params->counter);
			drift(indices_left_outside, rf_params->counter + 1);

		}
		// update inside, right outside particles

		set_periodicity();

	} else {
		kick(rf_params->counter);
		drift(rf_params->counter + 1);
	}
// cut particles by zeroing their id
// this way they will not be considered again in an update
// TODO maybe using lists and actually removing them can be
// faster, I need to test it (although I don't believe it)
	// TODO correct to horizontal
	if (dE_max > 0)
		horizontal_cut();
	rf_params->counter++;
// TODO I think there is no need to make any updates, counter is enough
}

inline void RingAndRfSection::horizontal_cut() {
// TODO resizing the array would be very expensive
// We should think of a new data structure or way of doing this
// Maybe the use of lists for beam-dE,dt would be a good solution
// In order to cut a particle we will 0 its id
	for (int i = 0; i < beam->n_macroparticles; ++i) {
		if (beam->dE[i] < -dE_max || beam->dE[i] > dE_max)
			beam->id[i] = 0;
	}
}

RingAndRfSection::RingAndRfSection(GeneralParameters *_gp,
		RfParameters *_rf_params, Beams *_beam, solver_type _solver,
		ftype *_PhaseLoop, ftype * _NoiseFB, bool _periodicity, ftype _dE_max,
		bool _rf_kick_interp, ftype* _Slices, ftype * _TotalInducedVoltage,
		int _n_threads) {
	this->gp = _gp;
	this->rf_params = _rf_params;
	this->beam = _beam;
	this->solver = _solver;
	this->PhaseLoop = _PhaseLoop;
	this->NoiseFB = _NoiseFB;
	this->periodicity = _periodicity;
	this->dE_max = _dE_max;
	this->rf_kick_interp = _rf_kick_interp;
	this->Slices = _Slices;
	this->TotalInducedVoltage = _TotalInducedVoltage;
	this->n_threads = _n_threads;
	this->indices_left_outside = new bool[beam->n_macroparticles];
	this->indices_right_outside = new bool[beam->n_macroparticles];
	this->indices_inside_frame = new bool[beam->n_macroparticles];
	for (int i = 0; i < beam->n_macroparticles; ++i) {
		indices_inside_frame[i] = true;
		indices_left_outside[i] = indices_right_outside[i] = false;
	}

//TODO fill unused eta arrays with zeros

	this->acceleration_kick = new ftype[rf_params->n_rf * (gp->n_turns)];
	for (int i = 0; i < rf_params->n_rf * gp->n_turns; ++i) {
		acceleration_kick[i] = -rf_params->E_increment[i];
	}

	if (solver != simple && solver != full) {
		dprintf(
				"ERROR: Choice of longitudinal solver not recognized! Aborting...");
		exit(-1);
	}

	if (gp->alpha_order > 1) {
		solver = full;
	}

	if (periodicity) {
		// TODO why not using this function anyway?
		int a = 0;
		for (int i = 0; i < beam->n_macroparticles; i++)
			a += beam->dt[i] < 0;
		if (a > 0) {
			dprintf("ERROR: condition beam.dt >= 0 not true!");
			exit(-1);
		}

		set_periodicity();
		// TODO Here we have a bad technique, we have to deallocate t_rev and allocate it again
		// TODO we will either do this by using vectors or not do it at all if there is no actual need
		gp->t_rev.push_back(gp->t_rev.back());
	}

}

inline void RingAndRfSection::set_periodicity() {

	for (int i = 0; i < beam->n_macroparticles; i++) {
		if (beam->dt[i] > gp->t_rev[rf_params->counter + 1]) {
			indices_right_outside[i] = beam->id[i] > 0;
			indices_inside_frame[i] = false;
		} else {
			indices_inside_frame[i] = beam->id[i] > 0;
			indices_right_outside[i] = false;
		}
	}
}

#endif /* TRACKERS_TRACKER_H_ */

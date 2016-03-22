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

// !!!!WARNING!!!!
// we only use beam->dt, dE
// bool arrays are used every time to update only the right values!!

class RingAndRfSection {

private:
public:
	ftype elapsed_time;

	bool *indices_right_outside;
	bool *indices_inside_frame;
	bool *indices_left_outside;
	inline void set_periodicity(const int start, const int end, const int turn);
	inline void kick(const bool *update, const int index, const int start,
			const int end);
	inline void kick(const ftype * __restrict__ beam_dt,
			ftype * __restrict__ beam_dE, const int n_rf,
			const ftype * __restrict__ voltage,
			const ftype * __restrict__ omega_RF,
			const ftype * __restrict__ phi_RF, const int n_macroparticles,
			const ftype acc_kick, const bool * __restrict__ update,
			const int start, const int end);
	inline void kick(const int index, const int start, const int end);
	inline void kick(const ftype * __restrict__ beam_dt,
			ftype * __restrict__ beam_dE, const int n_rf,
			const ftype * __restrict__ voltage,
			const ftype * __restrict__ omega_RF,
			const ftype * __restrict__ phi_RF, const int n_macroparticles,
			const ftype acc_kick, const int start, const int end);
	inline void drift(const bool *update, const int index, const int start,
			const int end);
	inline void drift(ftype * __restrict__ beam_dt,
			const ftype * __restrict__ beam_dE, const solver_type solver,
			const ftype T0, const ftype length_ratio, const int alpha_order,
			const ftype eta_zero, const ftype eta_one, const ftype eta_two,
			const ftype beta, const ftype energy, const int n_macroparticles,
			const bool * __restrict__ update, const int start, const int end);
	inline void drift(const int index, const int start, const int end);
	inline void drift(ftype * __restrict__ beam_dt,
			const ftype * __restrict__ beam_dE, const solver_type solver,
			const ftype T0, const ftype length_ratio, const int alpha_order,
			const ftype eta_zero, const ftype eta_one, const ftype eta_two,
			const ftype beta, const ftype energy, const int n_macroparticles,
			const int start, const int end);
	inline void track(const int start, const int end, const int turn);
	inline void horizontal_cut(const int start, const int end);
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
// First go the versions without periodicity

// Kick without periodicity
inline void RingAndRfSection::kick(const ftype * __restrict__ beam_dt,
		ftype * __restrict__ beam_dE, const int n_rf,
		const ftype * __restrict__ voltage, const ftype * __restrict__ omega_RF,
		const ftype * __restrict__ phi_RF, const int n_macroparticles,
		const ftype acc_kick, const int start, const int end) {

	//beam_dE[0] += 1;
// KICK

	int k = 0;
	for (int j = 0; j < n_rf; j++) {
//#pragma omp parallel for
		for (int i = start; i < end; i++) {
			beam_dE[i] += voltage[k]
					* vdt::fast_sin(omega_RF[k] * beam_dt[i] + phi_RF[k]);
		}
		k += gp->n_turns;
	}

// SYNCHRONOUS ENERGY CHANGE
//#pragma omp parallel for
	for (int i = start; i < end; i++)
		beam_dE[i] += acc_kick;

}

// kick with periodicity
inline void RingAndRfSection::kick(const ftype * __restrict__ beam_dt,
		ftype * __restrict__ beam_dE, const int n_rf,
		const ftype * __restrict__ voltage, const ftype * __restrict__ omega_RF,
		const ftype * __restrict__ phi_RF, const int n_macroparticles,
		const ftype acc_kick, const bool * __restrict__ update, const int start,
		const int end) {

// KICK
	int k = 0;
	for (int j = 0; j < n_rf; j++) {
		for (int i = start; i < end; i++)
			beam_dE[i] +=
					update[i] ?
							voltage[k]
									* vdt::fast_sin(
											omega_RF[k] * beam_dt[i]
													+ phi_RF[k]) :
							0;
		// what will I do with this k??
		k += gp->n_turns;
	}

// SYNCHRONOUS ENERGY CHANGE
	for (int i = start; i < end; i++)
		beam_dE[i] += update[i] ? acc_kick : 0;

}

//drift without periodicity
inline void RingAndRfSection::drift(ftype * __restrict__ beam_dt,
		const ftype * __restrict__ beam_dE, const solver_type solver,
		const ftype T0, const ftype length_ratio, const int alpha_order,
		const ftype eta_zero, const ftype eta_one, const ftype eta_two,
		const ftype beta, const ftype energy, const int n_macroparticles,
		const int start, const int end) {

//beam_dt[0] += 0.000001;

	int i;
	ftype T = T0 * length_ratio;

	if (solver == simple) {
		ftype coeff = eta_zero / (beta * beta * energy);
//#pragma omp parallel for
		for (i = start; i < end; i++)
			beam_dt[i] += T * coeff * beam_dE[i];
	}

	else {
		const ftype coeff = 1. / (beta * beta * energy);
		const ftype eta0 = eta_zero * coeff;
		const ftype eta1 = eta_one * coeff * coeff;
		const ftype eta2 = eta_two * coeff * coeff * coeff;

		if (alpha_order == 1)
			for (i = start; i < end; i++)
				beam_dt[i] += T * (1. / (1. - eta0 * beam_dE[i]) - 1.);
		else if (alpha_order == 2)
			for (i = start; i < end; i++)
				beam_dt[i] += T
						* (1.
								/ (1. - eta0 * beam_dE[i]
										- eta1 * beam_dE[i] * beam_dE[i]) - 1.);
		else
			for (i = start; i < end; i++)
				beam_dt[i] += T
						* (1.
								/ (1. - eta0 * beam_dE[i]
										- eta1 * beam_dE[i] * beam_dE[i]
										- eta2 * beam_dE[i] * beam_dE[i]
												* beam_dE[i]) - 1.);
	}

}

// drift with periodicity
inline void RingAndRfSection::drift(ftype * __restrict__ beam_dt,
		const ftype * __restrict__ beam_dE, const solver_type solver,
		const ftype T0, const ftype length_ratio, const int alpha_order,
		const ftype eta_zero, const ftype eta_one, const ftype eta_two,
		const ftype beta, const ftype energy, const int n_macroparticles,
		const bool * __restrict__ update, const int start, const int end) {

	int i;
	ftype T = T0 * length_ratio;

	if (solver == simple) {
		ftype coeff = eta_zero / (beta * beta * energy);

		for (i = start; i < end; i++)
			beam_dt[i] += update[i] ? T * coeff * beam_dE[i] : 0;
	}

	else {
		const ftype coeff = 1. / (beta * beta * energy);
		const ftype eta0 = eta_zero * coeff;
		const ftype eta1 = eta_one * coeff * coeff;
		const ftype eta2 = eta_two * coeff * coeff * coeff;

		if (alpha_order == 1)
			for (i = start; i < end; i++)
				beam_dt[i] +=
						update[i] ?
								T * (1. / (1. - eta0 * beam_dE[i]) - 1.) : 0;
		else if (alpha_order == 2)
			for (i = start; i < end; i++)
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
			for (i = start; i < end; i++)
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

inline void RingAndRfSection::track(const int start, const int end,
		const int turn) {
//omp_set_num_threads(n_threads);
//timespec begin, fin;

	if (periodicity) {
		// Change reference of all the particles on the right of the current
		// frame; these particles skip one kick and drift
		//for (int i = 0; i < beam->n_macroparticles; ++i) {
		for (int i = start; i < end; ++i) {
			beam->dt[i] -= indices_right_outside[i] * gp->t_rev[turn + 1];
		}
		// Synchronize the bunch with the particles that are on the right of
		// the current frame applying kick and drift to the bunch; after that
		// all the particle are in the new updated frame

		kick(indices_inside_frame, turn, start, end);
		drift(indices_inside_frame, turn + 1, start, end);

		// find left outside particles and kick, drift them one more time
		int a = 0;
		for (int i = start; i < end; ++i) {
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
			for (int i = start; i < end; ++i) {
				beam->dt[i] += gp->t_rev[turn + 1] * indices_left_outside[i];
			}
			kick(indices_left_outside, turn, start, end);
			drift(indices_left_outside, turn + 1, start, end);

		}
		// update inside, right outside particles

		set_periodicity(start, end, turn);

	} else {
		//get_time(begin);
		//dprintf("before kick\n");
		kick(turn, start, end);
		//dprintf("before drift\n");
		drift(turn + 1, start, end);
		//dprintf("after drift\n");

		//get_time(end);
		//elapsed_time += time_diff(end, begin);
	}
// cut particles by zeroing their id
// this way they will not be considered again in an update

	if (dE_max > 0)
		horizontal_cut(start, end);
	//rf_params->counter++;
	//printf("rf_params->counter = %d\n", rf_params->counter);

	//printf("rf_params->counter = %d\n", rf_params->counter);
}

inline void RingAndRfSection::horizontal_cut(const int start, const int end) {
// In order to cut a particle we will 0 its id
//#pragma omp parallel for
	for (int i = start; i < end; ++i) {
		if (beam->dE[i] < -dE_max || beam->dE[i] > dE_max)
			beam->id[i] = 0;
	}
}

RingAndRfSection::RingAndRfSection(GeneralParameters *_gp,
		RfParameters *_rf_params, Beams *_beam, solver_type _solver,
		ftype *_PhaseLoop, ftype * _NoiseFB, bool _periodicity, ftype _dE_max,
		bool _rf_kick_interp, ftype* _Slices, ftype * _TotalInducedVoltage,
		int _n_threads) {
	this->elapsed_time = 0;
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
		int a = 0;
		for (int i = 0; i < beam->n_macroparticles; i++)
			a += beam->dt[i] < 0;
		if (a > 0) {
			dprintf("ERROR: condition beam.dt >= 0 not true!");
			exit(-1);
		}

		set_periodicity(0, beam->n_macroparticles, 0);
		// we will either do this by using vectors or not do it at all if there is no actual need
		gp->t_rev.push_back(gp->t_rev.back());
	}

}

inline void RingAndRfSection::set_periodicity(const int start, const int end,
		const int turn) {

	for (int i = start; i < end; i++) {
		if (beam->dt[i] > gp->t_rev[turn + 1]) {
			indices_right_outside[i] = beam->id[i] > 0;
			indices_inside_frame[i] = false;
		} else {
			indices_inside_frame[i] = beam->id[i] > 0;
			indices_right_outside[i] = false;
		}
	}
}

inline void RingAndRfSection::kick(const int index, const int start,
		const int end) {
	kick(beam->dt, beam->dE, rf_params->n_rf, &rf_params->voltage[index],
			&rf_params->omega_RF[index], &rf_params->phi_RF[index],
			beam->n_macroparticles, acceleration_kick[index], start, end);
}

inline void RingAndRfSection::kick(const bool* update, const int index,
		const int start, const int end) {
	kick(beam->dt, beam->dE, rf_params->n_rf, &rf_params->voltage[index],
			&rf_params->omega_RF[index], &rf_params->phi_RF[index],
			beam->n_macroparticles, acceleration_kick[index], update, start,
			end);
}

inline void RingAndRfSection::drift(const bool *update, const int index,
		const int start, const int end) {
	drift(beam->dt, beam->dE, solver, gp->t_rev[index], rf_params->length_ratio,
			gp->alpha_order, rf_params->eta_0(index), rf_params->eta_1(index),
			rf_params->eta_2(index), rf_params->beta(index),
			rf_params->energy(index), beam->n_macroparticles, update, start,
			end);
}

inline void RingAndRfSection::drift(const int index, const int start,
		const int end) {
	drift(beam->dt, beam->dE, solver, gp->t_rev[index], rf_params->length_ratio,
			gp->alpha_order, rf_params->eta_0(index), rf_params->eta_1(index),
			rf_params->eta_2(index), rf_params->beta(index),
			rf_params->energy(index), beam->n_macroparticles, start, end);
}

#endif /* TRACKERS_TRACKER_H_ */

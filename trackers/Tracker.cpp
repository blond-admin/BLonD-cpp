/*
 * Tracker.cpp
 *
 *  Created on: Mar 22, 2016
 *      Author: kiliakis
 */

#include "Tracker.h"
//#include "../includes/utilities.h"

RingAndRfSection::~RingAndRfSection() {
	//if (this->PhaseLoop != NULL) delete this->PhaseLoop;
	//if (this->NoiseFB != NULL) delete this->NoiseFB;
	//if (this->TotalInducedVoltage != NULL) delete this->TotalInducedVoltage;
	//if (this->acceleration_kick != NULL) delete this->acceleration_kick;
	delete_array(this->PL);
	delete_array(this->NoiseFB);
	delete_array(this->TotalInducedVoltage);
	delete_array(this->acceleration_kick);
	delete_array(this->Slices);
	delete_array(this->indices_left_outside);
	delete_array(this->indices_right_outside);
	delete_array(this->indices_inside_frame);

}

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
		k += GP->n_turns;
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
		k += GP->n_turns;
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

void RingAndRfSection::track(const int start, const int end) {
//omp_set_num_threads(n_threads);
//timespec begin, fin;

	if (periodicity) {
		// Change reference of all the particles on the right of the current
		// frame; these particles skip one kick and drift
		//for (int i = 0; i < Beam->n_macroparticles; ++i) {
		for (int i = start; i < end; ++i) {
			Beam->dt[i] -= indices_right_outside[i]
					* GP->t_rev[RfP->counter + 1];
		}
		// Synchronize the bunch with the particles that are on the right of
		// the current frame applying kick and drift to the bunch; after that
		// all the particle are in the new updated frame

		kick(indices_inside_frame, RfP->counter, start, end);
		drift(indices_inside_frame, RfP->counter + 1, start, end);

		// find left outside particles and kick, drift them one more time
		int a = 0;
		for (int i = start; i < end; ++i) {
			if (Beam->dt[i] < 0) {
				indices_left_outside[i] = Beam->id[i] > 0;
				a++;
			} else {
				indices_left_outside[i] = false;
			}
		}
		if (a > 0) {
			// This will update only the indices_left_outside values
			//  need to test this
			for (int i = start; i < end; ++i) {
				Beam->dt[i] += GP->t_rev[RfP->counter + 1]
						* indices_left_outside[i];
			}
			kick(indices_left_outside, RfP->counter, start, end);
			drift(indices_left_outside, RfP->counter + 1, start, end);

		}
		// update inside, right outside particles

		set_periodicity(start, end);

	} else {
		//get_time(begin);
		//dprintf("before kick\n");
		kick(RfP->counter, start, end);
		//dprintf("before drift\n");
		drift(RfP->counter + 1, start, end);
		//dprintf("after drift\n");

		//get_time(end);
		//elapsed_time += time_diff(end, begin);
	}
// cut particles by zeroing their id
// this way they will not be considered again in an update

	if (dE_max > 0)
		horizontal_cut(start, end);
	//RfP->counter++;
}

inline void RingAndRfSection::horizontal_cut(const int start, const int end) {
// In order to cut a particle we will 0 its id
//#pragma omp parallel for
	for (int i = start; i < end; ++i) {
		if (Beam->dE[i] < -dE_max || Beam->dE[i] > dE_max)
			Beam->id[i] = 0;
	}
}

RingAndRfSection::RingAndRfSection(solver_type _solver, PhaseLoop *_PhaseLoop,
		ftype * _NoiseFB, bool _periodicity, ftype _dE_max,
		bool _rf_kick_interp, ftype* _Slices, ftype * _TotalInducedVoltage) {
	this->elapsed_time = 0;
	this->solver = _solver;
	this->PL = _PhaseLoop;
	this->NoiseFB = _NoiseFB;
	this->periodicity = _periodicity;
	this->dE_max = _dE_max;
	this->rf_kick_interp = _rf_kick_interp;
	this->Slices = _Slices;
	this->TotalInducedVoltage = _TotalInducedVoltage;
	//this->n_threads = _n_threads;
	this->indices_left_outside = new bool[Beam->n_macroparticles];
	this->indices_right_outside = new bool[Beam->n_macroparticles];
	this->indices_inside_frame = new bool[Beam->n_macroparticles];
	for (int i = 0; i < Beam->n_macroparticles; ++i) {
		indices_inside_frame[i] = true;
		indices_left_outside[i] = indices_right_outside[i] = false;
	}

	this->acceleration_kick = new ftype[RfP->n_rf * (GP->n_turns)];
	for (int i = 0; i < RfP->n_rf * GP->n_turns; ++i) {
		acceleration_kick[i] = -RfP->E_increment[i];
	}

	if (solver != simple && solver != full) {
		dprintf(
				"ERROR: Choice of longitudinal solver not recognized! Aborting...");
		exit(-1);
	}

	if (GP->alpha_order > 1) {
		solver = full;
	}

	if (periodicity) {
		int a = 0;
		for (int i = 0; i < Beam->n_macroparticles; i++)
			a += Beam->dt[i] < 0;
		if (a > 0) {
			dprintf("ERROR: condition Beam.dt >= 0 not true!");
			exit(-1);
		}

		set_periodicity(0, Beam->n_macroparticles);
		// we will either do this by using vectors or not do it at all if there is no actual need
		GP->t_rev.push_back(GP->t_rev.back());
	}

}

inline void RingAndRfSection::set_periodicity(const int start, const int end) {

	for (int i = start; i < end; i++) {
		if (Beam->dt[i] > GP->t_rev[RfP->counter + 1]) {
			indices_right_outside[i] = Beam->id[i] > 0;
			indices_inside_frame[i] = false;
		} else {
			indices_inside_frame[i] = Beam->id[i] > 0;
			indices_right_outside[i] = false;
		}
	}
}

inline void RingAndRfSection::kick(const int index, const int start,
		const int end) {
	kick(Beam->dt, Beam->dE, RfP->n_rf, &RfP->voltage[index],
			&RfP->omega_RF[index], &RfP->phi_RF[index], Beam->n_macroparticles,
			acceleration_kick[index], start, end);
}

inline void RingAndRfSection::kick(const bool* update, const int index,
		const int start, const int end) {
	kick(Beam->dt, Beam->dE, RfP->n_rf, &RfP->voltage[index],
			&RfP->omega_RF[index], &RfP->phi_RF[index], Beam->n_macroparticles,
			acceleration_kick[index], update, start, end);
}

inline void RingAndRfSection::drift(const bool *update, const int index,
		const int start, const int end) {
	drift(Beam->dt, Beam->dE, solver, GP->t_rev[index], RfP->length_ratio,
			GP->alpha_order, RfP->eta_0(index), RfP->eta_1(index),
			RfP->eta_2(index), RfP->beta(index), RfP->energy(index),
			Beam->n_macroparticles, update, start, end);
}

inline void RingAndRfSection::drift(const int index, const int start,
		const int end) {
	drift(Beam->dt, Beam->dE, solver, GP->t_rev[index], RfP->length_ratio,
			GP->alpha_order, RfP->eta_0(index), RfP->eta_1(index),
			RfP->eta_2(index), RfP->beta(index), RfP->energy(index),
			Beam->n_macroparticles, start, end);
}

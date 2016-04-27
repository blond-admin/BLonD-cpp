/*
 * TC8_Phase_loop.cpp
 *
 *  Created on: Apr 11, 2016
 *      Author: kiliakis
 */

#include "globals.h"
#include "utilities.h"
#include "math_functions.h"
#include <omp.h>
#include <stdio.h>
#include "../input_parameters/GeneralParameters.h"
#include "../input_parameters/RfParameters.h"
#include "../beams/Beams.h"
#include "../beams/Slices.h"
#include "../beams/Distributions.h"
#include "../trackers/Tracker.h"
#include "../llrf/PhaseLoop.h"
// Simulation parameters --------------------------------------------------------

// Bunch parameters
const ftype tau_0 = 0.4e-9;          // Initial bunch length, 4 sigma [s]

// Machine and RF parameters
//const long p_i = 450e9;          // Synchronous momentum [eV/c]
//const ftype p_f = 460.005e9;          // Synchronous momentum, final
//const long h = 35640;          // Harmonic number
//const ftype dphi = 0;          // Phase modulation/offset

// Global variables
GeneralParameters *GP;
Beams *Beam;
Slices *Slice;
RfParameters *RfP;
int n_threads = 1;

// Simulation setup -------------------------------------------------------------
int main(int argc, char **argv) {

	// Beam parameters
	//particle_type particle_type = proton;
	int n_macroparticles = 100000;
	int n_particles = 0;
	int N_slices = 500;

	// Machine and RF params
	ftype radious = 25; // [m]
	ftype gamma_transition = 4.076750841; // [1]
	ftype alpha = 1 / (gamma_transition * gamma_transition); // [1]
	ftype C = 2 * pi * radious;
	int n_turns = 500;

	// Cavities parameters
	int n_rf_systems = 1;
	int harmonic_numbers_1 = 1;
	ftype voltage_1 = 8000; // [V]
	ftype phi_offset_1 = 0; //[rad]

	int alpha_order = 1;

	ftype *momentum = new ftype[n_turns + 1];
	std::fill_n(momentum, n_turns + 1, 310891054.809);

	ftype *alpha_array = new ftype[(alpha_order + 1) * n_rf_systems];
	std::fill_n(alpha_array, (alpha_order + 1) * n_rf_systems, alpha);

	ftype *C_array = new ftype[n_rf_systems];
	C_array[0] = C;

	// Environmental variables
	n_turns = atoi(GETENV("N_TURNS")) ? atoi(GETENV("N_TURNS")) : n_turns;
	n_macroparticles =
			atoi(GETENV("N_PARTICLES")) ?
					atoi(GETENV("N_PARTICLES")) : n_macroparticles;
	N_slices = atoi(GETENV("N_SLICES")) ? atoi(GETENV("N_SLICES")) : N_slices;
	n_threads =
			atoi(GETENV("N_THREADS")) ? atoi(GETENV("N_THREADS")) : n_threads;
	omp_set_num_threads(n_threads);

	printf("Setting up the simulation...\n\n");
	printf("Number of turns: %d\n", n_turns);
	printf("Number of macro-particles: %d\n", n_macroparticles);
	printf("Number of Slices: %d\n", N_slices);

#pragma omp parallel
	{
		if (omp_get_thread_num() == 0)
			printf("Number of openmp threads: %d\n", omp_get_num_threads());
	}

	timespec begin, end;
	get_time(begin);

	ftype *h_array = new ftype[n_rf_systems * (n_turns + 1)];
	std::fill_n(h_array, (n_turns + 1) * n_rf_systems, harmonic_numbers_1);

	ftype *V_array = new ftype[n_rf_systems * (n_turns + 1)];
	std::fill_n(V_array, (n_turns + 1) * n_rf_systems, voltage_1);

	ftype *dphi_array = new ftype[n_rf_systems * (n_turns + 1)];
	std::fill_n(dphi_array, (n_turns + 1) * n_rf_systems, phi_offset_1);

	// TODO variables must be in the correct format (arrays for all)
	// fix this with builder design pattern + multiple constructors

	GP = new GeneralParameters(n_turns, C_array, alpha_array, alpha_order,
			momentum, proton);
	//printf("ok\n");
	Beam = new Beams(n_macroparticles, n_particles);
	//printf("ok\n");

	Slice = new Slices(N_slices, -pi, pi, rad);

	RfP = new RfParameters(n_rf_systems, h_array, V_array, dphi_array);

	ftype RL_gain[2] = { 0, 0 };
	ftype *PL_gain = new ftype[n_turns];
	std::fill_n(PL_gain, n_turns, 1.0 / (25e-6));
	//dump(PL_gain, 100, "PL_gain\n");
	PhaseLoop *psb = new PSB(PL_gain, RL_gain, 10e-6, 7);
	//dump(psb->gain, 100, "psb->gain\n");

	RingAndRfSection *long_tracker = new RingAndRfSection(simple, psb, NULL,
			false);

	//longitudinal_bigaussian(tau_0 / 4, 0, 1, false);

	//dump(Slice->bin_centers, N_slices, "bin_centers\n");
	//dump(Beam->dE, 10, "dE\n");

	double slice_time = 0, track_time = 0;

	timespec begin_t, end_t;

#pragma omp parallel
	{
		int id = omp_get_thread_num();
		int threads = omp_get_num_threads();
		int tile = std::ceil(1.0 * n_macroparticles / threads);
		int start = id * tile;
		int end = std::min(start + tile, n_macroparticles);
		//printf("id, threads, tile, start, end = %d, %d, %d, %d, %d\n", id,
		//		threads, tile, start, end);
		for (int i = 0; i < n_turns; ++i) {
			printf("Turn %d\n", i);

#ifdef TIMING
			if (id == 0) get_time(begin_t);
#endif
			//dump(Beam->dE, 1, "dE\n");

			long_tracker->track(start, end);

#pragma omp barrier

#ifdef TIMING
			if (id == 0) track_time += time_elapsed(begin_t);
			if (id == 0) get_time(begin_t);
#endif
			Slice->track(start, end);

#pragma omp barrier

#ifdef TIMING
			if (id == 0) slice_time += time_elapsed(begin_t);
#endif

#pragma omp single
			{
				Slice->fwhm();
				psb->track();
#ifdef PRINT_RESULTS
				if (i % 1000 == 0) {
					printf("bl_fwhm\n%.15lf\n", Slice->bl_fwhm);
					printf("bp_fwhm\n%.15lf\n", Slice->bp_fwhm);
					//dump(Slice->n_macroparticles, N_slices,
					//		"n_macroparticles\n");
				}
#endif
				RfP->counter++;
			}
			//beam->losses_longitudinal_cut(beam->dt, 0, 2.5e-9, beam->id);
		}
	}

//printf("Total simulation time: %.10lf\n", long_tracker->elapsed_time);
//printf("Time/turn : %.10lf\n", long_tracker->elapsed_time / n_turns);
	ftype result = mymath::trapezoid(Slice->n_macroparticles,
			Slice->bin_centers, Slice->n_slices);
	printf("result = %e\n", result);

	get_time(end);
#ifdef TIMING
	print_time("Simulation Time", begin, end);
	double total_time = track_time + slice_time;
	printf("Track time : %.4lf ( %.2lf %% )\n", track_time,
			100 * track_time / total_time);
	printf("Slice time : %.4lf ( %.2lf %% )\n", slice_time,
			100 * slice_time / total_time);
#endif

	dump(Beam->dE, 10, "dE\n");
	dump(Beam->dt, 10, "dt\n");
	dump(Slice->n_macroparticles, 10, "n_macroparticles\n");
	delete Slice;
	delete long_tracker;
	delete RfP;
	delete GP;
	delete Beam;
	printf("Done!\n");

}


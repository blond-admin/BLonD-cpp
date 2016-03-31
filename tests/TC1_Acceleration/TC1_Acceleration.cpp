/*
 * TC1_Acceleration.cpp
 *
 *  Created on: Mar 9, 2016
 *      Author: kiliakis
 */
#include "globals.h"
#include "utilities.h"
#include "math_functions.h"
#include "../../input_parameters/GeneralParameters.h"
#include "../../input_parameters/RfParameters.h"
#include "../../beams/Beams.h"
#include "../../beams/Slices.h"
#include "../../beams/Distributions.h"
#include "../../trackers/Tracker.h"
#include <omp.h>
#include <stdio.h>
//#include <gsl/gsl_const_mksa.h>
//#include <gsl/gsl_math.h>

// Simulation parameters --------------------------------------------------------
// Bunch parameters
const long N_b = 1e9;           // Intensity
const ftype tau_0 = 0.4e-9;          // Initial bunch length, 4 sigma [s]

// Machine and RF parameters
const ftype C = 26658.883;          // Machine circumference [m]
const long p_i = 450e9;          // Synchronous momentum [eV/c]
const ftype p_f = 460.005e9;          // Synchronous momentum, final
const long h = 35640;          // Harmonic number
const ftype V = 6e6;          // RF voltage [V]
const ftype dphi = 0;          // Phase modulation/offset
const ftype gamma_t = 55.759505;          // Transition gamma
const ftype alpha = 1.0 / gamma_t / gamma_t;    // First order mom. comp. factor
const int alpha_order = 1;
const int n_sections = 1;
// Tracking details
int N_t = 50000;    // Number of turns to track
int N_p = 10000;         // Macro-particles

int n_threads = 1;

GeneralParameters *GP;
Beams *Beam;
Slices *Slice;
RfParameters *RfP;

// Simulation setup -------------------------------------------------------------
int main(int argc, char **argv) {

	N_t = atoi(GETENV("N_TURNS")) ? atoi(GETENV("N_TURNS")) : N_t;
	N_p = atoi(GETENV("N_PARTICLES")) ? atoi(GETENV("N_PARTICLES")) : N_p;
	n_threads =
			atoi(GETENV("N_THREADS")) ? atoi(GETENV("N_THREADS")) : n_threads;
	omp_set_num_threads(n_threads);

	// Number of tasks is either N_TASKS if specified or n_threads (1 task / thread) if not
	//int N_tasks = atoi(GETENV("N_TASKS")) ? atoi(GETENV("N_TASKS")) : n_threads;

	printf("Setting up the simulation...\n\n");
	printf("Number of turns: %d\n", N_t);
	printf("Number of macro-particles: %d\n", N_p);
	//printf("Number of Tasks: %d\n", N_tasks);

	/// initializations
#pragma omp parallel
	{
		if (omp_get_thread_num() == 0)
			printf("Number of openmp threads = %d\n", omp_get_num_threads());
	}

	timespec begin, end;
	double slice_time = 0, track_time = 0;
	get_time(begin);

	ftype *momentum = new ftype[N_t + 1];
	mymath::linspace(momentum, p_i, p_f, N_t + 1);

	ftype *alpha_array = new ftype[(alpha_order + 1) * n_sections];
	std::fill_n(alpha_array, (alpha_order + 1) * n_sections, alpha);

	ftype *C_array = new ftype[n_sections];
	C_array[0] = C;

	ftype *h_array = new ftype[n_sections * (N_t + 1)];
	std::fill_n(h_array, (N_t + 1) * n_sections, h);

	ftype *V_array = new ftype[n_sections * (N_t + 1)];
	//ftype * V_array = (ftype *) aligned_malloc(
	//		sizeof(ftype) * n_sections * (N_t + 1));
	std::fill_n(V_array, (N_t + 1) * n_sections, V);

	ftype *dphi_array = new ftype[n_sections * (N_t + 1)];
	std::fill_n(dphi_array, (N_t + 1) * n_sections, dphi);

// TODO variables must be in the correct format (arrays for all)

	GP = new GeneralParameters(N_t, C_array, alpha_array, alpha_order, momentum,
			proton);

	Beam = new Beams(N_p, N_b);

	RfP = new RfParameters(n_sections, h_array, V_array, dphi_array);

	RingAndRfSection *long_tracker = new RingAndRfSection();

	longitudinal_bigaussian(tau_0 / 4, 0, 1, false);

	Slice = new Slices(100);

#pragma omp parallel
	{
		int id = omp_get_thread_num();
		int threads = omp_get_num_threads();
		int tile = std::ceil(1.0 * N_p / threads);
		int start = id * tile;
		int end = std::min(start + tile, N_p);
		//printf("id, threads, tile, start, end = %d, %d, %d, %d, %d\n", id,
		//		threads, tile, start, end);
		for (int i = 0; i < N_t; ++i) {
			//get_time(begin);
			long_tracker->track(start, end);

#pragma omp barrier

			Slice->track(start, end);

#pragma omp barrier
#pragma omp single
			{
				RfP->counter++;
			}
			//slice_time += time_elapsed(begin);

			//beam->losses_longitudinal_cut(beam->dt, 0, 2.5e-9, beam->id);
		}
	}

//printf("Total simulation time: %.10lf\n", long_tracker->elapsed_time);
//printf("Time/turn : %.10lf\n", long_tracker->elapsed_time / N_t);
	get_time(end);
	print_time("Total Simulation Time", begin, end);
//print_time("Total Track Time", track_time);
//print_time("Total Slice Time", slice_time);

	dump(Beam->dE, 10, "beam->dE\n");
	dump(Beam->dt, 10, "beam->dt\n");
	dump(Slice->n_macroparticles, 10, "Slice->n_macroparticles\n");
	delete Slice;
	delete long_tracker;
	delete RfP;
	delete GP;
	delete Beam;
	printf("Done!\n");

}

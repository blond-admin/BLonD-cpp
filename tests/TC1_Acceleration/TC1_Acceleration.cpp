/*
 * TC1_Acceleration.cpp
 *
 *  Created on: Mar 9, 2016
 *      Author: kiliakis
 */

#include "../../includes/utilities.h"
#include "../../includes/math_functions.h"
#include "../../input_parameters/GeneralParameters.h"
#include "../../input_parameters/RfParameters.h"
#include "../../beams/Beams.h"
#include "../../beams/Distributions.h"
#include "../../trackers/Tracker.h"
#include <omp.h>
#include <stdio.h>
//#include "../../includes/globals.h"

//GeneralParameters *GP;
//Beams *Beam;
//RfParameters *RfP;

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
int N_t = 5000;    // Number of turns to track
int N_p = 10;         // Macro-particles
int n_threads = 1;

// Simulation setup -------------------------------------------------------------
int main(int argc, char **argv) {

	N_t = atoi(GETENV("N_TURNS")) ? atoi(GETENV("N_TURNS")) : N_t;
	N_p = atoi(GETENV("N_PARTICLES")) ? atoi(GETENV("N_PARTICLES")) : N_p;
	n_threads =
			atoi(GETENV("N_THREADS")) ? atoi(GETENV("N_THREADS")) : n_threads;
	omp_set_num_threads(n_threads);

	printf("Setting up the simulation...\n\n");
	printf("Number of turns: %d\n", N_t);
	printf("Number of macro-particles: %d\n", N_p);
	
	/// initializations
#pragma omp parallel
	{
		if (omp_get_thread_num() == 0)
			printf("Number of openmp threads = %d\n", omp_get_num_threads());
	}

	timespec begin, end;
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
	GeneralParameters *general_params = new GeneralParameters(N_t, C_array,
			alpha_array, alpha_order, momentum, proton);

// TODO maybe general_params, beam, and RfParameters could be global?

	Beams *beam = new Beams(general_params, N_p, N_b);

	RfParameters *rf_params = new RfParameters(general_params, beam, n_sections,
			h_array, V_array, dphi_array);

	RingAndRfSection *long_tracker = new RingAndRfSection(general_params,
			rf_params, beam);

	longitudinal_bigaussian(general_params, rf_params, beam, tau_0 / 4, 0, 1,
			false);
	//dump(general_params->beta, N_t+1, "beta\n");
	//dump(rf_params->harmonic, N_t+1, "harmonic\n");
	//dump(beam->dE, N_p, "beam->dE\n");
	//dump(beam->dt, N_p, "beam->dt\n");
	//dump(long_tracker->acceleration_kick, N_p, "acc_kick");
	//dump(rf_params->E_increment, n_sections * (N_t), "E_increment");
	//dump(rf_params->Qs, n_sections * (N_t + 1), "Qs");
	//dump(general_params.eta_0, N_t + 1, "eta_0");

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
			//dump(beam->dE, beam->n_macroparticles, "beam->dE");
			long_tracker->track(start, end, i);
			//beam->losses_longitudinal_cut(beam->dt, 0, 2.5e-9, beam->id);
		}
		//printf("rf_params->counter = %d\n", rf_params->counter);
	}
	//printf("Total simulation time: %.10lf\n", long_tracker->elapsed_time);
	printf("Time/turn: %.10lf\n", long_tracker->elapsed_time / N_t);
	get_time(end);
	print_time("Total Simulation Time", begin, end);
	
	dump(beam->dE, 10, "beam->dE\n");
	dump(beam->dt, 10, "beam->dt\n");
	
	delete long_tracker;
	delete rf_params;
	delete general_params;
	delete beam;
	printf("Done!\n");

}

/*
 * LHC_restart.cpp
 *
 *  Created on: Apr 12, 2016
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

const int N_b = 1.2e9;          // Intensity
int N_p = 100000;         // Macro-particles

// Machine and RF parameters
const float C = 26658.883;        // Machine circumference [m]
const int h = 35640;            // Harmonic number
const float dphi = 0.;            // Phase modulation/offset
const float gamma_t = 55.759505;  // Transition gamma
const float alpha = 1. / gamma_t / gamma_t;     // First order mom. comp. factor

// Tracking details
int N_t = 1000000;        // Number of turns to track; full ramp: 8700001
int dt_plt = 100000;      // Time steps between plots
int dt_mon = 1;           // Time steps between monitoring
int dt_save = 1000000;    // Time steps between saving coordinates
int bl_target = 1.25e-9;  // 4 sigma r.m.s. target bunch length in [s]

int N_slices = 151;
const std::string datafiles =
		"/afs/cern.ch/work/k/kiliakis/workspace/BLonD-minimal-cpp/datafiles/";

// Global variables
GeneralParameters *GP;
Beams *Beam;
Slices *Slice;
RfParameters *RfP;
int n_threads = 1;
//const int size = 14e6;
// Simulation setup -------------------------------------------------------------
int main(int argc, char **argv) {

	// Environmental variables
	N_t = atoi(GETENV("N_TURNS")) ? atoi(GETENV("N_TURNS")) : N_t;
	N_p = atoi(GETENV("N_PARTICLES")) ? atoi(GETENV("N_PARTICLES")) : N_p;
	N_slices = atoi(GETENV("N_SLICES")) ? atoi(GETENV("N_SLICES")) : N_slices;
	n_threads =
			atoi(GETENV("N_THREADS")) ? atoi(GETENV("N_THREADS")) : n_threads;
	omp_set_num_threads(n_threads);

	printf("Setting up the simulation...\n\n");
	printf("Number of turns: %d\n", N_t);
	printf("Number of macro-particles: %d\n", N_p);
	printf("Number of Slices: %d\n", N_slices);

#pragma omp parallel
	{
		if (omp_get_thread_num() == 0)
			printf("Number of openmp threads: %d\n", omp_get_num_threads());
	}

	printf("Setting up the simulation..\n");

	//ftype *ps = new ftype[N_t + 1];
	std::vector < ftype > v;
	read_vector_from_file(v, datafiles + "LHC_momentum_programme_6.5TeV.dat");
	//for (int i = 0; i < 100; ++i) {
	//	std::cout << "v[" << i << "] = " << v[i] << "\n";
	//}
	v.erase(v.begin(), v.begin() + 13e6);
	//std::cout << "vector size is " << v.size() << "\n";
	int remaining = N_t + 1 - v.size();
	for (int i = 0; i < remaining; ++i) {
		v.push_back(6.5e12);
	}
	assert((int) v.size() == N_t + 1);
	ftype *ps = &v[0];	//new ftype[v.size()];
	printf("Length of ps is %lu\n", v.size());
	printf("Flat top momentum %.4e eV\n", ps[N_t]);

	ftype *V_array = new ftype[N_t + 1];
	mymath::linspace(V_array, 6e6, 10e6, 13563374, 13e6);
	std::fill_n(&V_array[563374], 436627, 10e6);
	printf("Length of V is %d\n", N_t + 1);
	printf("Flat top voltage %.4e eV\n", V_array[N_t]);
	printf("Momentum and voltage loaded\n");

	// Define general parameters
	int alpha_order = 1;
	int n_sections = 1;
	ftype *alpha_array = new ftype[(alpha_order + 1) * n_sections];
	std::fill_n(alpha_array, (alpha_order + 1) * n_sections, alpha);

	ftype *C_array = new ftype[n_sections];
	C_array[0] = C;

	GP = new GeneralParameters(N_t, C_array, alpha_array, alpha_order, ps,
			proton);
	printf("General parameters set...\n");

	// Define rf_params
	ftype *dphi_array = new ftype[n_sections * (N_t + 1)];
	std::fill_n(dphi_array, (N_t + 1) * n_sections, dphi);

	ftype *h_array = new ftype[n_sections * (N_t + 1)];
	std::fill_n(h_array, (N_t + 1) * n_sections, h);

	RfP = new RfParameters(n_sections, h_array, V_array, dphi_array);

	printf("RF parameters set...\n");

	// Define beam and distribution: Load matched, filamented distribution
	Beam = new Beams(N_p, N_b);
	std::vector < ftype > v2;
	read_vector_from_file(v2, datafiles + "coords_13000001.dat");
	printf("v2 size is %lu\n", v2.size());
	int k = 0;
	for (unsigned int i = 0; i < v2.size(); i += 3) {
		Beam->dt[k] = v2[i];
		Beam->dE[k] = v2[i + 1];
		Beam->id[k] = v2[i + 2];
		k++;
	}

	Slice = new Slices(N_slices, 0, -0.5e-9, 3e-9);
	printf("Beam generated, slices set...\n");
	//std::cout << "vector size is " << v.size() << "\n";
	//dump(ps, size, "ps\n");
	// Beam parameters
	//particle_type particle_type = proton;

	/*
	 // Machine and RF params
	 ftype radious = 25; // [m]
	 ftype gamma_transition = 4.076750841; // [1]
	 ftype alpha = 1 / (gamma_transition * gamma_transition); // [1]
	 ftype C = 2 * pi * radious;
	 int N_t = 500;

	 // Cavities parameters
	 int n_sections = 1;
	 int harmonic_numbers_1 = 1;
	 ftype voltage_1 = 8000; // [V]
	 ftype phi_offset_1 = 0; //[rad]

	 int alpha_order = 1;

	 ftype *momentum = new ftype[N_t + 1];
	 std::fill_n(momentum, N_t + 1, 310891054.809);

	 ftype *alpha_array = new ftype[(alpha_order + 1) * n_sections];
	 std::fill_n(alpha_array, (alpha_order + 1) * n_sections, alpha);

	 ftype *C_array = new ftype[n_sections];
	 C_array[0] = C;

	 timespec begin, end;
	 get_time(begin);

	 ftype *h_array = new ftype[n_sections * (N_t + 1)];
	 std::fill_n(h_array, (N_t + 1) * n_sections, harmonic_numbers_1);

	 ftype *V_array = new ftype[n_sections * (N_t + 1)];
	 std::fill_n(V_array, (N_t + 1) * n_sections, voltage_1);

	 ftype *dphi_array = new ftype[n_sections * (N_t + 1)];
	 std::fill_n(dphi_array, (N_t + 1) * n_sections, phi_offset_1);

	 // TODO variables must be in the correct format (arrays for all)
	 // fix this with builder design pattern + multiple constructors

	 GP = new GeneralParameters(N_t, C_array, alpha_array, alpha_order, momentum,
	 proton);
	 //printf("ok\n");
	 Beam = new Beams(N_p, n_particles);
	 //printf("ok\n");

	 Slice = new Slices(N_slices, -pi, pi, rad);

	 RfP = new RfParameters(n_sections, h_array, V_array, dphi_array);

	 ftype RL_gain[2] = { 0, 0 };
	 ftype *PL_gain = new ftype[N_t];
	 std::fill_n(PL_gain, N_t, 1.0 / (25e-6));
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
	 int tile = std::ceil(1.0 * N_p / threads);
	 int start = id * tile;
	 int end = std::min(start + tile, N_p);
	 //printf("id, threads, tile, start, end = %d, %d, %d, %d, %d\n", id,
	 //		threads, tile, start, end);
	 for (int i = 0; i < N_t; ++i) {
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
	 //dump(Slice->N_p, N_slices,
	 //		"N_p\n");
	 }
	 #endif
	 RfP->counter++;
	 }
	 //beam->losses_longitudinal_cut(beam->dt, 0, 2.5e-9, beam->id);
	 }
	 }

	 //printf("Total simulation time: %.10lf\n", long_tracker->elapsed_time);
	 //printf("Time/turn : %.10lf\n", long_tracker->elapsed_time / N_t);
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
	 dump(Slice->n_macroparticles, 10, "N_p\n");
	 
	 delete Slice;
	 delete long_tracker;
	 delete RfP;
	 delete GP;
	 delete Beam;
	 */
	printf("Done!\n");

}


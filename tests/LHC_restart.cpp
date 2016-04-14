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
		"/afs/cern.ch/work/k/kiliakis/testcases/synchroLoop/";

// Global variables
GeneralParameters *GP;
Beams *Beam;
Slices *Slice;
RfParameters *RfP;
int n_threads = 1;
//const int size = 14e6;
const int from_line = 436626;
// Simulation setup -------------------------------------------------------------
int main(int argc, char **argv) {
	timespec begin, end;
	get_time(begin);

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
	read_vector_from_file(v, datafiles + "LHC_momentum_programme_folder/xan");
	//for (int i = 0; i < 100; ++i) {
	//	std::cout << "v[" << i << "] = " << v[i] << "\n";
	//}

	// optional
	v.erase(v.begin(), v.begin() + from_line);

	std::cout << "vector size is " << v.size() << "\n";
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

	/*
	 dump(RfP->harmonic, 3, "harmonic\n");
	 dprintf("charge %lf\n", GP->charge);
	 dump(RfP->voltage, 3, "voltage\n");
	 dump(RfP->phi_s, 3, "phi_s\n");
	 dump(GP->eta_0, 3, "eta_0\n");
	 dump(RfP->Qs, 3, "Qs\n");
	 dump(GP->omega_rev, 3, "omega_rev\n");
	 dump(RfP->omega_s0, 3, "omega_s0\n");
	 */

	// Define beam and distribution: Load matched, filamented distribution
	Beam = new Beams(N_p, N_b);
	std::vector < ftype > v2;
	read_vector_from_file(v2, datafiles + "coords_13000001.dat");
	//printf("v2 size is %lu\n", v2.size());
	int k = 0;
	for (unsigned int i = 0; i < v2.size(); i += 3) {
		Beam->dt[k] = v2[i] * 1e-9; // [s]
		Beam->dE[k] = v2[i + 1] * 1e6; // [eV]
		Beam->id[k] = v2[i + 2];
		k++;
	}

	Slice = new Slices(N_slices, 0, -0.5e-9, 3e-9);
	printf("Beam generated, slices set...\n");
	//dump(Slice->bin_centers, 10, "bin_centers\n");
	// Define phase loop and frequency loop gain
	ftype PL_gain = 1 / (5 * GP->t_rev[0]);
	ftype SL_gain = PL_gain / 10;
	ftype *PL_gain_array = new ftype[N_t + 1];
	std::fill_n(PL_gain_array, N_t + 1, PL_gain);
	LHC *PL = new LHC(PL_gain_array, SL_gain);

	printf("\tPL gain is %.4e 1/s for initial turn T0 = %.4e s\n", PL->gain[0],
			GP->t_rev[0]);
	printf("\tSL gain is %.4e turns\n", PL->gain2);
	printf("\tOmega_s0 = %.4e s at flat bottom, %.4e s at flat top\n",
			RfP->omega_s0[0], RfP->omega_s0[N_t]);
	printf("\tSL a_i = %.4f a_f = %.4f\n", PL->lhc_a[0], PL->lhc_a[N_t]);
	printf("\tSL t_i = %.4f t_f = %.4f\n", PL->lhc_t[0], PL->lhc_t[N_t]);

	// Injecting noise in the cavity, PL on
	RingAndRfSection *long_tracker = new RingAndRfSection(simple, PL);
	printf("PL, SL, and tracker set...\n");
	//dump(Slice->bin_centers, 10, "bin_centers\n");
	double slice_time = 0, track_time = 0;

	timespec begin_t;
	printf("Map set\n");

	printf("Initial mean bunch position %.4e s\n", Beam->mean_dt);
	printf("Initial four-times r.m.s. bunch length %.4e s\n",
			4. * Beam->sigma_dt);
	//print("Initial Gaussian bunch length %.4e ns" %slices.bl_gauss

	printf("Ready for tracking!\n");

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
			printf("   Beam energy %.6e eV\n", GP->energy[0]);
			printf("   RF phase %.6e rad\n", RfP->dphi_RF[0]);
			printf("   PL phase correction %.6e rad\n", PL->dphi);
			//printf("Turn %d\n", i);
			Slice->track(start, end);

#pragma omp barrier
			long_tracker->track(start, end);
#pragma omp barrier
#pragma omp single
			{
				RfP->counter++;
			}
			/*
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

			 */

		}
	}

//printf("Total simulation time: %.10lf\n", long_tracker->elapsed_time);
//printf("Time/turn : %.10lf\n", long_tracker->elapsed_time / N_t);
//ftype result = mymath::trapezoid(Slice->n_macroparticles,
//Slice->bin_centers, Slice->n_slices);
//printf("result = %e\n", result);

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

//delete_array(ps);
	v.clear();
	v2.clear();
	delete PL;
	delete Slice;
	delete long_tracker;
	delete RfP;
	delete GP;
	delete Beam;

	printf("Done!\n");

}


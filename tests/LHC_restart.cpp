/*
 * LHC_restart.cpp
 *
 *  Created on: Apr 12, 2016
 *      Author: kiliakis
 */

#include <blond/globals.h>
#include <blond/utilities.h>
#include <blond/math_functions.h>
#include <omp.h>
#include <stdio.h>
#include <blond/input_parameters/GeneralParameters.h>
#include <blond/input_parameters/RfParameters.h>
#include <blond/beams/Beams.h>
#include <blond/beams/Slices.h>
#include <blond/trackers/Tracker.h>
#include <blond/llrf/PhaseLoop.h>
#include <blond/optionparser.h>
#include <algorithm>
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
   "../tests/input_files/LHC_restart/";

// Global variables
GeneralParameters *GP;
Beams *Beam;
Slices *Slice;
RfParameters *RfP;
int n_threads = 1;
//const int size = 14e6;
const int from_line = 0;

void parse_args(int argc, char **argv);

// Simulation setup -------------------------------------------------------------
int main(int argc, char **argv)
{

   parse_args(argc, argv);
   // Environmental variables
   /*
   N_t = atoi(util::GETENV("N_TURNS")) ? atoi(util::GETENV("N_TURNS")) : N_t;
   N_p = atoi(util::GETENV("N_PARTICLES")) ? atoi(util::GETENV("N_PARTICLES")) : N_p;
   N_slices = atoi(util::GETENV("N_SLICES")) ? atoi(util::GETENV("N_SLICES")) : N_slices;
   n_threads =
       atoi(util::GETENV("N_THREADS")) ? atoi(util::GETENV("N_THREADS")) : n_threads;
   */
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

   //printf("Setting up the simulation..\n");
   timespec begin, end;
   util::get_time(begin);

   std::vector < ftype > v;
   util::read_vector_from_file(v, datafiles + "LHC_momentum_programme");

   // optional
   v.erase(v.begin(), v.begin() + from_line);

   //std::cout << "vector size is " << v.size() << "\n";
   int remaining = N_t + 1 - v.size();
   for (int i = 0; i < remaining; ++i) {
      v.push_back(6.5e12);
   }
   assert((int) v.size() == N_t + 1);
   ftype *ps = &v[0];   //new ftype[v.size()];

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
   util::read_vector_from_file(v2, datafiles + "coords_13000001.dat");
   int k = 0;
   for (unsigned int i = 0; i < v2.size(); i += 3) {
      Beam->dt[k] = v2[i] * 1e-9; // [s]
      Beam->dE[k] = v2[i + 1] * 1e6; // [eV]
      Beam->id[k] = v2[i + 2];
      k++;
   }

   Slice = new Slices(N_slices, 0, -0.5e-9, 3e-9);
   printf("Beam generated, slices set...\n");
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
      //    threads, tile, start, end);
      for (int i = 0; i < N_t; ++i) {

         #pragma omp single
         {
            printf("\nTurn %d\n", i);

            if (RfP->counter < 570000)
               PL->reference = 0.5236;
            else
               PL->reference = 1.0472;

         }

         if (id == 0) util::get_time(begin_t);

         Slice->track(start, end);

         #pragma omp barrier
         if (id == 0) slice_time += util::time_elapsed(begin_t);
         if (id == 0) util::get_time(begin_t);

         long_tracker->track(start, end);

         #pragma omp barrier
         if (id == 0) track_time += util::time_elapsed(begin_t);

         #pragma omp single
         {
            //printf("   Beam energy %.6e eV\n", GP->energy[0]);
            printf("   RF phase %.6e rad\n", RfP->dphi_RF[0]);
            printf("   PL phase correction %.6e rad\n", PL->dphi);
            RfP->counter++;
         }

      }
   }

   util::get_time(end);
   util::print_time("Simulation Time", begin, end);

   double total_time = track_time + slice_time;
   printf("Track time : %.4lf ( %.2lf %% )\n", track_time,
          100 * track_time / total_time);
   printf("Slice time : %.4lf ( %.2lf %% )\n", slice_time,
          100 * slice_time / total_time);
   util::dump(Beam->dE, 10, "dE\n");
   util::dump(Beam->dt, 10, "dt\n");
   util::dump(Slice->n_macroparticles, 10, "n_macroparticles\n");

   delete PL;
   delete Slice;
   delete long_tracker;
   delete RfP;
   delete GP;
   delete Beam;

   printf("Done!\n");

}

void parse_args(int argc, char **argv)
{
   using namespace std;
   using namespace option;

   enum optionIndex {UNKNOWN, HELP, N_THREADS, N_TURNS,
                     N_PARTICLES, N_SLICES, OPTIONS_NUM
                    };

   const option::Descriptor usage[] = {
      {
         UNKNOWN, 0, "", "", Arg::None,                  "USAGE: ./LHC_restart [options]\n\n"
         "Options:"
      },
      {  HELP, 0, "h", "help", Arg::None,                "  --help,              -h        Print usage and exit." },
      {N_TURNS, 0, "t", "turns", util::Arg::Numeric,        "  --turns=<num>,       -t <num>  Number of turns (default: 10k)" },
      {N_PARTICLES, 0, "p", "particles", util::Arg::Numeric,   "  --particles=<num>,   -p <num>  Number of particles (default: 10k)" },
      {N_SLICES, 0, "s", "slices", util::Arg::Numeric,      "  --slices=<num>,      -s <num>  Number of slices (default: 100)" },
      {N_THREADS, 0, "m", "threads", util::Arg::Numeric,       "  --threads=<num>,     -m <num>  Number of threads (default: 1)" },
      {
         UNKNOWN, 0, "", "", Arg::None,                  "\nExamples:\n"
         "\t./LHC_restart\n"
         "\t./LHC_restart -t 1000 -p 10000 -m 4\n"
      },
      {0, 0, 0, 0, 0, 0}
   };

   argc -= (argc > 0);
   argv += (argc > 0); // skip program name argv[0] if present
   Stats stats(usage, argc, argv);
   vector<Option> options(stats.options_max);
   vector<Option> buffer(stats.buffer_max);
   Parser parse(usage, argc, argv, &options[0], &buffer[0]);

   if (options[HELP]) {
      printUsage(cout, usage);
      exit(0);
   }

   for (int i = 0; i < parse.optionsCount(); ++i) {
      Option &opt = buffer[i];
      //fprintf(stdout, "Argument #%d is ", i);
      switch (opt.index()) {
         case HELP:
         // not possible, because handled further above and exits the program
         case N_TURNS:
            N_t = atoi(opt.arg);
            //fprintf(stdout, "--numeric with argument '%s'\n", opt.arg);
            break;
         case N_THREADS:
            n_threads = atoi(opt.arg);
            //fprintf(stdout, "--numeric with argument '%s'\n", opt.arg);
            break;
         case N_SLICES:
            N_slices = atoi(opt.arg);
            //fprintf(stdout, "--numeric with argument '%s'\n", opt.arg);
            break;
         case N_PARTICLES:
            N_p = atoi(opt.arg);
            //fprintf(stdout, "--numeric with argument '%s'\n", opt.arg);
            break;
         case UNKNOWN:
            // not possible because Arg::Unknown returns ARG_ILLEGAL
            // which aborts the parse with an error
            break;
      }
   }


}


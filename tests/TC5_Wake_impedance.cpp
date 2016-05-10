/*
 * TC5_Wake_impedance.cpp
 *
 *  Created on: May 9, 2016
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
#include "../impedances/Intensity.h"

const std::string datafiles =
   "../tests/input_files/TC5_Wake_impedance/";

// Simulation parameters --------------------------------------------------------
// Bunch parameters
const int N_b = (int) 1e10;                          // Intensity
const ftype tau_0 = 2e-9;                       // Initial bunch length, 4 sigma [s]
// const particle_type particle = proton;
// Machine and RF parameters
const ftype C = 6911.56;                        // Machine circumference [m]
const ftype p_i = 25.92e9;                      // Synchronous momentum [eV/c]
//const ftype p_f = 460.005e9;                  // Synchronous momentum, final
const long h = 4620;                            // Harmonic number
const ftype V = 0.9e6;                          // RF voltage [V]
const ftype dphi = 0;                           // Phase modulation/offset
const ftype gamma_t = 1 / std::sqrt(0.00192);   // Transition gamma
const ftype alpha = 1.0 / gamma_t / gamma_t;    // First order mom. comp. factor
const int alpha_order = 1;
const int n_sections = 1;
// Tracking details

int N_t = 2;    // Number of turns to track
int N_p = 5000000;         // Macro-particles

int n_threads = 1;
int N_slices = 1 << 8; // = (2^8)

GeneralParameters *GP;
Beams *Beam;
Slices *Slice;
RfParameters *RfP;


void parse_args(int argc, char **argv);

// Simulation setup -------------------------------------------------------------
int main(int argc, char **argv)
{

   parse_args(argc, argv);

   omp_set_num_threads(n_threads);

   /// initializations

   printf("Setting up the simulation...\n\n");
   printf("Number of turns: %d\n", N_t);
   printf("Number of macro-particles: %d\n", N_p);
   printf("Number of Slices: %d\n", N_slices);

   #pragma omp parallel
   {
      if (omp_get_thread_num() == 0)
         printf("Number of openmp threads: %d\n", omp_get_num_threads());
   }

   timespec begin, end;
   util::get_time(begin);

   ftype *momentum = new ftype[N_t + 1];
   std::fill_n(momentum, N_t + 1, p_i);

   ftype *alpha_array = new ftype[(alpha_order + 1) * n_sections];
   std::fill_n(alpha_array, (alpha_order + 1) * n_sections, alpha);

   ftype *C_array = new ftype[n_sections];
   std::fill_n(C_array, n_sections, C);

   ftype *h_array = new ftype[n_sections * (N_t + 1)];
   std::fill_n(h_array, (N_t + 1) * n_sections, h);

   ftype *V_array = new ftype[n_sections * (N_t + 1)];
   std::fill_n(V_array, (N_t + 1) * n_sections, V);

   ftype *dphi_array = new ftype[n_sections * (N_t + 1)];
   std::fill_n(dphi_array, (N_t + 1) * n_sections, dphi);

// TODO variables must be in the correct format (arrays for all)

   GP = new GeneralParameters(N_t, C_array, alpha_array, alpha_order, momentum,
                              proton);

   Beam = new Beams(N_p, N_b);

   RfP = new RfParameters(n_sections, h_array, V_array, dphi_array);
   //printf("omega_rf = %lf\n",RfP->omega_RF[0]);

   RingAndRfSection *long_tracker = new RingAndRfSection();
   //util::dump(Beam->dE, 10, "dE\n");

   longitudinal_bigaussian(tau_0 / 4, 0, 1, false);

   Slice = new Slices(N_slices, 0, 0, 2 * constant::pi, rad);


   std::vector<ftype> v;
   util::read_vector_from_file(v, datafiles +
                               "TC5_new_HQ_table.dat");
   assert(v.size() % 3 == 0);
   //ftype *p = &v[0];
   //util::dump(p, 30, "input array\n");
   std::vector<ftype> R_shunt, f_res, Q_factor;

   R_shunt.reserve(v.size() / 3);
   f_res.reserve(v.size() / 3);
   Q_factor.reserve(v.size() / 3);

   for (uint i = 0; i < v.size(); i += 3) {
      f_res.push_back(v[i] * 1e9);
      Q_factor.push_back(v[i + 1]);
      R_shunt.push_back(v[i + 2] * 1e6);
   }

   Resonators *resonator = new Resonators(R_shunt, f_res, Q_factor);

   //ftype *p = &R_shunt[0];
   //util::dump(p, R_shunt.size(), "R_shunt\n");

   /*
   double slice_time = 0, track_time = 0;

   timespec begin_t;

   
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

         if (id == 0) util::get_time(begin_t);

         long_tracker->track(start, end);

         #pragma omp barrier

         if (id == 0) track_time += util::time_elapsed(begin_t);
         if (id == 0) util::get_time(begin_t);

         Slice->track(start, end);

         #pragma omp barrier

         if (id == 0) slice_time += util::time_elapsed(begin_t);

         #pragma omp single
         {
            Slice->fwhm();

            if (i % 1000 == 0) {
               util::dump(Slice->bl_fwhm, "bl_fwhm");
               util::dump(Slice->bp_fwhm, "bp_fwhm");
            }

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
   */

   delete Slice;
   delete long_tracker;
   delete RfP;
   delete GP;
   delete Beam;
   delete resonator;
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
         UNKNOWN, 0, "", "", Arg::None,                  "USAGE: ./TC5_Wake_impedance [options]\n\n"
         "Options:"
      },
      {  HELP, 0, "h", "help", Arg::None,                "  --help,              -h        Print usage and exit." },
      {N_TURNS, 0, "t", "turns", util::Arg::Numeric,        "  --turns=<num>,       -t <num>  Number of turns (default: 10k)" },
      {N_PARTICLES, 0, "p", "particles", util::Arg::Numeric,   "  --particles=<num>,   -p <num>  Number of particles (default: 10k)" },
      {N_SLICES, 0, "s", "slices", util::Arg::Numeric,      "  --slices=<num>,      -s <num>  Number of slices (default: 100)" },
      {N_THREADS, 0, "m", "threads", util::Arg::Numeric,       "  --threads=<num>,     -m <num>  Number of threads (default: 1)" },
      {
         UNKNOWN, 0, "", "", Arg::None,                  "\nExamples:\n"
         "\t./TC5_Wake_impedance\n"
         "\t./TC5_Wake_impedance -t 1000 -p 10000 -m 4\n"
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


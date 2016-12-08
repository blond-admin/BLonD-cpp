/*
 * TC1_Acceleration.cpp
 *
 *  Created on: Mar 9, 2016
 *      Author: kiliakis
 */
#include <blond/beams/Beams.h>
#include <blond/beams/Distributions.h>
#include <blond/beams/Slices.h>
#include <blond/globals.h>
#include <blond/input_parameters/GeneralParameters.h>
#include <blond/input_parameters/RfParameters.h>
#include <blond/math_functions.h>
#include <blond/optionparser.h>
#include <blond/trackers/Tracker.h>
#include <blond/utilities.h>
#include <stdio.h>



using namespace std;
// Simulation parameters
// --------------------------------------------------------
// Bunch parameters
const long long N_b = 1e9;  // Intensity
const double tau_0 = 0.4e-9; // Initial bunch length, 4 sigma [s]

// Machine and RF parameters
const double C = 26658.883;                   // Machine circumference [m]
const long long p_i = 450e9;                 // Synchronous momentum [eV/c]
const double p_f = 460.005e9;                 // Synchronous momentum, final
const long long h = 35640;                   // Harmonic number
const double V = 6e6;                         // RF voltage [V]
const double dphi = 0;                        // Phase modulation/offset
const double gamma_t = 55.759505;             // Transition gamma
const double alpha = 1.0 / gamma_t / gamma_t; // First order mom. comp. factor
const int alpha_order = 1;
const int n_sections = 1;
// Tracking details

int N_t = 10000; // Number of turns to track
int N_p = 10000; // Macro-particles

int N_slices = 100;

void parse_args(int argc, char **argv);


// Simulation setup
// -------------------------------------------------------------
int main(int argc, char **argv)
{
    parse_args(argc, argv);

    omp_set_num_threads(Context::n_threads);

    cout << "Setting up the simulation...\n\n";
    cout << "Number of turns: " <<  N_t << "\n";
    cout << "Number of macro-particles: " << N_p << "\n";
    cout << "Number of Slices: " <<  N_slices << "\n";
    cout << "Number of openmp threads: " <<  Context::n_threads << "\n";

    // timespec begin, end;
    // util::get_time(begin);

    f_vector_2d_t momentumVec(n_sections, f_vector_t(N_t + 1));
    for (auto &v : momentumVec)
        mymath::linspace(v.data(), p_i, p_f, N_t + 1);

    f_vector_2d_t alphaVec(n_sections, f_vector_t(alpha_order, alpha));

    f_vector_t CVec(n_sections, C);

    f_vector_2d_t hVec(n_sections, f_vector_t(N_t + 1, h));

    f_vector_2d_t voltageVec(n_sections, f_vector_t(N_t + 1, V));

    f_vector_2d_t dphiVec(n_sections, f_vector_t(N_t + 1, dphi));

    auto GP = Context::GP = new GeneralParameters(N_t, CVec, alphaVec,
             momentumVec,
            GeneralParameters::particle_t::proton);

    // auto Beam = Context::Beam = new Beams(N_p, N_b);
    auto Beam = Context::Beam = new Beams(GP, N_p, N_b);

    auto RfP = Context::RfP = new RfParameters(GP, n_sections, hVec,
            voltageVec, dphiVec);


    auto long_tracker = RingAndRfSection();

    longitudinal_bigaussian(GP, RfP, Beam, tau_0 / 4, 0, 1, false);

    auto Slice = Context::Slice = new Slices(RfP, Beam, N_slices);

    double slice_time = 0, track_time = 0;
    timespec begin_t;

    for (int i = 0; i < N_t; ++i) {
        // if (i % 10000 == 0)
        //     cout << "Turn number " << i << " of " << N_t << "\n";
        // util::get_time(begin_t);
        // long_tracker->track();
        // track_time += util::time_elapsed(begin_t);

        util::get_time(begin_t);
        long_tracker.track();
        track_time += util::time_elapsed(begin_t);

        util::get_time(begin_t);
        Slice->track();
        slice_time += util::time_elapsed(begin_t);
        // Slice->fwhm();

        // if (i % 1000 == 0) {
        //    util::dump(Slice->bl_fwhm, "bl_fwhm");
        //    util::dump(Slice->bp_fwhm, "bp_fwhm");
        // }

        // RfP->counter++;
    }

    cout << scientific;
    cout << "Average Turn Time : " << (slice_time + track_time) / N_t
         << "\n";
    cout << "Average Tracker Track Time : " << track_time / N_t << "\n";
    cout << "Average Slice Track Time : " << slice_time / N_t << "\n";

    // util::dump(Beam->dE.data(), 10, "dE\n");
    // util::dump(Beam->dt.data(), 10, "dt\n");
    // util::dump(Slice->n_macroparticles, 10, "n_macroparticles\n");
    delete Slice;
    delete RfP;
    delete GP;
    delete Beam;

    printf("Done!\n");
}

void parse_args(int argc, char **argv)
{
    using namespace std;
    using namespace option;

    enum optionIndex {
        UNKNOWN,
        HELP,
        N_THREADS,
        N_TURNS,
        N_PARTICLES,
        N_SLICES,
        OPTIONS_NUM
    };

    const option::Descriptor usage[] = {
        {
            UNKNOWN, 0, "", "", Arg::None,
            "USAGE: ./TC1_Acceleration [options]\n\n"
            "Options:"
        },
        {
            HELP, 0, "h", "help", Arg::None,
            "  --help,              -h        Print usage and exit."
        },
        {
            N_TURNS, 0, "t", "turns", util::Arg::Numeric,
            "  --turns=<num>,       -t <num>  Number of turns (default: 10k)"
        },
        {
            N_PARTICLES, 0, "p", "particles", util::Arg::Numeric,
            "  --particles=<num>,   -p <num>  Number of particles (default: 10k)"
        },
        {
            N_SLICES, 0, "s", "slices", util::Arg::Numeric,
            "  --slices=<num>,      -s <num>  Number of slices (default: 100)"
        },
        {
            N_THREADS, 0, "m", "threads", util::Arg::Numeric,
            "  --threads=<num>,     -m <num>  Number of threads (default: 1)"
        },
        {
            UNKNOWN, 0, "", "", Arg::None,
            "\nExamples:\n"
            "\t./TC1_Acceleration\n"
            "\t./TC1_Acceleration -t 1000 -p 10000 -m 4\n"
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
        // fprintf(stdout, "Argument #%d is ", i);
        switch (opt.index()) {
            case HELP:
            // not possible, because handled further above and exits the program
            case N_TURNS:
                N_t = atoi(opt.arg);
                // fprintf(stdout, "--numeric with argument '%s'\n", opt.arg);
                break;
            case N_THREADS:
                Context::n_threads = atoi(opt.arg);
                // fprintf(stdout, "--numeric with argument '%s'\n", opt.arg);
                break;
            case N_SLICES:
                N_slices = atoi(opt.arg);
                // fprintf(stdout, "--numeric with argument '%s'\n", opt.arg);
                break;
            case N_PARTICLES:
                N_p = atoi(opt.arg);
                // fprintf(stdout, "--numeric with argument '%s'\n", opt.arg);
                break;
            case UNKNOWN:
                // not possible because Arg::Unknown returns ARG_ILLEGAL
                // which aborts the parse with an error
                break;
        }
    }
}

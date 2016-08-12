/*
 * TC8_Phase_Loop.cpp
 *
 *  Created on: July 11, 2016
 *      Author: kiliakis
 */
#include <blond/beams/Distributions.h>
#include <blond/globals.h>
#include <blond/llrf/PhaseLoop.h>
#include <blond/math_functions.h>
#include <blond/trackers/Tracker.h>
#include <blond/utilities.h>

// Simulation parameters
// --------------------------------------------------------

// Bunch parameters
const uint N_b = 0; // Intensity

// Machine and RF parameters
const ftype radius = 25;
const ftype C = 2 * constant::pi * radius;   // Machine circumference [m]
const ftype p_i = 310891054.809;             // Synchronous momentum [eV/c]
const uint h = 1;                            // Harmonic number
const ftype V = 8000;                        // RF voltage [V]
const ftype dphi = -constant::pi;            // Phase modulation/offset
const ftype gamma_t = 4.076750841;           // Transition gamma
const ftype alpha = 1.0 / gamma_t / gamma_t; // First order mom. comp. factor
const uint alpha_order = 1;
const uint n_sections = 1;
// Tracking details

uint N_t = 500;    // Number of turns to track
uint N_p = 100000; // Macro-particles

uint N_slices = 200; // = (2^8)

void parse_args(int argc, char** argv);

// Simulation setup
// -------------------------------------------------------------
int main(int argc, char** argv) {

    parse_args(argc, argv);

    omp_set_num_threads(Context::n_threads);

    /// initializations

    printf("Setting up the simulation...\n\n");
    printf("Number of turns: %d\n", N_t);
    printf("Number of macro-particles: %d\n", N_p);
    printf("Number of Slices: %d\n", N_slices);
    printf("Number of openmp threads: %d\n", Context::n_threads);

    f_vector_2d_t momentumVec(n_sections, f_vector_t(N_t + 1, p_i));

    f_vector_2d_t alphaVec(n_sections, f_vector_t(alpha_order + 1, alpha));

    f_vector_t CVec(n_sections, C);

    f_vector_2d_t hVec(n_sections, f_vector_t(N_t + 1, h));

    f_vector_2d_t voltageVec(N_t + 1, f_vector_t(n_sections, V));

    f_vector_2d_t dphiVec(N_t + 1, f_vector_t(n_sections, dphi));

    Context::GP = new GeneralParameters(N_t, CVec, alphaVec, alpha_order,
                                        momentumVec, proton);

    Context::Beam = new Beams(N_p, N_b);

    Context::RfP = new RfParameters(n_sections, hVec, voltageVec, dphiVec);

    // long_tracker = new RingAndRfSection(RfP);

    Context::Slice = new Slices(N_slices, 0, -constant::pi, constant::pi,
                                cuts_unit_type::rad);

    longitudinal_bigaussian(200e-9, 1e6, 1, false);

    auto psb =
        new PSB(f_vector_t(N_t, 1.0 / 25e-6), f_vector_t{0, 0}, 10e-6, 7);

    auto long_tracker =
        new RingAndRfSection(Context::RfP, simple, NULL, NULL, false);

    for (auto& v : Context::Beam->dE)
        v += 90.0e3;

    Context::Slice->track();

    timespec begin;
    double track_time = 0.0;
    double slice_time = 0.0;
    double pl_time = 0.0;

    for (uint i = 0; i < N_t; ++i) {

        util::get_time(begin);
        psb->track();
        pl_time += util::time_elapsed(begin);

        util::get_time(begin);
        long_tracker->track();
        track_time += util::time_elapsed(begin);

        util::get_time(begin);
        Context::Slice->track();
        slice_time += util::time_elapsed(begin);
    }

    std::cout << std::scientific;
    std::cout << "Average Turn Time : "
              << (track_time + slice_time + pl_time) / N_t << std::endl;
    std::cout << "Average Tracker Track Time : " << track_time / N_t
              << std::endl;
    std::cout << "Average PhaseLoop Time : " << pl_time / N_t << std::endl;
    std::cout << "Average Slice Track Time : " << slice_time / N_t << std::endl;

    // util::dump(Beam->dE, "dE\n", 10);
    // util::dump(Beam->dt, "dt\n", 10);
    // util::dump(Slice->n_macroparticles, "n_macroparticles\n", 10);

    delete Context::Slice;
    delete long_tracker;
    delete Context::RfP;
    delete Context::GP;
    delete Context::Beam;
    delete psb;

    printf("Done!\n");
}

void parse_args(int argc, char** argv) {
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
        {UNKNOWN, 0, "", "", Arg::None, "USAGE: ./TC8_Phase_Loop [options]\n\n"
                                        "Options:"},
        {HELP, 0, "h", "help", Arg::None,
         "  --help,              -h        Print usage and exit."},
        {N_TURNS, 0, "t", "turns", util::Arg::Numeric,
         "  --turns=<num>,       -t <num>  Number of turns (default: 500)"},
        {N_PARTICLES, 0, "p", "particles", util::Arg::Numeric,
         "  --particles=<num>,   -p <num>  Number of particles (default: "
         "100k)"},
        {N_SLICES, 0, "s", "slices", util::Arg::Numeric,
         "  --slices=<num>,      -s <num>  Number of slices (default: 200)"},
        {N_THREADS, 0, "m", "threads", util::Arg::Numeric,
         "  --threads=<num>,     -m <num>  Number of threads (default: 1)"},
        {UNKNOWN, 0, "", "", Arg::None,
         "\nExamples:\n"
         "\t./TC8_Phase_Loop\n"
         "\t./TC8_Phase_Loop -t 1000 -p 10000 -m 4\n"},
        {0, 0, 0, 0, 0, 0}};

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
        Option& opt = buffer[i];
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

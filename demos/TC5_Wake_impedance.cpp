/*
 * TC5_Wake_impedance.cpp
 *
 *  Created on: May 9, 2016
 *      Author: kiliakis
 */
#include <blond/beams/Beams.h>
#include <blond/beams/Distributions.h>
#include <blond/beams/Slices.h>
#include <blond/globals.h>
#include <blond/input_parameters/GeneralParameters.h>
#include <blond/input_parameters/RfParameters.h>
#include <blond/math_functions.h>
#include <blond/trackers/Tracker.h>
#include <blond/utilities.h>
#include <stdio.h>
#include <blond/impedances/InducedVoltage.h>

using namespace std;

const string datafiles = DEMO_FILES "/TC5_Wake_impedance/";

// Simulation parameters
// --------------------------------------------------------
// Bunch parameters
const long long N_b = 1e10; // Intensity
const double tau_0 = 2e-9;                 // Initial bunch length, 4 sigma [s]
// const particle_type particle = proton;
// Machine and RF parameters
const double C = 6911.56;   // Machine circumference [m]
const double p_i = 25.92e9; // Synchronous momentum [eV/c]
// const double p_f = 460.005e9;                  // Synchronous momentum, final
const double h = 4620;                     // Harmonic number
const double V = 0.9e6;                        // RF voltage [V]
const double dphi = 0;                         // Phase modulation/offset
const double gamma_t = 1 / sqrt(0.00192); // Transition gamma
const double alpha = 1.0 / gamma_t / gamma_t;  // First order mom. comp. factor
const int alpha_order = 1;
const int n_sections = 1;
// Tracking details

int N_t = 1000;    // Number of turns to track
int N_p = 5000000; // Macro-particles

int N_slices = 1 << 8; // = (2^8)

void parse_args(int argc, char **argv);

// Simulation setup
// -------------------------------------------------------------
int main(int argc, char **argv)
{
    parse_args(argc, argv);

    omp_set_num_threads(Context::n_threads);

    /// initializations

    printf("Setting up the simulation...\n\n");
    printf("Number of turns: %d\n", N_t);
    printf("Number of macro-particles: %d\n", N_p);
    printf("Number of Slices: %d\n", N_slices);
    printf("Number of openmp threads: %d\n", Context::n_threads);

    timespec begin;
    // timespec end;
    f_vector_2d_t momentumVec(n_sections, f_vector_t(N_t + 1, p_i));

    f_vector_2d_t alphaVec(n_sections, f_vector_t(alpha_order + 1, alpha));

    f_vector_t CVec(n_sections, C);

    f_vector_2d_t hVec(n_sections, f_vector_t(N_t + 1, h));

    f_vector_2d_t voltageVec(n_sections, f_vector_t(N_t + 1, V));

    f_vector_2d_t dphiVec(n_sections, f_vector_t(N_t + 1, dphi));

    Context::GP = new GeneralParameters(N_t, CVec, alphaVec,
                                        alpha_order, momentumVec,
                                        GeneralParameters::particle_t::proton);

    // Context::Beam = new Beams(N_p, N_b);

    auto GP = Context::GP;
    auto Beam = Context::Beam = new Beams(GP, N_p, N_b);

    auto RfP = Context::RfP = new RfParameters(GP, n_sections, hVec,
            voltageVec, dphiVec);


    RingAndRfSection *long_tracker = new RingAndRfSection();

    longitudinal_bigaussian(GP, RfP, Beam, tau_0 / 4, 0, 1, false);

    auto Slice = Context::Slice = new Slices(RfP, Beam, N_slices,
            0, 0, 2 * constant::pi, Slices::cuts_unit_t::rad);
    // util::dump(Slice->bin_centers, 10, "bin_centers\n");

    f_vector_t v;
    util::read_vector_from_file(v, datafiles + "TC5_new_HQ_table.dat");
    assert(v.size() % 3 == 0);

    f_vector_t R_shunt, f_res, Q_factor;

    R_shunt.reserve(v.size() / 3);
    f_res.reserve(v.size() / 3);
    Q_factor.reserve(v.size() / 3);

    for (uint i = 0; i < v.size(); i += 3) {
        f_res.push_back(v[i] * 1e9);
        Q_factor.push_back(v[i + 1]);
        R_shunt.push_back(v[i + 2] * 1e6);
    }

    auto resonator = new Resonators(R_shunt, f_res, Q_factor);

    auto indVoltTime = new InducedVoltageTime(Slice, {resonator},
            InducedVoltageTime::freq_domain);

    auto totVol = new TotalInducedVoltage(Beam, Slice, {indVoltTime});

    auto indTrack = 0.0, longTrack = 0.0, sliceTrack = 0.0;
    for (int i = 0; i < N_t; ++i) {

        util::get_time(begin);
        totVol->track(Beam);
        // util::print_time_elapsed("Induced Voltage Track", begin);
        indTrack += util::time_elapsed(begin);

        util::get_time(begin);
        long_tracker->track();
        // util::print_time_elapsed("Tracker Track", begin);
        longTrack += util::time_elapsed(begin);

        util::get_time(begin);
        Slice->track();
        sliceTrack += util::time_elapsed(begin);

        // util::print_time_elapsed("Slice Track", begin);
    }

    cout << scientific;
    cout << "Average Turn Time : "
         << (indTrack + longTrack + sliceTrack) / N_t << endl;
    cout << "Average Induced Voltage Track Time : " << indTrack / N_t
         << endl;
    cout << "Average Tracker Track Time : " << longTrack / N_t
         << endl;
    cout << "Average Slice Track Time : " << sliceTrack / N_t << endl;

    // util::dump(Beam->dE, 10, "dE\n");
    // util::dump(Beam->dt, 10, "dt\n");
    // util::dump(Slice->n_macroparticles, 10, "n_macroparticles\n");

    delete Slice;
    delete long_tracker;
    delete RfP;
    delete GP;
    delete Beam;
    delete resonator;
    delete indVoltTime;
    delete totVol;

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
            "USAGE: ./TC5_Wake_impedance [options]\n\n"
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

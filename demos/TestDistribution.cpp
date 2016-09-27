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
#include <blond/plots/plot_beams.h>
#include <blond/plots/plot_parameters.h>
#include <blond/plots/plot_slices.h>
#include <blond/plots/plot_llrf.h>
#include <blond/python.h>
#include <blond/monitors/Monitors.h>
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

int N_t = 500;    // Number of turns to track
int N_p = 100000; // Macro-particles

int N_slices = 200; // = (2^8)

void parse_args(int argc, char **argv);

// Simulation setup
// -------------------------------------------------------------
int main(int argc, char **argv)
{

    parse_args(argc, argv);
    python::initialize();
    omp_set_num_threads(Context::n_threads);

    /// initializations

    printf("Setting up the simulation...\n\n");
    printf("Number of turns: %d\n", N_t);
    printf("Number of macro-particles: %d\n", N_p);
    printf("Number of Slices: %d\n", N_slices);
    printf("Number of openmp threads: %d\n", Context::n_threads);

    f_vector_2d_t momentumVec(n_sections, f_vector_t(N_t + 1, p_i));

    // f_vector_2d_t momentumVec(n_sections, f_vector_t(N_t + 1));
    // for (auto& v : momentumVec)
    //     mymath::linspace(v.data(), p_i, p_i*1.001, N_t + 1);

    f_vector_2d_t alphaVec(n_sections, f_vector_t(alpha_order + 1, alpha));

    f_vector_t CVec(n_sections, C);

    f_vector_2d_t hVec(n_sections, f_vector_t(N_t + 1, h));

    f_vector_2d_t voltageVec(n_sections, f_vector_t(N_t + 1, V));

    f_vector_2d_t dphiVec(n_sections, f_vector_t(N_t + 1, dphi));

    auto GP = Context::GP = new GeneralParameters(N_t, CVec, alphaVec, alpha_order,
            momentumVec, proton);

    auto Beam = Context::Beam = new Beams(N_p, N_b);


    auto RfP = Context::RfP = new RfParameters(n_sections, hVec, voltageVec, dphiVec);

    // f_vector_t time(N_t + 1);
    // for (int i = 0; i < time.size(); i++)
    //     time[i] = 1.0 * (i + 1);
    // plot_voltage_programme(time, voltageVec[0]);

    // auto bunchmonitor = BunchMonitor(GP, RfP, Beam, "bunch.h5", 10);
    auto tracker = RingAndRfSection();

    longitudinal_bigaussian(200e-9, 1e6, 1, false);

    auto slice = Context::Slice = new Slices(N_slices, 0, -constant::pi, constant::pi,
            cuts_unit_type::rad);


    // for (int i = 0; i < N_t; i++) {
    //     tracker.track();
    //     slice->track();
    //     bunchmonitor.track();
    // }
    // bunchmonitor.track();
    // bunchmonitor.close();
    plot_PL_bunch_phase(RfP, "output_data_full.h5");
    plot_PL_RF_phase(RfP, "output_data_full.h5");
    plot_PL_phase_corr(RfP, "output_data_full.h5");
    plot_PL_RF_freq(RfP, "output_data_full.h5");
    plot_PL_freq_corr(RfP, "output_data_full.h5");
    plot_RF_phase_error(RfP, "output_data_full.h5");
    plot_COM_motion(RfP, "output_data_full.h5");
    plot_LHCNoiseFB(RfP, "output_data_full.h5");
    plot_LHCNoiseFB_FWHM(RfP, "output_data_full.h5");
    plot_LHCNoiseFB_FWHM_bbb(RfP, "output_data_full.h5");
    // auto psb =
    //     new PSB(f_vector_t(N_t, 1.0 / 25e-6), f_vector_t{0, 0}, 10e-6, 7);

    // auto long_tracker = new RingAndRfSection(RfP, simple);

    // std::vector<RingAndRfSection *> trackerVector{long_tracker};
    // auto full_ring = new FullRingAndRf(trackerVector);

    // std::map<std::string, std::string> line_density_opt;
    // line_density_opt["type"] = "gaussian";
    // line_density_opt["bunch_length"] = "200e-9";
    // line_density_opt["density_variable"] = "density_from_J";

    // // longitudinal_bigaussian(200e-9, 1e6, 1, false);

    // matched_from_line_density(full_ring, line_density_opt, "lowest_freq", "savefig");
    // util::dump(Beam->dt, "Beam->dt\n");
    // util::dump(Beam->dE, "Beam->dE\n");
    // util::dump(Beam->id, "Beam->id\n");

    // plot_long_phase_space(GP, RfP, Beam, 5e-7, 1.2e-6, -3e5, 3e5, "s", 1, true);
    // Context::Slice->track();
    // Context::Slice->beam_spectrum_generation(100, false);

    // plot_bunch_length_evol(RfP, "bunch.h5");
    // plot_position_evol(RfP, "bunch.h5");
    // plot_energy_evol(RfP, "bunch.h5");
    // plot_transmitted_particles(RfP, "bunch.h5");

    // auto f = mymath::arange(0., 5.6227612455e+03, 1.12455000e-02);
    // auto spectrum = f_vector_t(4980, 1.111e-07);//np.concatenate((1.11100000e-07 * np.ones(4980), np.zeros(495021)))
    // spectrum.resize(495021 + 4980, 0);

    // // auto rfnoise = PhaseNoise();

    // // rfnoise.spectrum_to_phase_noise();
    // plot_noise_spectrum(f, spectrum, 100);
    // plot_phase_noise(f, spectrum, 100);
    // plot_bunch_length_evol_gaussian(RfP, "bunch.h5", 1, "./");
    // plot_beam_profile(Context::Slice, 0);
    // plot_beam_spectrum(Context::Slice, 0);
    // plot_beam_profile_derivative(Context::Slice, 100, "-",
    //                              "fig", {"gradient", "diff", "filter1d"});
    // std::map<std::string, std::string> distribution_opt;
    // distribution_opt["type"] = "binomial";
    // distribution_opt["exponent"] = "1.5";
    // distribution_opt["bunch_length"] = "100.0e-9";
    // distribution_opt["density_variable"] = "density_from_J";

    // matched_from_distribution_density(full_ring, distribution_opt);



    // util::dump(Beam->dt, "Beam->dt\n");
    // util::dump(Beam->dE, "Beam->dE\n");

    // for (auto& v : Context::Beam->dE)
    //     v += 90.0e3;

    // Context::Slice->track();

    // timespec begin;
    // double track_time = 0.0;
    // double slice_time = 0.0;
    // double pl_time = 0.0;

    // for (uint i = 0; i < N_t; ++i) {

    //     util::get_time(begin);
    //     psb->track();
    //     pl_time += util::time_elapsed(begin);

    //     util::get_time(begin);
    //     long_tracker->track();
    //     track_time += util::time_elapsed(begin);

    //     util::get_time(begin);
    //     Context::Slice->track();
    //     slice_time += util::time_elapsed(begin);
    // }

    // std::cout << std::scientific;
    // std::cout << "Average Turn Time : "
    //           << (track_time + slice_time + pl_time) / N_t << std::endl;
    // std::cout << "Average Tracker Track Time : " << track_time / N_t
    //           << std::endl;
    // std::cout << "Average PhaseLoop Time : " << pl_time / N_t << std::endl;
    // std::cout << "Average Slice Track Time : " << slice_time / N_t << std::endl;

    // // util::dump(Beam->dE, "dE\n", 10);
    // util::dump(Beam->dt, "dt\n", 10);
    // util::dump(Slice->n_macroparticles, "n_macroparticles\n", 10);

    // delete Context::Slice;
    // delete long_tracker;
    // delete Context::RfP;
    delete GP;
    delete Beam;
    delete slice;
    python::finalize();
    // delete psb;

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
            UNKNOWN, 0, "", "", Arg::None, "USAGE: ./TC8_Phase_Loop [options]\n\n"
            "Options:"
        },
        {
            HELP, 0, "h", "help", Arg::None,
            "  --help,              -h        Print usage and exit."
        },
        {
            N_TURNS, 0, "t", "turns", util::Arg::Numeric,
            "  --turns=<num>,       -t <num>  Number of turns (default: 500)"
        },
        {
            N_PARTICLES, 0, "p", "particles", util::Arg::Numeric,
            "  --particles=<num>,   -p <num>  Number of particles (default: "
            "100k)"
        },
        {
            N_SLICES, 0, "s", "slices", util::Arg::Numeric,
            "  --slices=<num>,      -s <num>  Number of slices (default: 200)"
        },
        {
            N_THREADS, 0, "m", "threads", util::Arg::Numeric,
            "  --threads=<num>,     -m <num>  Number of threads (default: 1)"
        },
        {
            UNKNOWN, 0, "", "", Arg::None,
            "\nExamples:\n"
            "\t./TC8_Phase_Loop\n"
            "\t./TC8_Phase_Loop -t 1000 -p 10000 -m 4\n"
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

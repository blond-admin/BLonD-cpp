/*
 * LHC_restart.cpp
 *
 *  Created on: Apr 12, 2016
 *      Author: kiliakis
 */
#include <blond/beams/Distributions.h>
#include <blond/globals.h>
#include <blond/llrf/PhaseLoop.h>
#include <blond/math_functions.h>
#include <blond/trackers/Tracker.h>
#include <blond/utilities.h>
#include <blond/llrf/LHCNoiseFB.h>
#include <blond/llrf/PhaseNoise.h>
#include <blond/impedances/InducedVoltage.h>
#include <blond/monitors/Monitors.h>

using namespace std;

// Simulation parameters
// --------------------------------------------------------
const long long N_b = 1.2e9;                // Intensity
int N_p = 600000;                           // Macro-particles

// Machine and RF parameters
const double C = 26658.883;                  // Machine circumference [m]
const int h = 35640;                        // Harmonic number
const double dphi = 0.;                      // Phase modulation/offset
const double gamma_t = 55.759505;            // Transition gamma
const double alpha = 1. / gamma_t / gamma_t; // First order mom. comp. factor

int alpha_order = 1;
int n_sections = 1;

// Tracking details
uint N_t = 9000001;       // Number of turns to track; full ramp: 8700001
ftype bl_target = 0.9e-9; // 4 sigma r.m.s. target bunch length in [s]

int N_slices = 2200;

int dt_save = 500000;
int dt_plot = 500000;

const std::string datafiles = "/afs/cern.ch/user/h/htimko/public/LHC/input/";




void parse_args(int argc, char **argv);

// Simulation setup
// -------------------------------------------------------------
int main(int argc, char **argv)
{

    parse_args(argc, argv);
    omp_set_num_threads(Context::n_threads);

    cout << "Setting up the simulation...\n\n";
    cout << "Number of turns: " << N_t << "\n";
    cout << "Number of macro-particles: " << N_p << "\n";
    cout << "Number of Slices: " << N_slices << "\n";
    cout << "Number of openmp threads: " << Context::n_threads << "\n";

    // printf("Setting up the simulation..\n");
    timespec begin, end;
    util::get_time(begin);

    f_vector_2d_t momentumVec(1, f_vector_t());
    util::read_vector_from_file(momentumVec[0],
                                datafiles + "LHC_momentum_programme.dat");
    momentumVec[0].resize(N_t + 1, momentumVec[0].back());

    f_vector_t V;
    util::read_vector_from_file(V, datafiles + "LHC_voltage_programme.dat");
    V.resize(N_t + 1, V.back());

    f_vector_2d_t voltageVec(1, V);
    printf("Flat top momentum %.4e eV\n", momentumVec[0][N_t]);
    printf("Flat top voltage %.4e eV\n", voltageVec[0][N_t]);
    cout << "Momentum and voltage loaded...\n";

    // Define general parameters

    f_vector_2d_t alphaVec(n_sections, f_vector_t(alpha_order + 1, alpha));

    f_vector_t CVec(n_sections, C);

    auto GP = Context::GP = new GeneralParameters(N_t, CVec, alphaVec, alpha_order,
            momentumVec, proton);

    cout << "General parameters set...\n";
    // Define rf_params
    f_vector_2d_t dphiVec(n_sections, f_vector_t(N_t + 1, dphi));

    f_vector_2d_t hVec(n_sections, f_vector_t(N_t + 1, h));

    auto RfP = Context::RfP = new RfParameters(n_sections, hVec, voltageVec, dphiVec);
    cout << "RF parameters set...\n";

    // Define beam and distribution: Load matched, filamented distribution
    auto Beam = Context::Beam = new Beams(N_p, N_b);
    f_vector_t v2;
    util::read_vector_from_file(v2, "/afs/cern.ch/work/k/kiliakis/testcases/"
                                "htimko/LHC/re_3b_extremes_batch/out/"
                                "initial_coords.dat");
    int k = 0;
    for (unsigned int i = 0; i < v2.size(); i += 3) {
        Context::Beam->dt[k] = v2[i];     // [s]
        Context::Beam->dE[k] = v2[i + 1]; // [eV]
        Context::Beam->id[k] = v2[i + 2];
        k++;
    }

    auto Slice = Context::Slice = new Slices(N_slices, 0, -1.0e-9, 54.0e-9);
    cout << "Beam generated, slices set...\n";

    // auto phaseNoise = new LHCFlatSpectrum(f_vector_t(10, 0), f_vector_t(10,
    // 0));

    auto phaseNoise = new LHCFlatSpectrum(1000000, 10000, 0.8571, 1.1, 1.0e-5, 1234, 7564,
                                          LHCFlatSpectrum::predistortion_t::weightfunction);

    phaseNoise->fDphi.clear();
    util::read_vector_from_file(
        phaseNoise->fDphi,
        datafiles + "LHCNoise_fmin0.8571_fmax1.001_ampl1e-5_weightfct.dat");

    auto noiseFB = new LHCNoiseFB(bl_target, 0.1e9, 0.93, 22500, true,
                                  f_vector_t{0, 10, 20});
    cout << "Phase noise feedback set...\n";

    // Define phase loop and frequency loop gain
    ftype PL_gain = 1 / (5 * Context::GP->t_rev[0]);
    ftype SL_gain = PL_gain / 10;

    f_vector_t PL_gainVec(N_t + 1, PL_gain);

    auto PL = new LHC(PL_gainVec, SL_gain, 0, phaseNoise, noiseFB);

    printf("\tPL gain is %.4e 1/s for initial turn T0 = %.4e s\n", PL->gain[0],
           Context::GP->t_rev[0]);
    printf("\tSL gain is %.4e turns\n", PL->gain2);
    printf("\tOmega_s0 = %.4e s at flat bottom, %.4e s at flat top\n",
           Context::RfP->omega_s0[0], Context::RfP->omega_s0[N_t]);
    printf("\tSL a_i = %.4f a_f = %.4f\n", PL->lhc_a[0], PL->lhc_a[N_t]);
    printf("\tSL t_i = %.4f t_f = %.4f\n", PL->lhc_t[0], PL->lhc_t[N_t]);

    // Injecting noise in the cavity, PL on
    auto long_tracker = new RingAndRfSection(Context::RfP, simple, PL);


    printf("PL, SL, and tracker set...\n");
    f_vector_t ZTot;
    util::read_vector_from_file(
        ZTot,
        datafiles + "Zlong_Allthemachine_450GeV_B1_LHC_inj_450GeV_B1.dat");
    assert(ZTot.size() % 3 == 0);

    f_vector_t freq, ReZ, ImZ;

    for (uint i = 0; i < ZTot.size(); i += 3) {
        freq.push_back(ZTot[i]);
        ReZ.push_back(ZTot[i + 1]);
        ImZ.push_back(ZTot[i + 2]);
    }

    auto ZTable = new InputTable(freq, ReZ, ImZ);
    vector<Intensity *> ZTableV{ZTable};

    auto indVoltage = new InducedVoltageFreq(ZTableV, 1.0e7);
    vector<InducedVoltage *> indVoltageV{indVoltage};

    auto totVoltage = new TotalInducedVoltage(indVoltageV);

    string h5file = "output_data-"
                    + to_string(Context::n_threads)
                    + "threads.h5";
    auto monitor = new BunchMonitor(GP, RfP, Beam, h5file,
                                    dt_save, Slice, PL, noiseFB);
    
    // monitor->fCompressionLevel = 0;
    // double slice_time = 0, track_time = 0;
    double turn_time = 0.0;

    timespec begin_t;
    printf("Map set\n");

    printf("Initial mean bunch position %.4e s\n", Context::Beam->mean_dt);
    printf("Initial four-times r.m.s. bunch length %.4e s\n",
           4. * Context::Beam->sigma_dt);
    // print("Initial Gaussian bunch length %.4e ns" %slices.bl_gauss

    printf("Ready for tracking!\n");

    auto totVoltageTime = 0.0;
    auto separatrixTime = 0.0;
    auto sliceTime = 0.0;
    auto trackerTime = 0.0;
    auto noiseTime = 0.0;
    auto monitorTime = 0.0;
    timespec start;

    for (uint i = 0; i < N_t; ++i) {

        util::get_time(begin_t);

        util::get_time(start);
        totVoltage->track();
        totVoltageTime += util::time_elapsed(start);

        // if (i % 10 == 0) {
        util::get_time(start);
        Beam->losses_separatrix(GP, RfP);
        separatrixTime += util::time_elapsed(start);
        // }

        util::get_time(start);
        Context::Slice->track();
        sliceTime += util::time_elapsed(start);

        util::get_time(start);
        long_tracker->track();
        trackerTime += util::time_elapsed(start);

        util::get_time(start);
        monitor->track();
        monitorTime += util::time_elapsed(start);

        util::get_time(start);
        noiseFB->track();
        noiseTime += util::time_elapsed(start);


        turn_time += util::time_elapsed(begin_t);

        if (i % dt_plot == 0) {
            printf("   Outputting at time step %d, tracking time %.4e s...\n", i, turn_time);
            printf("   Mean totVoltage track time %.4e\n", totVoltageTime / (i + 1));
            printf("   Mean Slices track time %.4e\n", sliceTime / (i + 1));
            printf("   Mean Tracker track time %.4e\n", trackerTime / (i + 1));
            printf("   Mean noiseFB track time %.4e\n", noiseTime / (i + 1));
            printf("   Mean Separatrix track time %.4e\n", separatrixTime / (i + 1));
            printf("   Mean Monitor track time %.4e\n", monitorTime / (i + 1));
            printf("   RF tracker counter is %d\n", RfP->counter);
            printf("   Beam momentum %0.6e eV\n", RfP->momentum(RfP->counter));
            printf("   Beam gamma %4.3f\n", RfP->gamma(RfP->counter));
            printf("   Beam beta %.8f\n", RfP->beta(RfP->counter));
            printf("   Beam energy %.6e eV\n", RfP->energy(RfP->counter));
            printf("   Design RF revolution frequency %.10e Hz\n", RfP->omega_RF_d[0][i]);
            printf("   RF revolution frequency %.10e Hz\n", RfP->omega_RF[0][i]);
            printf("   RF phase %.4f rad\n", RfP->phi_RF[0][i]);
            printf("   Beam phase %.4f rad\n", PL->phi_beam);
            printf("   Phase noise %.4f rad\n", noiseFB->fX * phaseNoise->fDphi[i]);
            printf("   PL phase error %.4f rad\n", PL->RFnoise->fDphi[i]);
            printf("   Synchronous phase %.4f rad\n", RfP->phi_s[i]);
            printf("   PL phase correction %.4f rad\n", PL->dphi);
            printf("   SL recursion variable %.4e\n", PL->lhc_y);
            printf("   Mean bunch position %.4e s\n", Beam->mean_dt);
            printf("   Four-times r.m.s. bunch length %.4e s\n", 4.0 * Beam->sigma_dt);
            // printf("   Gaussian bunch length %.4e s\n", Slice->bl_gauss);
            printf("   Slices min %.4e ns max %.4e s\n", Slice->cut_left, Slice->cut_right);
            printf("   Bunch 1 FWHM bunch length %.4e s\n", noiseFB->fBlMeasBBB[0]);
            printf("   Bunch 2 FWHM bunch length %.4e s\n", noiseFB->fBlMeasBBB[1]);
            printf("   Bunch 3 FWHM bunch length %.4e s\n", noiseFB->fBlMeasBBB[2]);
            printf("\n");
            printf("   Mean turn time:\t\t %.4lf s\n", turn_time / (i + 1));
            fflush(stdout);
        }

        // if (i % dt_save == 0) {
        //     ofstream out;
        //     out.open("out/coords_" + to_string(RfP->counter) + ".dat");
        //     cout.precision(10);
        //     cout << scientific << showpos;
        //     for (int j = 0; j < N_p; ++j) {
        //         out << Beam->dt[j] << "\t" << Beam->dE[j] << "\n";
        //     }
        //     out.close();
        // }

    }
    // ofstream out;
    // out.open("out/coords_" + to_string(RfP->counter) + ".dat");
    // cout.precision(10);
    // cout << scientific << showpos;
    // for (int j = 0; j < N_p; ++j) {
    //     out << Beam->dt[j] << "\t" << Beam->dE[j] << "\n";
    // }
    // out.close();

    util::get_time(end);
    util::print_time("Simulation Time in sec", begin, end);

    delete PL;
    delete Slice;
    delete long_tracker;
    delete RfP;
    delete GP;
    delete Beam;
    delete monitor;

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
            UNKNOWN, 0, "", "", Arg::None, "USAGE: ./_LHC_acc [options]\n\n"
            "Options:"
        },
        {
            HELP, 0, "h", "help", Arg::None,
            "  --help,              -h        Print usage and exit."
        },
        {
            N_TURNS, 0, "t", "turns", util::Arg::Numeric,
            "  --turns=<num>,       -t <num>  Number of turns (default: 1M)"
        },
        {
            N_PARTICLES, 0, "p", "particles", util::Arg::Numeric,
            "  --particles=<num>,   -p <num>  Number of particles (default: "
            "100k)"
        },
        {
            N_SLICES, 0, "s", "slices", util::Arg::Numeric,
            "  --slices=<num>,      -s <num>  Number of slices (default: 151)"
        },
        {
            N_THREADS, 0, "m", "threads", util::Arg::Numeric,
            "  --threads=<num>,     -m <num>  Number of threads (default: 1)"
        },
        {
            UNKNOWN, 0, "", "", Arg::None,
            "\nExamples:\n"
            "\t./_LHC_acc\n"
            "\t./_LHC_acc -t 1000000 -p 100000 -m 4\n"
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
        switch (opt.index()) {
            case HELP:
            case N_TURNS:
                N_t = atoi(opt.arg);
                break;
            case N_THREADS:
                Context::n_threads = atoi(opt.arg);
                break;
            case N_SLICES:
                N_slices = atoi(opt.arg);
                break;
            case N_PARTICLES:
                N_p = atoi(opt.arg);
                break;
            case UNKNOWN:
                break;
        }
    }
}

#include <blond/beams/Distributions.h>
#include <blond/globals.h>
#include <blond/llrf/PhaseLoop.h>
#include <blond/math_functions.h>
#include <blond/trackers/Tracker.h>
#include <blond/utilities.h>
#include <gtest/gtest.h>

// Simulation parameters
// --------------------------------------------------------

// Bunch parameters
const uint N_b = 0; // Intensity

// Machine and RF parameters
const double radius = 25;
const double C = 2 * constant::pi * radius;   // Machine circumference [m]
const double p_i = 310891054.809;             // Synchronous momentum [eV/c]
const uint h = 1;                            // Harmonic number
const double V = 8000;                        // RF voltage [V]
const double dphi = -constant::pi;            // Phase modulation/offset
const double gamma_t = 4.076750841;           // Transition gamma
const double alpha = 1.0 / gamma_t / gamma_t; // First order mom. comp. factor
const uint alpha_order = 1;
const uint n_sections = 1;
// Tracking details

uint N_t = 1000;   // Number of turns to track
uint N_p = 100000; // Macro-particles

uint N_slices = 200; // = (2^8)

// RingAndRFSection *long_tracker;

class testPLSPS_RL : public ::testing::Test {

  protected:
    virtual void SetUp() {
        omp_set_num_threads(1);
        
        f_vector_2d_t momentumVec(n_sections, f_vector_t(N_t + 1, p_i));

        f_vector_2d_t alphaVec(n_sections, f_vector_t(alpha_order, alpha));

        f_vector_t CVec(n_sections, C);

        f_vector_2d_t hVec(n_sections, f_vector_t(N_t + 1, h));

        f_vector_2d_t voltageVec(n_sections, f_vector_t(N_t + 1, V));

        f_vector_2d_t dphiVec(n_sections, f_vector_t(N_t + 1, dphi));

        Context::GP = new GeneralParameters(N_t, CVec, alphaVec,
                                             momentumVec,
                                            GeneralParameters::particle_t::proton);

        // Context::Beam = new Beams(N_p, N_b);

        auto GP = Context::GP;
        auto Beam = Context::Beam = new Beams(GP, N_p, N_b);
        
        auto RfP = Context::RfP = new RfParameters(GP, n_sections, hVec,
                                        voltageVec, dphiVec);


        // long_tracker = new RingAndRfSection();

        Context::Slice = new Slices(RfP, Beam, N_slices, 0, -constant::pi, constant::pi,
                                    Slices::cuts_unit_t::rad);
    }

    virtual void TearDown() {
        // Code here will be called immediately after each test
        // (right before the destructor).
        delete Context::GP;
        delete Context::Beam;
        delete Context::RfP;
        delete Context::Slice;
        // delete long_tracker;
    }
};

TEST_F(testPLSPS_RL, track1) {
    auto GP = Context::GP;
    auto Beam = Context::Beam;
    auto RfP = Context::RfP;
    longitudinal_bigaussian(GP, RfP, Beam, 200e-9, 100e6, 42321, false);

    auto sps = new SPS_RL(25e-6, 0, 5e-6);

    auto params = std::string(TEST_FILES "/") + "PL/SPS_RL/track1/";

    Context::Slice->track();
    f_vector_t domega_rf;

    for (uint i = 0; i < N_t; ++i) {
        sps->track();
        domega_rf.push_back(sps->domega_rf);
        Context::RfP->counter++;
    }

    f_vector_t v;

    util::read_vector_from_file(v, params + "domega_RF_mean.txt");
    // util::dump(domega_RF, "domega_RF");
    auto epsilon = 1e-2;
    auto ref = v[0];
    auto real = mymath::mean(domega_rf.data(), domega_rf.size());
    // std::cout << real << "\n";
    ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
        << "Testing of domega_RF_mean failed\n";

    v.clear();
    util::read_vector_from_file(v, params + "domega_RF_std.txt");
    epsilon = 1e-2;
    ref = v[0];
    real = mymath::standard_deviation(domega_rf.data(), domega_rf.size());
    // std::cout << real << "\n";
    ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
        << "Testing of domega_RF_std failed\n";

    delete sps;
}

TEST_F(testPLSPS_RL, track2) {
    auto GP = Context::GP;
    auto Beam = Context::Beam;
    auto RfP = Context::RfP;

    longitudinal_bigaussian(GP, RfP, Beam, 200e-9, 100e6, 42321, false);

    auto sps = new SPS_RL(25e-6, 0, 5e-6);
    sps->reference = 1e-6;

    auto params = std::string(TEST_FILES "/") + "PL/SPS_RL/track2/";

    Context::Slice->track();
    f_vector_t domega_rf;

    for (uint i = 0; i < N_t; ++i) {
        sps->track();
        domega_rf.push_back(sps->domega_rf);
        Context::RfP->counter++;
    }

    f_vector_t v;

    util::read_vector_from_file(v, params + "domega_RF_mean.txt");
    // util::dump(domega_RF, "domega_RF");
    auto epsilon = 1e-2;
    auto ref = v[0];
    auto real = mymath::mean(domega_rf.data(), domega_rf.size());
    // std::cout << real << "\n";
    ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
        << "Testing of domega_RF_mean failed\n";

    v.clear();
    util::read_vector_from_file(v, params + "domega_RF_std.txt");
    epsilon = 1e-2;
    ref = v[0];
    real = mymath::standard_deviation(domega_rf.data(), domega_rf.size());
    // std::cout << real << "\n";
    ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
        << "Testing of domega_RF_std failed\n";

    delete sps;
}

int main(int ac, char* av[]) {
    ::testing::InitGoogleTest(&ac, av);
    return RUN_ALL_TESTS();
}

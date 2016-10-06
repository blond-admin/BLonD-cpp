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

// RingAndRFSection *long_tracker;

class testPLPSB : public ::testing::Test {

  protected:
    virtual void SetUp() {
        omp_set_num_threads(1);
        
        f_vector_2d_t momentumVec(n_sections, f_vector_t(N_t + 1, p_i));

        f_vector_2d_t alphaVec(n_sections, f_vector_t(alpha_order + 1, alpha));

        f_vector_t CVec(n_sections, C);

        f_vector_2d_t hVec(n_sections, f_vector_t(N_t + 1, h));

        f_vector_2d_t voltageVec(n_sections, f_vector_t(N_t + 1, V));

        f_vector_2d_t dphiVec(n_sections, f_vector_t(N_t + 1, dphi));

        Context::GP = new GeneralParameters(N_t, CVec, alphaVec, alpha_order,
                                            momentumVec, proton);

        Context::Beam = new Beams(N_p, N_b);

        Context::RfP = new RfParameters(n_sections, hVec, voltageVec, dphiVec);

        // long_tracker = new RingAndRfSection();

        Context::Slice = new Slices(N_slices, 0, -constant::pi, constant::pi,
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

TEST_F(testPLPSB, constructor1) {

    // longitudinal_bigaussian(100e-9, 0.05e6, 1, false);

    auto psb =
        new PSB(f_vector_t(N_t, 1.0 / 25e-6), f_vector_t{0, 0}, 10e-6, 7);

    auto params = std::string(TEST_FILES "/") + "PL/PSB/constructor/test1/";
    f_vector_t v;

    util::read_vector_from_file(v, params + "gain2.txt");
    ASSERT_EQ(v.size(), psb->gain2.size());

    auto epsilon = 1e-8;
    for (uint i = 0; i < v.size(); ++i) {
        auto ref = v[i];
        auto real = psb->gain2[i];

        ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
            << "Testing of gain2 failed on i " << i << std::endl;
    }

    v.clear();
    util::read_vector_from_file(v, params + "gain.txt");
    ASSERT_EQ(v.size(), psb->gain.size());
    epsilon = 1e-8;
    for (uint i = 0; i < v.size(); ++i) {
        uint ref = v[i];
        uint real = psb->gain[i];

        ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
            << "Testing of gain failed on i " << i << std::endl;
    }

    v.clear();
    util::read_vector_from_file(v, params + "on_time.txt");
    ASSERT_EQ(v.size(), psb->on_time.size());

    for (uint i = 0; i < v.size(); ++i) {
        uint ref = v[i];
        uint real = psb->on_time[i];

        ASSERT_EQ(ref, real) << "Testing of on_time failed on i " << i
                             << std::endl;
    }

    delete psb;
}

TEST_F(testPLPSB, track1) {

    longitudinal_bigaussian(100e-9, 0.05e6, 1, false);

    auto psb =
        new PSB(f_vector_t(N_t, 1.0 / 25e-6), f_vector_t{0, 0}, 10e-6, 7);

    Context::Slice->track();
    f_vector_t dphi_av, t_accum, domega_PL, drho, domega_RL, domega_RF;

    for (uint i = 0; i < 500; ++i) {
        psb->track();
        dphi_av.push_back(psb->dphi_av);
        t_accum.push_back(psb->t_accum);
        domega_PL.push_back(psb->domega_PL);
        drho.push_back(psb->drho);
        domega_RL.push_back(psb->domega_RL);
        domega_RF.push_back(psb->domega_RF);
        Context::RfP->counter++;
    }

    auto params = std::string(TEST_FILES "/") + "PL/PSB/track/test1/";
    f_vector_t v;

    util::read_vector_from_file(v, params + "dphi_av_mean.txt");
    // util::dump(dphi_av, "dphi_av");
    auto epsilon = 1e-2;
    auto ref = v[0];
    auto real = mymath::mean(dphi_av.data(), dphi_av.size());
    ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
        << "Testing of dphi_av_mean failed\n";

    v.clear();
    util::read_vector_from_file(v, params + "dphi_av_std.txt");
    epsilon = 1e-2;
    ref = v[0];
    real = mymath::standard_deviation(dphi_av.data(), dphi_av.size());
    ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
        << "Testing of dphi_av_std failed\n";

    v.clear();
    util::read_vector_from_file(v, params + "t_accum_mean.txt");
    epsilon = 1e-2;
    ref = v[0];
    real = mymath::mean(t_accum.data(), t_accum.size());
    ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
        << "Testing of t_accum_mean failed\n";

    v.clear();
    util::read_vector_from_file(v, params + "t_accum_std.txt");
    epsilon = 1e-2;
    ref = v[0];
    real = mymath::standard_deviation(t_accum.data(), t_accum.size());
    ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
        << "Testing of t_accum_std failed\n";

    v.clear();
    util::read_vector_from_file(v, params + "domega_PL_mean.txt");
    epsilon = 1e-2;
    ref = v[0];
    real = mymath::mean(domega_PL.data(), domega_PL.size());
    ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
        << "Testing of domega_PL_mean failed\n";

    v.clear();
    util::read_vector_from_file(v, params + "domega_PL_std.txt");
    epsilon = 1e-2;
    ref = v[0];
    real = mymath::standard_deviation(domega_PL.data(), domega_PL.size());
    ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
        << "Testing of domega_PL_std failed\n";

    v.clear();
    util::read_vector_from_file(v, params + "drho_mean.txt");
    epsilon = 1e-2;
    ref = v[0];
    real = mymath::mean(drho.data(), drho.size());
    ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
        << "Testing of drho_mean failed\n";

    v.clear();
    util::read_vector_from_file(v, params + "drho_std.txt");
    epsilon = 1e-2;
    ref = v[0];
    real = mymath::standard_deviation(drho.data(), drho.size());
    ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
        << "Testing of drho_std failed\n";

    v.clear();
    util::read_vector_from_file(v, params + "domega_RL_mean.txt");
    epsilon = 1e-2;
    ref = v[0];
    real = mymath::mean(domega_RL.data(), domega_RL.size());
    ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
        << "Testing of domega_RL_mean failed\n";

    v.clear();
    util::read_vector_from_file(v, params + "domega_RL_std.txt");
    epsilon = 1e-2;
    ref = v[0];
    real = mymath::standard_deviation(domega_RL.data(), domega_RL.size());
    ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
        << "Testing of domega_RL_std failed\n";

    v.clear();
    util::read_vector_from_file(v, params + "domega_RF_mean.txt");
    epsilon = 1e-2;
    ref = v[0];
    real = mymath::mean(domega_RF.data(), domega_RF.size());
    ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
        << "Testing of domega_RF_mean failed\n";

    v.clear();
    util::read_vector_from_file(v, params + "domega_RF_std.txt");
    epsilon = 1e-2;
    ref = v[0];
    real = mymath::standard_deviation(domega_RF.data(), domega_RF.size());
    ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
        << "Testing of domega_RF_std failed\n";

    delete psb;
}

int main(int ac, char* av[]) {
    ::testing::InitGoogleTest(&ac, av);
    return RUN_ALL_TESTS();
}

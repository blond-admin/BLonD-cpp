#include <iostream>

#include <blond/beams/Distributions.h>
#include <blond/input_parameters/GeneralParameters.h>
#include <blond/llrf/PhaseLoop.h>
#include <blond/math_functions.h>
#include <blond/trackers/Tracker.h>
#include <blond/utilities.h>
#include <gtest/gtest.h>

const ftype epsilon = 1e-7;
const std::string params = TEST_FILES "/PL/PL_params/";

LHC* PL;
RingAndRfSection* long_tracker;

class testPL : public ::testing::Test {

  protected:
    // const ftype tau_0 = 0.4e-9;          // Initial bunch length, 4 sigma [s]

    virtual void SetUp() {
        f_vector_2d_t momentumVec(1, f_vector_t());
        util::read_vector_from_file(momentumVec[0],
                                    datafiles + "LHC_momentum_programme");

        // optional
        momentumVec[0].erase(momentumVec[0].begin(),
                             momentumVec[0].begin() + from_line);

        // std::cout << "vector size is " << v.size() << "\n";
        int remaining = N_t + 1 - momentumVec[0].size();
        for (int i = 0; i < remaining; ++i) {
            momentumVec[0].push_back(6.5e12);
        }
        assert(momentumVec[0].size() == N_t + 1);
        // ftype *V_array = new ftype[N_t + 1];
        f_vector_2d_t voltageVec(1, f_vector_t(N_t + 1));
        mymath::linspace(voltageVec[0].data(), 6e6, 10e6, 13563374, 13e6);
        std::fill_n(&voltageVec[0][563374], 436627, 10e6);

        // Define general parameters
        f_vector_2d_t alphaVec(alpha_order + 1, f_vector_t(n_sections, alpha));

        f_vector_t CVec(n_sections, C);

        Context::GP = new GeneralParameters(N_t, CVec, alphaVec, alpha_order,
                                            momentumVec, proton);
        auto GP = Context::GP;
        f_vector_2d_t dphiVec(n_sections, f_vector_t(N_t + 1, dphi));

        f_vector_2d_t hVec(n_sections, f_vector_t(N_t + 1, h));

        Context::RfP = new RfParameters(n_sections, hVec, voltageVec, dphiVec);

        // Define beam and distribution: Load matched, filamented distribution
        Context::Beam = new Beams(N_p, N_b);
        auto Beam = Context::Beam;
        f_vector_t v2;
        util::read_vector_from_file(v2, datafiles + "coords_13000001.dat");
        int k = 0;
        for (unsigned int i = 0; i < v2.size(); i += 3) {
            Beam->dt[k] = v2[i] * 1e-9;    // [s]
            Beam->dE[k] = v2[i + 1] * 1e6; // [eV]
            Beam->id[k] = v2[i + 2];
            k++;
        }

        Context::Slice = new Slices(N_slices, 0, -0.5e-9, 3e-9);
        // Define phase loop and frequency loop gain
        ftype PL_gain = 1 / (5 * GP->t_rev[0]);
        ftype SL_gain = PL_gain / 10;

        f_vector_t PL_gainVec(N_t + 1, PL_gain);

        PL = new LHC(PL_gainVec, SL_gain);

        // Injecting noise in the cavity, PL on
        long_tracker = new RingAndRfSection(Context::RfP, simple, PL);
    }

    virtual void TearDown() {
        // Code here will be called immediately after each test
        // (right before the destructor).
        delete Context::GP;
        delete Context::Beam;
        delete Context::RfP;
        delete Context::Slice;
        delete PL;
        delete long_tracker;
    }

  private:
    // Machine and RF parameters
    const float C = 26658.883;                  // Machine circumference [m]
    const int h = 35640;                        // Harmonic number
    const float dphi = 0.;                      // Phase modulation/offset
    const float gamma_t = 55.759505;            // Transition gamma
    const float alpha = 1. / gamma_t / gamma_t; // First order mom. comp. factor

    // Tracking details
    const int N_p = 100000;    // Macro-particles
    const long long N_b = 1e9; // Intensity
    const int alpha_order = 1;
    const int n_sections = 1;
    const uint N_t = 1000000; // Number of turns to track; full ramp: 8700001
    const ftype bl_target =
        1.25e-9; // 4 sigma r.m.s. target bunch length in [s]

    const int N_slices = 151;

    const std::string datafiles = DEMO_FILES "/LHC_restart/";

    const int from_line = 0;
};

class testPL2 : public ::testing::Test {

    virtual void SetUp() {
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
                                    cuts_unit_type::rad);
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

  protected:
    // Bunch parameters
    const uint N_b = 0; // Intensity

    // Machine and RF parameters
    const ftype radius = 25;
    const ftype C = 2 * constant::pi * radius; // Machine circumference [m]
    const ftype p_i = 310891054.809;           // Synchronous momentum [eV/c]
    const uint h = 1;                          // Harmonic number
    const ftype V = 8000;                      // RF voltage [V]
    const ftype dphi = -constant::pi;          // Phase modulation/offset
    const ftype gamma_t = 4.076750841;         // Transition gamma
    const ftype alpha =
        1.0 / gamma_t / gamma_t; // First order mom. comp. factor
    const uint alpha_order = 1;
    const uint n_sections = 1;
    // Tracking details

    uint N_t = 1000;   // Number of turns to track
    uint N_p = 500000; // Macro-particles

    uint N_slices = 200; // = (2^8)
};

TEST_F(testPL2, radial_difference1) {
    longitudinal_bigaussian(1e-9, 10e6, 42, false);
    auto long_tracker = RingAndRfSection();

    auto sps = new SPS_RL(1.0 / 25e-6, 0, 1e-6);

    auto params = std::string(TEST_FILES "/") + "PL/radial_difference/test1/";

    Context::Slice->track();
    f_vector_t drho;

    for (uint i = 0; i < N_t; ++i) {
        sps->radial_difference();
        long_tracker.track();
        drho.push_back(sps->drho);
    }

    f_vector_t v;

    util::read_vector_from_file(v, params + "drho_mean.txt");
    auto epsilon = 1e-2;
    auto ref = v[0];
    auto real = mymath::mean(drho.data(), drho.size());
    // ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
    //       << "Testing of drho_mean failed\n";

    v.clear();
    util::read_vector_from_file(v, params + "drho_std.txt");

    epsilon = 1e-2;
    ref = v[0];
    real = mymath::standard_deviation(drho.data(), drho.size());

    ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
        << "Testing of drho_std failed\n";

    delete sps;
}

TEST_F(testPL2, radial_steering_from_freq1) {
    auto RfP = Context::RfP;
    longitudinal_bigaussian(100e-9, 1e6, 1, false);
    auto long_tracker = RingAndRfSection();

    auto sps = new SPS_RL(1.0 / 25e-6, 0, 0);

    auto params =
        std::string(TEST_FILES "/") + "PL/radial_steering_from_freq/test1/";

    Context::Slice->track();

    N_t = 100;

    for (uint i = 0; i < N_t + 1; ++i)
        RfP->omega_RF[0][i] = std::sqrt(i);

    for (uint i = 0; i < N_t; ++i) {
        sps->radial_steering_from_freq();
        RfP->counter++;
    }

    f_vector_t v;
    util::read_vector_from_file(v, params + "omega_RF.txt");

    // ASSERT_EQ(v.size(), RfP->omega_RF[i].size());

    auto epsilon = 1e-8;
    for (unsigned int i = 0; i < v.size(); ++i) {
        auto ref = v[i];
        auto real = RfP->omega_RF[0][i];
        // std::cout << real << "\n";
        ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
            << "Testing of RfP->omega_RF failed on i " << i << std::endl;
    }

    v.clear();

    util::read_vector_from_file(v, params + "dphi_RF_steering.txt");

    ASSERT_EQ(v.size(), RfP->dphi_RF_steering.size());

    epsilon = 1e-8;
    for (unsigned int i = 0; i < v.size(); ++i) {
        auto ref = v[i];
        auto real = RfP->dphi_RF_steering[i];
        // std::cout << real << "\n";

        ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
            << "Testing of RfP->dphi_RF_steering failed on i " << i
            << std::endl;
    }

    v.clear();

    util::read_vector_from_file(v, params + "phi_RF.txt");

    // ASSERT_EQ(v.size(), RfP->phi_RF.size());

    epsilon = 1e-8;
    for (unsigned int i = 0; i < v.size(); ++i) {
        auto ref = v[i];
        auto real = RfP->phi_RF[0][i];
        // std::cout << real << "\n";

        ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
            << "Testing of RfP->phi_RF failed on i " << i << std::endl;
    }

    v.clear();

    delete sps;
}

TEST_F(testPL, lhc_a) {

    std::vector<ftype> v;
    util::read_vector_from_file(v, params + "lhc_a");
    // Only check 1 out of 10 elements
    // otherwise refernece file too big
    ASSERT_EQ(v.size(), Context::GP->n_turns / 10 + 1);
    for (unsigned int i = 0; i < v.size(); ++i) {
        // printf("ok here \n");

        ftype ref = v[i];
        ftype real = PL->lhc_a[i * 10];
        ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)));
    }

    // printf("ok here\n");
}

TEST_F(testPL, lhc_t) {

    std::vector<ftype> v;
    util::read_vector_from_file(v, params + "lhc_t");
    // Only check 1 out of 10 elements
    // otherwise refernece file too big
    ASSERT_EQ(v.size(), Context::GP->n_turns / 10 + 1);
    for (unsigned int i = 0; i < v.size(); ++i) {
        // printf("%d\n", i);
        ftype ref = v[i];
        ftype real = PL->lhc_t[i * 10];
        ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)));
    }
}

TEST_F(testPL, phi_beam) {
    Context::Slice->track();
    PL->beam_phase();

    std::vector<ftype> v;
    util::read_vector_from_file(v, params + "phi_beam");
    ftype ref = v[0];
    ftype real = PL->phi_beam;
    ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)));
}

TEST_F(testPL, dphi) {
    Context::Slice->track();
    PL->beam_phase();
    PL->phase_difference();

    std::vector<ftype> v;
    util::read_vector_from_file(v, params + "dphi");
    ftype ref = v[0];
    ftype real = PL->dphi;
    ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)));
}

TEST_F(testPL, domega_RF) {
    Context::Slice->track();
    PL->beam_phase();
    PL->phase_difference();
    long_tracker->track();
    std::vector<ftype> v;
    util::read_vector_from_file(v, params + "domega_RF");
    ftype ref = v[0];
    ftype real = PL->domega_RF;
    ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)));
}

TEST_F(testPL, lhc_y) {
    Context::Slice->track();
    PL->beam_phase();
    PL->phase_difference();
    long_tracker->track();
    std::vector<ftype> v;
    util::read_vector_from_file(v, params + "lhc_y");
    ftype ref = v[0];
    ftype real = PL->lhc_y;
    ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)));
}

TEST_F(testPL, omega_RF) {
    auto RfP = Context::RfP;
    Context::Slice->track();

    PL->beam_phase();
    PL->phase_difference();
    long_tracker->track();
    // RfP->counter++;
    int counter = RfP->counter;

    std::vector<ftype> v;
    util::read_vector_from_file(v, params + "omega_RF");
    ftype ref = v[0];
    ftype real = RfP->omega_RF[0][counter];
    ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)));
}

TEST_F(testPL, dphi_RF) {
    Context::Slice->track();
    PL->beam_phase();
    PL->phase_difference();
    long_tracker->track();
    // RfP->counter++;

    std::vector<ftype> v;
    util::read_vector_from_file(v, params + "dphi_RF");
    ftype ref = v[0];
    ftype real = Context::RfP->dphi_RF[0];
    ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)));
}

TEST_F(testPL, phi_RF) {
    auto RfP = Context::RfP;

    Context::Slice->track();
    PL->beam_phase();
    PL->phase_difference();
    long_tracker->track();
    // RfP->counter++;

    int counter = RfP->counter;

    std::vector<ftype> v;
    util::read_vector_from_file(v, params + "phi_RF");
    ftype ref = v[0];
    ftype real = RfP->phi_RF[0][counter];
    ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)));
}

int main(int ac, char* av[]) {
    ::testing::InitGoogleTest(&ac, av);
    return RUN_ALL_TESTS();
}

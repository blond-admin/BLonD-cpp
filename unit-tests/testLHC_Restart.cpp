
#include <iostream>
#include <list>
#include <string>

#include <blond/beams/Distributions.h>
#include <blond/input_parameters/GeneralParameters.h>
#include <blond/llrf/PhaseLoop.h>
#include <blond/math_functions.h>
#include <blond/trackers/Tracker.h>
#include <blond/utilities.h>
#include <gtest/gtest.h>
#include <omp.h>

// const ftype epsilon = 1e-3;
const std::string params = TEST_FILES "/PL/LHC_restart_params/";

LHC* PL;
RingAndRfSection* long_tracker;

class testLHC_Restart : public ::testing::Test {

  protected:
    const int N_p = 100000; // Macro-particles
    // const ftype tau_0 = 0.4e-9;          // Initial bunch length, 4 sigma [s]
    const uint N_t = 1000000; // Number of turns to track; full ramp: 8700001

    const int N_slices = 151;

    virtual void SetUp() {
        // printf("ok here\n");
        omp_set_num_threads(1);

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

        f_vector_2d_t voltageVec(1, f_vector_t(N_t + 1));
        mymath::linspace(voltageVec[0].data(), 6e6, 10e6, 13563374, 13e6);
        std::fill_n(&voltageVec[0][563374], 436627, 10e6);

        // Define general parameters
        f_vector_2d_t alphaVec(n_sections, f_vector_t(alpha_order + 1, alpha));

        f_vector_t CVec(n_sections, C);

        Context::GP = new GeneralParameters(N_t, CVec, alphaVec, alpha_order,
                                            momentumVec, proton);
        auto GP = Context::GP;
        // Define rf_params
        f_vector_2d_t dphiVec(n_sections, f_vector_t(N_t + 1, dphi));

        f_vector_2d_t hVec(n_sections, f_vector_t(N_t + 1, h));

        Context::RfP = new RfParameters(n_sections, hVec, voltageVec, dphiVec);

        // Define beam and distribution: Load matched, filamented distribution
        Context::Beam = new Beams(N_p, N_b);
        auto Beam = Context::Beam;
        std::vector<ftype> v2;
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

    const long long N_b = 1e9; // Intensity
    const int alpha_order = 1;
    const int n_sections = 1;
    const int bl_target = 1.25e-9; // 4 sigma r.m.s. target bunch length in [s]

    const std::string datafiles = DEMO_FILES "/LHC_restart/";

    const int from_line = 0;
};

TEST_F(testLHC_Restart, dphi_RF_and_dphi) {
    auto RfP = Context::RfP;
    auto Beam = Context::Beam;
    auto Slice = Context::Slice;

    f_vector_t real1, real2;

    for (int i = 0; i < 1000; ++i) {

        if (RfP->counter < 570000)
            PL->reference = 0.5236;
        else
            PL->reference = 1.0472;

        Context::Slice->track();

        long_tracker->track();

        // printf("   Beam energy %.6e eV\n", GP->energy[0]);
        // printf("   RF phase %.6e rad\n", RfP->dphi_RF[0]);
        // printf("   PL phase correction %.6e rad\n", PL->dphi);
        real1.push_back(RfP->dphi_RF[0]);
        real2.push_back(PL->dphi);
    }

    f_vector_t v;

    util::read_vector_from_file(v, params + "dphi_RF[0]");

    ASSERT_EQ(v.size(), real1.size());
    ftype epsilon = 1e-6;
    for (unsigned int i = 0; i < v.size(); ++i) {
        // printf("%d\n", i);
        ftype ref = v[i];
        ftype real = real1[i];
        ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
            << "Testing of dphi_RF failed on i " << i << std::endl;
    }

    epsilon = 1e-3;
    v.clear();
    util::read_vector_from_file(v, params + "dphi");
    ASSERT_EQ(v.size(), real2.size());
    for (unsigned int i = 0; i < v.size(); ++i) {
        // printf("%d\n", i);
        ftype ref = v[i];
        ftype real = real2[i];
        ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
            << "Testing of dphi failed on i " << i << std::endl;
    }

    epsilon = 5e-1;
    v.clear();
    util::read_vector_from_file(v, params + "dE");
    ASSERT_EQ(v.size(), Beam->n_macroparticles);
    for (unsigned int i = 0; i < v.size(); ++i) {
        // printf("%d\n", i);
        ftype ref = v[i];
        ftype real = Beam->dE[i];
        ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
            << "Testing of dE failed on i " << i << std::endl;
    }

    epsilon = 1e-6;
    v.clear();
    util::read_vector_from_file(v, params + "dt");
    ASSERT_EQ(v.size(), Beam->n_macroparticles);
    for (unsigned int i = 0; i < v.size(); ++i) {
        // printf("%d\n", i);
        ftype ref = v[i];
        ftype real = Beam->dt[i];
        ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
            << "Testing of dt failed on i " << i << std::endl;
    }

    epsilon = 1e-8;
    v.clear();
    util::read_vector_from_file(v, params + "n_macroparticles");
    ASSERT_EQ(v.size(), Slice->n_slices);
    for (unsigned int i = 0; i < v.size(); ++i) {
        // printf("%d\n", i);
        ftype ref = v[i];
        ftype real = Slice->n_macroparticles[i];
        ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
            << "Testing of n_macroparticles failed on i " << i << std::endl;
    }
}

int main(int ac, char* av[]) {
    ::testing::InitGoogleTest(&ac, av);
    return RUN_ALL_TESTS();
}

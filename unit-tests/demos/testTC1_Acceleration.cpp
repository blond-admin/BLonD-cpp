#include <iostream>

#include <blond/beams/Distributions.h>
#include <blond/input_parameters/GeneralParameters.h>
#include <blond/math_functions.h>
#include <blond/trackers/Tracker.h>
#include <blond/utilities.h>
#include <gtest/gtest.h>

using namespace std;

class testTC1 : public ::testing::Test {

protected:
    // Machine and RF parameters
    const double C = 26658.883;       // Machine circumference [m]
    const double p_i = 450e9;     // Synchronous momentum [eV/c]
    const double p_f = 460.005e9;     // Synchronous momentum, final
    const double h = 35640;       // Harmonic number
    const double V = 6e6;             // RF voltage [V]
    const double dphi = 0;            // Phase modulation/offset
    const double gamma_t = 55.759505; // Transition gamma
    const double alpha =
        1.0 / gamma_t / gamma_t; // First order mom. comp. factor
    const int alpha_order = 1;
    const int n_sections = 1;
    // Tracking details
    const long long N_b = 1e9;  // Intensity
    const double tau_0 = 0.4e-9; // Initial bunch length, 4 sigma [s]

    const int N_t = 2000; // Number of turns to track
    const int N_p = 10000;  // Macro-particles
    const int N_slices = 100;


    virtual void SetUp()
    {
        omp_set_num_threads(4);


        f_vector_2d_t momentumVec(n_sections, f_vector_t(N_t + 1));
        for (auto &v : momentumVec)
            mymath::linspace(v.data(), p_i, p_f, N_t + 1);

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

        Context::RfP = new RfParameters(GP, n_sections, hVec,
                                        voltageVec, dphiVec);



        // Context::Slice = new Slices(N_slices);
    }

    virtual void TearDown()
    {
        // Code here will be called immediately after each test
        // (right before the destructor).
        delete Context::GP;
        delete Context::Beam;
        delete Context::RfP;
        delete Context::Slice;
    }

};

TEST_F(testTC1, phaseSpace1)
{
    auto Beam = Context::Beam;
    auto RfP = Context::RfP;
    auto GP = Context::GP;
    longitudinal_bigaussian(GP, RfP, Beam, tau_0 / 4, 0, -1, false);
    auto epsilon = 1e-6;
    omp_set_num_threads(Context::n_threads);
    string params = TEST_FILES "/TC1_final/phaseSpace1/";

    auto slices = new Slices(RfP, Beam, N_slices);

    auto long_tracker = new RingAndRfSection();
    for (int i = 0; i < N_t; ++i) {
        long_tracker->track();
        slices->track();
    }

    f_vector_t v;
    util::read_vector_from_file(v, params + "dE.txt");
    ASSERT_EQ(v.size(), Beam->dE.size());
    for (uint i = 0; i < v.size(); ++i) {
        auto ref = v[i];
        auto real = Beam->dE[i];
        ASSERT_NEAR(ref, real, epsilon * max(abs(ref), abs(real)))
                << "Testing of dE failed on i " << i << '\n';

    }


    util::read_vector_from_file(v, params + "dt.txt");
    ASSERT_EQ(v.size(), Beam->dt.size());
    for (uint i = 0; i < v.size(); ++i) {
        auto ref = v[i];
        auto real = Beam->dt[i];
        ASSERT_NEAR(ref, real, epsilon * max(abs(ref), abs(real)))
                << "Testing of dt failed on i " << i << '\n';

    }

    util::read_vector_from_file(v, params + "n_macroparticles.txt");
    ASSERT_EQ(v.size(), slices->n_macroparticles.size());
    for (uint i = 0; i < v.size(); ++i) {
        auto ref = v[i];
        auto real = 1.0 * slices->n_macroparticles[i];
        ASSERT_EQ(ref, real)
                << "Testing of n_macroparticles failed on i " << i << '\n';

    }

    delete long_tracker;
    delete slices;
}


TEST_F(testTC1, phaseSpace2)
{
    auto Beam = Context::Beam;
    auto GP = Context::GP;

    auto RfP = Context::RfP;
    longitudinal_bigaussian(GP, RfP, Beam, tau_0 / 4, 0, -1, false);
    auto epsilon = 1e-6;

    Context::n_threads = 3;
    auto slices = new Slices(RfP, Beam, N_slices);
    omp_set_num_threads(Context::n_threads);
    string params = TEST_FILES "/TC1_final/phaseSpace1/";

    auto long_tracker = new RingAndRfSection();
    for (int i = 0; i < N_t; ++i) {
        long_tracker->track();
        slices->track();
    }

    f_vector_t v;
    util::read_vector_from_file(v, params + "dE.txt");
    ASSERT_EQ(v.size(), Beam->dE.size());
    for (uint i = 0; i < v.size(); ++i) {
        auto ref = v[i];
        auto real = Beam->dE[i];
        ASSERT_NEAR(ref, real, epsilon * max(abs(ref), abs(real)))
                << "Testing of dE failed on i " << i << '\n';

    }


    util::read_vector_from_file(v, params + "dt.txt");
    ASSERT_EQ(v.size(), Beam->dt.size());
    for (uint i = 0; i < v.size(); ++i) {
        auto ref = v[i];
        auto real = Beam->dt[i];
        ASSERT_NEAR(ref, real, epsilon * max(abs(ref), abs(real)))
                << "Testing of dt failed on i " << i << '\n';

    }

    util::read_vector_from_file(v, params + "n_macroparticles.txt");
    ASSERT_EQ(v.size(), slices->n_macroparticles.size());
    for (uint i = 0; i < v.size(); ++i) {
        int ref = v[i];
        int real = slices->n_macroparticles[i];
        ASSERT_EQ(ref, real)
                << "Testing of n_macroparticles failed on i " << i << '\n';

    }

    delete long_tracker;
    delete slices;
}

int main(int ac, char *av[])
{
    ::testing::InitGoogleTest(&ac, av);
    return RUN_ALL_TESTS();
}

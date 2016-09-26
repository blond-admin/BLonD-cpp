#include <iostream>
#include <list>
#include <string>

#include <blond/beams/Distributions.h>
#include <blond/constants.h>
#include <blond/input_parameters/GeneralParameters.h>
#include <blond/math_functions.h>
#include <blond/trackers/Tracker.h>
#include <blond/utilities.h>
#include <gtest/gtest.h>

using namespace std;

class testBeam : public ::testing::Test {

protected:
    const long long N_b = 1e9;      // Intensity
    const ftype tau_0 = 0.4e-9;     // Initial bunch length, 4 sigma [s]
    // Machine and RF parameters
    const ftype C = 26658.883;       // Machine circumference [m]
    const long long p_i = 450e9;     // Synchronous momentum [eV/c]
    const ftype p_f = 460.005e9;     // Synchronous momentum, final
    const long long h = 35640;       // Harmonic number
    const ftype V = 6e6;             // RF voltage [V]
    const ftype dphi = 0;            // Phase modulation/offset
    const ftype gamma_t = 55.759505; // Transition gamma
    const ftype alpha =
        1.0 / gamma_t / gamma_t;     // First order mom. comp. factor
    const int alpha_order = 1;
    const int n_sections = 1;
    // Tracking details

    const int N_t = 2000;            // Number of turns to track
    const int N_p = 1000;            // Macro-particles
    const int N_slices = 10;
    virtual void SetUp()
    {
        omp_set_num_threads(1);
        
        f_vector_2d_t momentumVec(n_sections, f_vector_t(N_t + 1, 0));
        for (auto &v : momentumVec)
            mymath::linspace(v.data(), p_i, p_f, N_t + 1);

        f_vector_2d_t alphaVec(n_sections, f_vector_t(alpha_order + 1, alpha));

        f_vector_t CVec(n_sections, C);

        f_vector_2d_t hVec(n_sections, f_vector_t(N_t + 1, h));

        f_vector_2d_t voltageVec(n_sections, f_vector_t(N_t + 1, V));

        f_vector_2d_t dphiVec(n_sections, f_vector_t(N_t + 1, dphi));

        Context::GP = new GeneralParameters(N_t, CVec, alphaVec, alpha_order,
                                            momentumVec, proton);

        Context::Beam = new Beams(N_p, N_b);

        Context::RfP = new RfParameters(n_sections, hVec, voltageVec, dphiVec);


    }

    virtual void TearDown()
    {
        // Code here will be called immediately after each test
        // (right before the destructor).
        delete Context::GP;
        delete Context::Beam;
        delete Context::RfP;
    }
};


TEST_F(testBeam, statistics1)
{
    std::string params =
        TEST_FILES "/Beam/statistics1/";
    longitudinal_bigaussian(tau_0 / 4, 0, -1, false);
    Context::Beam->statistics();

    auto epsilon = 1e-8;
    f_vector_t v;

    util::read_vector_from_file(v, params + "sigma_dE.txt");
    auto ref = v[0];
    auto real = Context::Beam->sigma_dE;
    ASSERT_NEAR(ref, real, epsilon * std::max(std::abs(ref), std::abs(real)))
            << "Testing of sigma_dE failed\n";

    util::read_vector_from_file(v, params + "sigma_dt.txt");
    ref = v[0];
    real = Context::Beam->sigma_dt;
    ASSERT_NEAR(ref, real, epsilon * std::max(std::abs(ref), std::abs(real)))
            << "Testing of sigma_dt failed\n";

    util::read_vector_from_file(v, params + "mean_dE.txt");
    ref = v[0];
    real = Context::Beam->mean_dE;
    ASSERT_NEAR(ref, real, epsilon * std::max(std::abs(ref), std::abs(real)))
            << "Testing of mean_dE failed\n";

    util::read_vector_from_file(v, params + "mean_dt.txt");
    ref = v[0];
    real = Context::Beam->mean_dt;
    ASSERT_NEAR(ref, real, epsilon * std::max(std::abs(ref), std::abs(real)))
            << "Testing of mean_dt failed\n";

    util::read_vector_from_file(v, params + "epsn_rms_l.txt");
    ref = v[0];
    real = Context::Beam->epsn_rms_l;
    ASSERT_NEAR(ref, real, epsilon * std::max(std::abs(ref), std::abs(real)))
            << "Testing of epsn_rms_l failed\n";

    util::read_vector_from_file(v, params + "n_macroparticles_lost.txt");
    ref = v[0];
    real = Context::Beam->n_macroparticles_lost;
    ASSERT_NEAR(ref, real, epsilon * std::max(std::abs(ref), std::abs(real)))
            << "Testing of n_macroparticles_lost failed\n";
}

TEST_F(testBeam, losses_long_cut1)
{
    auto epsilon = 1e-8;

    std::string params =
        TEST_FILES "/Beam/losses_long_cut1/";
    auto Beam = Context::Beam;
    longitudinal_bigaussian(tau_0 / 4, 0, -1, false);
    Beam->statistics();
    Beam->losses_longitudinal_cut(Beam->mean_dt, 10 * std::fabs(Beam->mean_dt));

    f_vector_t v;
    util::read_vector_from_file(v, params + "id.txt");
    ASSERT_EQ(v.size(), Beam->id.size());

    for (uint i = 0; i < v.size(); ++i) {
        auto ref = v[i];
        auto real = Beam->id[i];
        ASSERT_NEAR(ref, real, epsilon * std::max(std::fabs(ref), std::fabs(real)))
                << "Testing of id failed on i " << i << '\n';
    }
}

TEST_F(testBeam, losses_energy_cut1)
{
    auto epsilon = 1e-8;

    std::string params =
        TEST_FILES "/Beam/losses_energy_cut1/";
    auto Beam = Context::Beam;
    longitudinal_bigaussian(tau_0 / 4, 0, -1, false);
    Beam->statistics();

    Beam->losses_energy_cut(Beam->mean_dE, 10 * std::fabs(Beam->mean_dE));

    f_vector_t v;
    util::read_vector_from_file(v, params + "id.txt");
    ASSERT_EQ(v.size(), Beam->id.size());

    for (uint i = 0; i < v.size(); ++i) {
        auto ref = v[i];
        auto real = Beam->id[i];
        ASSERT_NEAR(ref, real, epsilon * std::max(std::fabs(ref), std::fabs(real)))
                << "Testing of id failed on i " << i << '\n';
    }
}



class testBeam2 : public ::testing::Test {

protected:
    const long long N_b = 1e9;      // Intensity
    const ftype tau_0 = 0.4e-9;     // Initial bunch length, 4 sigma [s]
    // Machine and RF parameters
    const ftype C = 26658.883;       // Machine circumference [m]
    const ftype p_i = 450e9;     // Synchronous momentum [eV/c]
    const ftype p_f = 478e9;     // Synchronous momentum, final
    const long long h = 35640;       // Harmonic number
    const ftype V = 6e6;             // RF voltage [V]
    const ftype dphi = 0;            // Phase modulation/offset
    const ftype gamma_t = 55.759505; // Transition gamma
    const ftype alpha =
        1.0 / gamma_t / gamma_t;     // First order mom. comp. factor
    const int alpha_order = 1;
    const int n_sections = 1;
    // Tracking details

    const int N_t = 5000;            // Number of turns to track
    const int N_p = 10000;            // Macro-particles
    const int N_slices = 200;
    virtual void SetUp()
    {
        f_vector_2d_t momentumVec(n_sections, f_vector_t(N_t + 1, 0));
        for (auto &v : momentumVec)
            mymath::linspace(v.data(), p_i, p_f, N_t + 1);

        f_vector_2d_t alphaVec(n_sections, f_vector_t(alpha_order + 1, alpha));

        f_vector_t CVec(n_sections, C);

        f_vector_2d_t hVec(n_sections, f_vector_t(N_t + 1, h));

        f_vector_2d_t voltageVec(n_sections, f_vector_t(N_t + 1, V));

        f_vector_2d_t dphiVec(n_sections, f_vector_t(N_t + 1, dphi));

        Context::GP = new GeneralParameters(N_t, CVec, alphaVec, alpha_order,
                                            momentumVec, proton);

        Context::Beam = new Beams(N_p, N_b);

        Context::RfP = new RfParameters(n_sections, hVec, voltageVec, dphiVec);


    }

    virtual void TearDown()
    {
        // Code here will be called immediately after each test
        // (right before the destructor).
        delete Context::GP;
        delete Context::Beam;
        delete Context::RfP;
    }
};


TEST_F(testBeam2, losses_separatrix1)
{
    auto params = string(TEST_FILES) +  "/Beam/losses_separatrix1/";
    auto GP = Context::GP;
    auto RfP = Context::RfP;
    auto Beam = Context::Beam;

    longitudinal_bigaussian(tau_0 / 4, 0, -1, false);
    Beam->losses_separatrix(GP, RfP);

    f_vector_t v;
    util::read_vector_from_file(v, params + "id.txt");
    ASSERT_EQ(v.size(), Beam->id.size());

    for (uint i = 0; i < v.size(); ++i) {
        int ref = v[i];
        int real = Beam->id[i];
        ASSERT_EQ(ref, real)
                << "Testing of id failed on i " << i << '\n';
    }
}


TEST_F(testBeam2, losses_separatrix2)
{
    auto params = string(TEST_FILES) +  "/Beam/losses_separatrix2/";
    auto GP = Context::GP;
    auto RfP = Context::RfP;
    auto Beam = Context::Beam;

    longitudinal_bigaussian(tau_0 / 3, 0, -1, false);
    auto tracker = RingAndRfSection();

    for (int i = 0; i < 500; i++) tracker.track();
    Beam->losses_separatrix(GP, RfP);

    f_vector_t v;
    util::read_vector_from_file(v, params + "id.txt");
    ASSERT_EQ(v.size(), Beam->id.size());

    for (uint i = 0; i < v.size(); ++i) {
        int ref = v[i];
        int real = Beam->id[i];
        ASSERT_EQ(ref, real)
                << "Testing of id failed on i " << i << '\n';
    }
}


int main(int ac, char *av[])
{
    ::testing::InitGoogleTest(&ac, av);
    return RUN_ALL_TESTS();
}

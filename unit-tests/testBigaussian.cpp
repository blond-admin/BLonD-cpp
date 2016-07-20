#include <iostream>
#include <list>
#include <string>

#include <blond/beams/Distributions.h>
#include <blond/constants.h>
#include <blond/input_parameters/GeneralParameters.h>
#include <blond/math_functions.h>
#include <blond/utilities.h>
#include <gtest/gtest.h>

const ftype epsilon = 1e-8;
const std::string fixed_params =
    "../unit-tests/references/Bigaussian/Bigaussian_fixed_params/";

class testBigaussian : public ::testing::Test {

  protected:
    const long long N_b = 1e9;       // Intensity
    const ftype tau_0 = 0.4e-9; // Initial bunch length, 4 sigma [s]

    virtual void SetUp() {
        f_vector_2d_t momentumVec(n_sections, f_vector_t(N_t + 1));
        for (auto& v : momentumVec)
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

    virtual void TearDown() {
        // Code here will be called immediately after each test
        // (right before the destructor).
        delete Context::GP;
        delete Context::Beam;
        delete Context::RfP;
    }

  private:
    // Machine and RF parameters
    const ftype C = 26658.883;       // Machine circumference [m]
    const long long p_i = 450e9;          // Synchronous momentum [eV/c]
    const ftype p_f = 460.005e9;     // Synchronous momentum, final
    const long long h = 35640;            // Harmonic number
    const ftype V = 6e6;             // RF voltage [V]
    const ftype dphi = 0;            // Phase modulation/offset
    const ftype gamma_t = 55.759505; // Transition gamma
    const ftype alpha =
        1.0 / gamma_t / gamma_t; // First order mom. comp. factor
    const int alpha_order = 1;
    const int n_sections = 1;
    // Tracking details

    const int N_t = 2000; // Number of turns to track
    const int N_p = 100;  // Macro-particles
    const int N_slices = 10;
};

class testBigaussianRandom : public ::testing::Test {

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
    }

    virtual void TearDown() {
        // Code here will be called immediately after each test
        // (right before the destructor).
        delete Context::GP;
        delete Context::Beam;
        delete Context::RfP;
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
    const ftype tau_0 = 0.4e-9; // Initial bunch length, 4 sigma [s]

    // Tracking details

    const int N_t = 100;    // Number of turns to track
    const int N_p = 500000; // Macro-particles
    const int N_slices = 10;
};

TEST_F(testBigaussian, test_sigma_dE) {

    longitudinal_bigaussian(tau_0 / 4, 0, -1, false);

    std::vector<ftype> v;

    util::read_vector_from_file(v, fixed_params + "sigma_dE");
    ftype ref = v[0];
    ftype real = Context::Beam->sigma_dE;
    ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)));
}

TEST_F(testBigaussian, test_sigma_dt) {

    longitudinal_bigaussian(tau_0 / 4, 0, -1, false);

    std::vector<ftype> v;

    util::read_vector_from_file(v, fixed_params + "sigma_dt");
    ftype ref = v[0];
    ftype real = Context::Beam->sigma_dt;
    ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)));
}

TEST_F(testBigaussian, test_dE) {

    longitudinal_bigaussian(tau_0 / 4, 0, -1, false);

    std::vector<ftype> v;
    util::read_vector_from_file(v, fixed_params + "dE");
    for (unsigned int i = 0; i < v.size(); ++i) {
        ftype ref = v[i];
        ftype real = Context::Beam->dE[i];
        ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)));
    }
}

TEST_F(testBigaussian, test_dt) {

    longitudinal_bigaussian(tau_0 / 4, 0, -1, false);

    std::vector<ftype> v;
    util::read_vector_from_file(v, fixed_params + "dt");
    for (unsigned int i = 0; i < v.size(); ++i) {
        ftype ref = v[i];
        ftype real = Context::Beam->dt[i];
        ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)));
    }
}

TEST_F(testBigaussianRandom, test_dE) {
    auto Beam = Context::Beam;

    auto params =
        std::string("../unit-tests/references/") + "Bigaussian/random/";

    longitudinal_bigaussian(100e-9, 0.05e6, 1, false);

    f_vector_t v;

    util::read_vector_from_file(v, params + "dE_mean.txt");
    auto epsilon = 1e-3;
    auto ref = v[0];
    auto real = mymath::mean(Beam->dE.data(), Beam->dE.size());
    // ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
    //       << "Testing of Beam->dE_mean failed\n";

    v.clear();
    util::read_vector_from_file(v, params + "dE_std.txt");
    epsilon = 1e-2;
    ref = v[0];
    real = mymath::standard_deviation(Beam->dE.data(), Beam->dE.size());
    ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
        << "Testing of Beam->dE_std failed\n";
}

TEST_F(testBigaussianRandom, test_dt) {
    auto Beam = Context::Beam;

    auto params =
        std::string("../unit-tests/references/") + "Bigaussian/random/";

    longitudinal_bigaussian(100e-9, 0.05e6, 1, false);

    f_vector_t v;

    util::read_vector_from_file(v, params + "dt_mean.txt");
    auto epsilon = 1e-2;
    auto ref = v[0];
    auto real = mymath::mean(Beam->dt.data(), Beam->dt.size());
    ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
        << "Testing of Beam->dt_mean failed\n";

    v.clear();
    util::read_vector_from_file(v, params + "dt_std.txt");
    epsilon = 1e-2;
    ref = v[0];
    real = mymath::standard_deviation(Beam->dt.data(), Beam->dt.size());
    ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
        << "Testing of Beam->dt_std failed\n";
}

int main(int ac, char* av[]) {
    ::testing::InitGoogleTest(&ac, av);
    return RUN_ALL_TESTS();
}

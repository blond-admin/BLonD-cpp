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
const std::string statistics_params =
    TEST_FILES "/Beam/Beam_statistics_params/";
const std::string long_cut_params = TEST_FILES "/Beam/Beam_long_cut_params/";
const std::string energy_cut_params =
    TEST_FILES "/Beam/Beam_energy_cut_params/";

class testBeam : public ::testing::Test {

  protected:
    const long long N_b = 1e9;  // Intensity
    const ftype tau_0 = 0.4e-9; // Initial bunch length, 4 sigma [s]

    virtual void SetUp() {
        f_vector_2d_t momentumVec(n_sections, f_vector_t(N_t + 1, 0));
        for (auto& v : momentumVec)
            mymath::linspace(v.data(), p_i, p_f, N_t + 1);

        // f_vector_2d_t alphaVec(alpha_order + 1, f_vector_t(n_sections,
        // alpha));
        f_vector_2d_t alphaVec(n_sections, f_vector_t(alpha_order + 1, alpha));

        f_vector_t CVec(n_sections, C);

        f_vector_2d_t hVec(n_sections, f_vector_t(N_t + 1, h));

        f_vector_2d_t voltageVec(N_t + 1, f_vector_t(n_sections, V));

        f_vector_2d_t dphiVec(n_sections, f_vector_t(N_t + 1, dphi));

        Context::GP = new GeneralParameters(N_t, CVec, alphaVec, alpha_order,
                                            momentumVec, proton);

        Context::Beam = new Beams(N_p, N_b);

        Context::RfP = new RfParameters(n_sections, hVec, voltageVec, dphiVec);

        longitudinal_bigaussian(tau_0 / 4, 0, -1, false);
        Context::Beam->statistics();
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
    const long long p_i = 450e9;     // Synchronous momentum [eV/c]
    const ftype p_f = 460.005e9;     // Synchronous momentum, final
    const long long h = 35640;       // Harmonic number
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

TEST_F(testBeam, test_sigma_dE) {

    std::vector<ftype> v;

    util::read_vector_from_file(v, statistics_params + "sigma_dE");
    ftype ref = v[0];
    ftype real = Context::Beam->sigma_dE;
    ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)));
}

TEST_F(testBeam, test_sigma_dt) {
    std::vector<ftype> v;

    util::read_vector_from_file(v, statistics_params + "sigma_dt");
    ftype ref = v[0];
    ftype real = Context::Beam->sigma_dt;
    ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)));
}

TEST_F(testBeam, test_mean_dE) {
    std::vector<ftype> v;

    util::read_vector_from_file(v, statistics_params + "mean_dE");
    ftype ref = v[0];
    ftype real = Context::Beam->mean_dE;
    ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)));
}

TEST_F(testBeam, test_mean_dt) {
    std::vector<ftype> v;

    util::read_vector_from_file(v, statistics_params + "mean_dt");
    ftype ref = v[0];
    ftype real = Context::Beam->mean_dt;
    ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)));
}

TEST_F(testBeam, test_epsn_rms_l) {
    std::vector<ftype> v;

    util::read_vector_from_file(v, statistics_params + "epsn_rms_l");
    ftype ref = v[0];
    ftype real = Context::Beam->epsn_rms_l;
    ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)));
}

TEST_F(testBeam, test_macroparticles_lost) {
    std::vector<ftype> v;

    util::read_vector_from_file(v, statistics_params + "n_macroparticles_lost");
    ftype ref = v[0];
    ftype real = Context::Beam->n_macroparticles_lost;
    ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)));
}

TEST_F(testBeam, test_losses_long_cut) {
    auto Beam = Context::Beam;
    Beam->losses_longitudinal_cut(Beam->dt.data(), Beam->mean_dt,
                                  10 * fabs(Beam->mean_dt), Beam->id.data());

    std::vector<ftype> v;
    util::read_vector_from_file(v, long_cut_params + "id");
    for (unsigned int i = 0; i < v.size(); ++i) {
        ftype ref = v[i];
        ftype real = Beam->id[i];
        ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)));
    }
}

TEST_F(testBeam, test_losses_energy_cut) {
    auto Beam = Context::Beam;

    Beam->losses_longitudinal_cut(Beam->dE.data(), Beam->mean_dE,
                                  10 * fabs(Beam->mean_dE), Beam->id.data());

    std::vector<ftype> v;
    util::read_vector_from_file(v, energy_cut_params + "id");
    for (unsigned int i = 0; i < v.size(); ++i) {
        ftype ref = v[i];
        ftype real = Beam->id[i];
        ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)));
    }
}

int main(int ac, char* av[]) {
    ::testing::InitGoogleTest(&ac, av);
    return RUN_ALL_TESTS();
}

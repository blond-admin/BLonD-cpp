#include <blond/globals.h>
#include <blond/llrf/PhaseNoise.h>
#include <blond/math_functions.h>
#include <blond/utilities.h>
#include <gtest/gtest.h>
#include <stdio.h>

// Simulation parameters
// --------------------------------------------------------

// Bunch parameters
const long long int N_b = (long int)1e9; // Intensity
const ftype tau_0 = 0.4e-9;         // Initial bunch length, 4 sigma [s]

// Machine and RF parameters
const ftype C = 26658.883;                   // Machine circumference [m]
const ftype p_i = 450e9;                     // Synchronous momentum [eV/c]
const long long h = 35640;                        // Harmonic number
const ftype V = 6e6;                         // RF voltage [V]
const ftype dphi = 0;                        // Phase modulation/offset
const ftype gamma_t = 55.759505;             // Transition gamma
const ftype alpha = 1.0 / gamma_t / gamma_t; // First order mom. comp. factor
const int alpha_order = 1;
const int n_sections = 1;
// Tracking details

unsigned N_t = 5000;  // Number of turns to track
unsigned N_p = 10001; // Macro-particles

unsigned N_slices = 1 << 8; // = (2^8)

class testPSBPhaseNoiseInjection : public ::testing::Test {

  protected:
    virtual void SetUp() {
        f_vector_2d_t momentumVec(n_sections, f_vector_t(N_t + 1));
        for (auto& v : momentumVec)
            mymath::linspace(v.data(), p_i, 1.01 * p_i, N_t + 1);

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
};

TEST_F(testPSBPhaseNoiseInjection, constructor1) {

    auto psbNoise = new PSBPhaseNoiseInjection();

    auto params = std::string(TEST_FILES"/PhaseNoise/") +
                  "PSBPhaseNoiseInjection/constructor1/";
    f_vector_t v;

    util::read_vector_from_file(v, params + "fs.txt");
    // util::dump(psbNoise->fFs, "fs\n");
    // ASSERT_EQ(v.size(), res.size());

    auto epsilon = 1e-8;
    for (unsigned int i = 0; i < v.size(); ++i) {
        auto ref = v[i];
        auto real = psbNoise->fFs[i];

        ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
            << "Testing of fFs failed on i " << i << std::endl;
    }

    delete psbNoise;
}

/*
TEST_F(testPSBPhaseNoiseInjection, spectrum_to_phase_noise1)
{

   auto psbNoise = new PSBPhaseNoiseInjection();

   auto params = std::string(TEST_FILES"/PhaseNoise/")
                 + "PSBPhaseNoiseInjection/spectrum_to_phase_noise1/";
   f_vector_t v;

   auto f = mymath::arange(0, 5.6e03, 1.12);
   f_vector_t spectrum = (f.size() / 10, 1);
   spectrum.resize(f.size(), 0);

   f_vector_t t, dphi;
   psbNoise->spectrum_to_phase_noise(t, dphi, f, spectrum);

   util::read_vector_from_file(v, params + "fs.txt");
   //util::dump(psbNoise->fFs, "fs\n");
   //ASSERT_EQ(v.size(), res.size());

   auto epsilon = 1e-8;
   for (unsigned int i = 0; i < v.size(); ++i) {
      auto ref = v[i];
      auto real = psbNoise->fFs[i];

      ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
            << "Testing of fFs failed on i "
            << i << std::endl;
   }

   delete psbNoise;
}
*/

TEST_F(testPSBPhaseNoiseInjection, generate_lin1) {

    auto psbNoise = new PSBPhaseNoiseInjection(
        2, 10, 0.8, 1.2, 2e-6, 12, 75, PhaseNoise::predistortion_t::linear);
    psbNoise->generate();

    auto params = std::string(TEST_FILES"/PhaseNoise/") +
                  "PSBPhaseNoiseInjection/generate/linear1/";
    f_vector_t v;

    util::read_vector_from_file(v, params + "dphi.txt");

    auto epsilon = 0.1;

    // auto i = mymath::max(v.data(), v.size());
    // auto maxV = std::abs(v[i]);
    // i = mymath::min(v.data(), v.size());
    // maxV = std::abs(v[i]) > maxV ? std::abs(v[i]) : maxV;

    // maxV = maxV > 1 ? maxV : 1 / maxV;
    auto meanV = mymath::mean(v.data(), v.size());

    auto stdV = mymath::standard_deviation(v.data(), v.size(), meanV);

    // std::cout << "ref mean = " << meanV << "\n";
    // std::cout << "ref std = " << stdV << "\n";

    auto real = psbNoise->fDphi;

    auto meanR = mymath::mean(real.data(), real.size());
    auto stdR = mymath::standard_deviation(real.data(), real.size(), meanR);
    // std::cout << "real mean = " << meanR << "\n";
    // std::cout << "real std = " << stdR << "\n";

    // ASSERT_NEAR(meanV, meanR, maxV * epsilon * std::min(fabs(meanV),
    // fabs(meanR)));

    ASSERT_NEAR(stdV, stdR, epsilon * std::min(fabs(stdV), fabs(stdR)));

    delete psbNoise;
}

TEST_F(testPSBPhaseNoiseInjection, generate_none1) {

    auto psbNoise = new PSBPhaseNoiseInjection(2, 10, 0.8, 1.2, 2e-6, 12, 75);
    psbNoise->generate();

    auto params = std::string(TEST_FILES"/PhaseNoise/") +
                  "PSBPhaseNoiseInjection/generate/none1/";
    f_vector_t v;

    util::read_vector_from_file(v, params + "dphi.txt");

    auto epsilon = 0.1;

    // auto i = mymath::max(v.data(), v.size());
    // auto maxV = std::abs(v[i]);
    // i = mymath::min(v.data(), v.size());
    // maxV = std::abs(v[i]) > maxV ? std::abs(v[i]) : maxV;

    // maxV = maxV > 1 ? maxV : 1 / maxV;
    // auto meanV = mymath::mean(v.data(), v.size());

    auto stdV = mymath::standard_deviation(v.data(), v.size());

    // std::cout << "ref mean = " << meanV << "\n";
    // std::cout << "ref std = " << stdV << "\n";

    auto real = psbNoise->fDphi;

    // auto meanR = mymath::mean(real.data(), real.size());
    auto stdR = mymath::standard_deviation(real.data(), real.size());
    // std::cout << "real mean = " << meanR << "\n";
    // std::cout << "real std = " << stdR << "\n";

    // ASSERT_NEAR(meanV, meanR, maxV * epsilon * std::min(fabs(meanV),
    // fabs(meanR)));

    ASSERT_NEAR(stdV, stdR, epsilon * std::min(fabs(stdV), fabs(stdR)));

    delete psbNoise;
}

int main(int ac, char* av[]) {
    ::testing::InitGoogleTest(&ac, av);
    return RUN_ALL_TESTS();
}

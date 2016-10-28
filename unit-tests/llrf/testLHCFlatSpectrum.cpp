#include <blond/globals.h>
#include <blond/llrf/PhaseNoise.h>
#include <blond/math_functions.h>
#include <blond/utilities.h>
#include <gtest/gtest.h>
#include <stdio.h>

// Simulation parameters
// --------------------------------------------------------

// Bunch parameters
const long long int N_b = 1e9; // Intensity
const double tau_0 = 0.4e-9;              // Initial bunch length, 4 sigma [s]

// Machine and RF parameters
const double C = 26658.883;                   // Machine circumference [m]
const double p_i = 450e9;                     // Synchronous momentum [eV/c]
const long long h = 35640;                   // Harmonic number
const double V = 6e6;                         // RF voltage [V]
const double dphi = 0;                        // Phase modulation/offset
const double gamma_t = 55.759505;             // Transition gamma
const double alpha = 1.0 / gamma_t / gamma_t; // First order mom. comp. factor
const int alpha_order = 1;
const int n_sections = 1;
// Tracking details

unsigned N_t = 1000;  // Number of turns to track
unsigned N_p = 10001; // Macro-particles

unsigned N_slices = 1 << 8; // = (2^8)

class testLHCFlatSpectrum : public ::testing::Test {

  protected:
    virtual void SetUp() {
        omp_set_num_threads(1);
        
        f_vector_2d_t momentumVec(n_sections, f_vector_t(N_t + 1));
        for (auto& v : momentumVec)
            mymath::linspace(v.data(), p_i, 1.01 * p_i, N_t + 1);

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

    }

    virtual void TearDown() {
        // Code here will be called immediately after each test
        // (right before the destructor).
        delete Context::GP;
        delete Context::Beam;
        delete Context::RfP;
    }
};

TEST_F(testLHCFlatSpectrum, constructor1) {

    auto lhcfs = new LHCFlatSpectrum(100, 1);

    auto params =
        std::string(TEST_FILES "/") + "LHCFlatSpectrum/constructor/test1/";
    f_vector_t v;

    util::read_vector_from_file(v, params + "fs.txt");
    // util::dump(lhcfs->fFs, "fs\n");
    // ASSERT_EQ(v.size(), res.size());

    auto epsilon = 1e-8;
    for (unsigned int i = 0; i < v.size(); ++i) {
        auto ref = v[i];
        auto real = lhcfs->fFs[i];

        ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
            << "Testing of fFs failed on i " << i << std::endl;
    }

    delete lhcfs;
}

TEST_F(testLHCFlatSpectrum, generate_exp1) {

    auto lhcfs =
        new LHCFlatSpectrum(1000, 10, 0.1, 1, 0.1, 1, 2,
                            LHCFlatSpectrum::predistortion_t::exponential);
    lhcfs->generate();

    auto params =
        std::string(TEST_FILES "/") + "LHCFlatSpectrum/generate/exponential/";
    f_vector_t v;

    util::read_vector_from_file(v, params + "dphi.txt");

    auto epsilon = 0.1;

    auto meanV = mymath::mean(v.data(), v.size());
    auto stdV = mymath::standard_deviation(v.data(), v.size(), meanV);

    auto real = lhcfs->fDphi;

    auto meanR = mymath::mean(real.data(), real.size());
    auto stdR = mymath::standard_deviation(real.data(), real.size(), meanR);

    // ASSERT_NEAR(meanV, meanR, epsilon * std::min(fabs(meanV), fabs(meanR)));

    ASSERT_NEAR(stdV, stdR, epsilon * std::min(fabs(stdV), fabs(stdR)));

    delete lhcfs;
}

TEST_F(testLHCFlatSpectrum, generate_lin1) {

    auto lhcfs = new LHCFlatSpectrum(1000, 10, 0.1, 1, 0.1, 1, 2,
                                     LHCFlatSpectrum::predistortion_t::linear);
    lhcfs->generate();

    auto params =
        std::string(TEST_FILES "/") + "LHCFlatSpectrum/generate/linear/";
    f_vector_t v;

    util::read_vector_from_file(v, params + "dphi.txt");

    auto epsilon = 0.1;

    auto meanV = mymath::mean(v.data(), v.size());
    auto stdV = mymath::standard_deviation(v.data(), v.size(), meanV);

    auto real = lhcfs->fDphi;

    auto meanR = mymath::mean(real.data(), real.size());
    auto stdR = mymath::standard_deviation(real.data(), real.size(), meanR);

    // ASSERT_NEAR(meanV, meanR, epsilon * std::min(fabs(meanV), fabs(meanR)));

    ASSERT_NEAR(stdV, stdR, epsilon * std::min(fabs(stdV), fabs(stdR)));

    delete lhcfs;
}

TEST_F(testLHCFlatSpectrum, generate_weight1) {

    auto lhcfs =
        new LHCFlatSpectrum(1000, 10, 0.1, 1, 0.1, 1, 2,
                            LHCFlatSpectrum::predistortion_t::weightfunction);
    lhcfs->generate();

    auto params = std::string(TEST_FILES "/") +
                  "LHCFlatSpectrum/generate/weightfunction/";
    f_vector_t v;

    util::read_vector_from_file(v, params + "dphi.txt");

    auto epsilon = 0.1;

    auto meanV = mymath::mean(v.data(), v.size());
    auto stdV = mymath::standard_deviation(v.data(), v.size(), meanV);

    auto real = lhcfs->fDphi;

    auto meanR = mymath::mean(real.data(), real.size());
    auto stdR = mymath::standard_deviation(real.data(), real.size(), meanR);

    // ASSERT_NEAR(meanV, meanR, epsilon * std::min(fabs(meanV), fabs(meanR)));

    ASSERT_NEAR(stdV, stdR, epsilon * std::min(fabs(stdV), fabs(stdR)));

    delete lhcfs;
}

TEST_F(testLHCFlatSpectrum, generate_hyper1) {

    auto lhcfs =
        new LHCFlatSpectrum(1000, 10, 0.1, 1, 0.1, 1, 2,
                            LHCFlatSpectrum::predistortion_t::hyperbolic);
    lhcfs->generate();

    auto params =
        std::string(TEST_FILES "/") + "LHCFlatSpectrum/generate/hyperbolic/";
    f_vector_t v;

    util::read_vector_from_file(v, params + "dphi.txt");

    auto epsilon = 0.1;

    auto meanV = mymath::mean(v.data(), v.size());
    auto stdV = mymath::standard_deviation(v.data(), v.size(), meanV);

    auto real = lhcfs->fDphi;

    auto meanR = mymath::mean(real.data(), real.size());
    auto stdR = mymath::standard_deviation(real.data(), real.size(), meanR);

    // ASSERT_NEAR(meanV, meanR, epsilon * std::min(fabs(meanV), fabs(meanR)));

    ASSERT_NEAR(stdV, stdR, epsilon * std::min(fabs(stdV), fabs(stdR)));

    delete lhcfs;
}

TEST_F(testLHCFlatSpectrum, generate_none1) {

    auto lhcfs = new LHCFlatSpectrum(1000, 10, 0.1, 1, 0.1, 1, 2);
    lhcfs->generate();

    auto params =
        std::string(TEST_FILES "/") + "LHCFlatSpectrum/generate/none/";
    f_vector_t v;

    util::read_vector_from_file(v, params + "dphi.txt");

    auto epsilon = 0.05;

    auto meanV = mymath::mean(v.data(), v.size());
    auto stdV = mymath::standard_deviation(v.data(), v.size(), meanV);

    auto real = lhcfs->fDphi;

    auto meanR = mymath::mean(real.data(), real.size());
    auto stdR = mymath::standard_deviation(real.data(), real.size(), meanR);

    // ASSERT_NEAR(meanV, meanR, epsilon * std::min(fabs(meanV), fabs(meanR)));

    ASSERT_NEAR(stdV, stdR, epsilon * std::min(fabs(stdV), fabs(stdR)));

    delete lhcfs;
}

int main(int ac, char* av[]) {
    ::testing::InitGoogleTest(&ac, av);
    return RUN_ALL_TESTS();
}

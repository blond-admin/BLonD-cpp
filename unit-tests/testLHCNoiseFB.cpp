#include <blond/beams/Distributions.h>
#include <blond/globals.h>
#include <blond/llrf/LHCNoiseFB.h>
#include <blond/math_functions.h>
#include <blond/trackers/Tracker.h>
#include <blond/utilities.h>
#include <gtest/gtest.h>
#include <stdio.h>

// Simulation parameters
// --------------------------------------------------------
const long long N_b = 1e9;      // Intensity
const ftype tau_0 = 0.4e-9; // Initial bunch length, 4 sigma [s]

// Machine and RF parameters
const ftype C = 26658.883;                   // Machine circumference [m]
const ftype p_i = 450e9;                     // Synchronous momentum [eV/c]
const long long h = 35640;                   // Harmonic number
const ftype V = 6e6;                         // RF voltage [V]
const ftype dphi = 0;                        // Phase modulation/offset
const ftype gamma_t = 55.759505;             // Transition gamma
const ftype alpha = 1.0 / gamma_t / gamma_t; // First order mom. comp. factor
const int alpha_order = 1;
const int n_sections = 1;

class testLHCNoiseFB : public ::testing::Test {
  protected:
    unsigned N_t = 1000;     // Number of turns to track
    unsigned N_p = 10001;    // Macro-particles
    unsigned N_slices = 100; // = (2^8)

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

        Context::GP = new GeneralParameters(N_t, CVec, alphaVec, alpha_order,
                                            momentumVec, proton);

        Context::Beam = new Beams(N_p, N_b);

        Context::RfP = new RfParameters(n_sections, hVec, voltageVec, dphiVec);

        // RingAndRfSection *long_tracker = new RingAndRfSection();

        longitudinal_bigaussian(tau_0 / 4, 0, -1, false);
        Context::Slice = new Slices(N_slices);
    }

    virtual void TearDown() {
        // Code here will be called immediately after each test
        // (right before the destructor).
        delete Context::GP;
        delete Context::Beam;
        delete Context::RfP;
        delete Context::Slice;
    }
};

class testLHCNoiseFBMultiBunch : public ::testing::Test {
  protected:
    unsigned N_t = 1000;        // Number of turns to track
    unsigned N_p = 10001;       // Macro-particles
    unsigned N_slices = 1 << 8; // = (2^8)

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

        // RingAndRfSection *long_tracker = new RingAndRfSection();

        // longitudinal_bigaussian(tau_0 / 4, 0, -1, false);
        Context::Slice = new Slices(N_slices);
    }

    virtual void TearDown() {
        // Code here will be called immediately after each test
        // (right before the destructor).
        delete Context::GP;
        delete Context::Beam;
        delete Context::RfP;
        delete Context::Slice;
    }
};

TEST_F(testLHCNoiseFB, constructor1) {

    auto lhcnfb = new LHCNoiseFB(1.0);

    auto params = std::string(TEST_FILES) +
                  "/LHCNoiseFB/constructor/test1/";
    f_vector_t v;

    util::read_vector_from_file(v, params + "g.txt");
    // util::dump(lhcnfb->fG, "fG\n");
    ASSERT_EQ(v.size(), lhcnfb->fG.size());

    auto epsilon = 1e-8;
    for (unsigned int i = 0; i < v.size(); ++i) {
        auto ref = v[i];
        auto real = lhcnfb->fG[i];

        ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
            << "Testing of fG failed on i " << i << std::endl;
    }

    delete lhcnfb;
}

TEST_F(testLHCNoiseFB, constructor2) {

    auto lhcnfb = new LHCNoiseFB(1.0, 0.1, 0.9, 100, false);

    auto params = std::string(TEST_FILES) +
                  "/LHCNoiseFB/constructor/test2/";
    f_vector_t v;

    util::read_vector_from_file(v, params + "g.txt");
    // util::dump(lhcnfb->fG, "fG\n");
    ASSERT_EQ(v.size(), lhcnfb->fG.size());

    auto epsilon = 1e-8;
    for (unsigned int i = 0; i < v.size(); ++i) {
        auto ref = v[i];
        auto real = lhcnfb->fG[i];

        ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
            << "Testing of fG failed on i " << i << std::endl;
    }

    delete lhcnfb;
}

TEST_F(testLHCNoiseFB, constructor3) {

    f_vector_t a = {1, 2, 3};
    auto lhcnfb = new LHCNoiseFB(1.0, 0.1, 0.9, 100, false, a);

    auto params = std::string(TEST_FILES) +
                  "/LHCNoiseFB/constructor/test3/";
    f_vector_t v;

    util::read_vector_from_file(v, params + "bl_meas_bbb.txt");
    // util::dump(lhcnfb->fG, "fG\n");
    ASSERT_EQ(v.size(), lhcnfb->fBlMeasBBB.size());

    delete lhcnfb;
}

TEST_F(testLHCNoiseFB, fwhm_interpolation1) {
    auto Slice = Context::Slice;
    auto lhcnfb = new LHCNoiseFB(1.0);
    for (uint i = 0; i < Slice->n_slices; i++) {
        Slice->n_macroparticles[i] = 50 * (i % 4);
        Slice->bin_centers[i] = 1e8 * (i + 1) / Slice->n_slices;
    }

    auto params = std::string(TEST_FILES) +
                  "/LHCNoiseFB/fwhm_interpolation/test1/";
    f_vector_t v;

    util::read_vector_from_file(v, params + "return.txt");
    // util::dump(lhcnfb->fG, "fG\n");
    auto index = mymath::arange<uint>(10, 20);
    auto epsilon = 1e-8;
    auto ref = v[0];
    auto real = lhcnfb->fwhm_interpolation(index, 100);
    ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)));

    delete lhcnfb;
}

TEST_F(testLHCNoiseFB, fwhm_interpolation2) {
    auto Slice = Context::Slice;

    auto lhcnfb = new LHCNoiseFB(1.0);
    for (uint i = 0; i < Slice->n_slices; i++) {
        Slice->n_macroparticles[i] = 100 * i;
        Slice->bin_centers[i] = 1e10 * (i + 1) / Slice->n_slices;
    }

    auto params = std::string(TEST_FILES) +
                  "/LHCNoiseFB/fwhm_interpolation/test2/";
    f_vector_t v;

    util::read_vector_from_file(v, params + "return.txt");
    // util::dump(lhcnfb->fG, "fG\n");
    auto index = mymath::arange<uint>(0, 99);
    auto epsilon = 1e-8;
    auto ref = v[0];
    auto real = lhcnfb->fwhm_interpolation(index, 1);
    ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)));

    delete lhcnfb;
}

TEST_F(testLHCNoiseFB, fwhm_single_bunch1) {
    auto Slice = Context::Slice;

    auto lhcnfb = new LHCNoiseFB(1.0);
    for (uint i = 0; i < Slice->n_slices; i++) {
        Slice->n_macroparticles[i] = (Slice->n_slices - i);
        Slice->bin_centers[i] = 1e10 * (i + 1) / Slice->n_slices;
    }
    lhcnfb->fwhm_single_bunch();

    auto params = std::string(TEST_FILES) +
                  "/LHCNoiseFB/fwhm_single_bunch/test1/";
    f_vector_t v;

    util::read_vector_from_file(v, params + "bl_meas.txt");
    // util::dump(lhcnfb->fG, "fG\n");
    auto epsilon = 1e-8;
    auto ref = v[0];
    auto real = lhcnfb->fBlMeas;
    ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)));

    delete lhcnfb;
}

TEST_F(testLHCNoiseFBMultiBunch, fwhm_multi_bunch1) {
    auto Slice = Context::Slice;

    auto long_tracker = new RingAndRfSection();

    for (uint i = 0; i < Slice->n_slices; i++) {
        Slice->n_macroparticles[i] = (Slice->n_slices - i);
        Slice->bin_centers[i] = 1e8 * (i + 1) / Slice->n_slices;
    }
    // lhcnfb->fwhm_multi_bunch();

    f_vector_t a = {0, 10, 20};
    auto lhcnfb = new LHCNoiseFB(0.9e-9, 0.1e9, 0.9, 1, true, a);

    f_vector_t realV;
    for (uint i = 0; i < N_t; ++i) {
        long_tracker->track();
        Slice->track();
        lhcnfb->track();
        realV.push_back(lhcnfb->fX);
        // std::cout << "x: " << lhcnfb->fX << '\n';
    }

    auto params = std::string(TEST_FILES) +
                  "/LHCNoiseFB/fwhm_multi_bunch/test1/";
    f_vector_t v;

    util::read_vector_from_file(v, params + "x.txt");
    // util::dump(lhcnfb->fG, "fG\n");
    ASSERT_EQ(v.size(), realV.size());

    auto epsilon = 1e-8;
    for (uint i = 0; i < v.size(); ++i) {
        auto ref = v[i];
        auto real = realV[i];

        ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
            << "Testing of fX failed on i " << i << std::endl;
    }

    delete lhcnfb;
}

TEST_F(testLHCNoiseFBMultiBunch, fwhm_multi_bunch2) {
    auto Slice = Context::Slice;

    auto long_tracker = new RingAndRfSection();

    for (uint i = 0; i < Slice->n_slices; i++) {
        Slice->n_macroparticles[i] = (Slice->n_slices - i);
        Slice->bin_centers[i] = 1e6 * (i + 1) / Slice->n_slices;
    }
    // lhcnfb->fwhm_multi_bunch();

    f_vector_t a = {0, 10, 20, 30, 40};
    auto lhcnfb = new LHCNoiseFB(1e-9, 1e8, 0.5, 10, false, a);

    f_vector_t realV;
    for (uint i = 0; i < N_t; ++i) {
        long_tracker->track();
        Slice->track();
        lhcnfb->track();
        realV.push_back(lhcnfb->fX);
        // std::cout << "x: " << lhcnfb->fX << '\n';
    }

    auto params = std::string(TEST_FILES) +
                  "/LHCNoiseFB/fwhm_multi_bunch/test2/";
    f_vector_t v;

    util::read_vector_from_file(v, params + "x.txt");
    // util::dump(lhcnfb->fG, "fG\n");
    ASSERT_EQ(v.size(), realV.size());

    auto epsilon = 1e-8;
    for (uint i = 0; i < v.size(); ++i) {
        auto ref = v[i];
        auto real = realV[i];

        ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
            << "Testing of fX failed on i " << i << std::endl;
    }

    delete lhcnfb;
}

TEST_F(testLHCNoiseFB, track1) {

    auto long_tracker = new RingAndRfSection();

    auto lhcnfb = new LHCNoiseFB(1.0, 0.1, 0.9, 1);

    f_vector_t res;
    for (uint i = 0; i < 100; ++i) {
        long_tracker->track();
        Context::Slice->track();
        lhcnfb->track();
        res.push_back(lhcnfb->fX);
    }

    auto params =
        std::string(TEST_FILES) + "/LHCNoiseFB/track/test1/";
    f_vector_t v;

    util::read_vector_from_file(v, params + "x.txt");
    // util::dump(lhcnfb->fG, "fG\n");

    ASSERT_EQ(v.size(), res.size());

    auto epsilon = 1e-8;
    for (unsigned int i = 0; i < v.size(); ++i) {
        auto ref = v[i];
        auto real = res[i];

        ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
            << "Testing of fX failed on i " << i << std::endl;
    }

    delete lhcnfb;
    delete long_tracker;
}

int main(int ac, char* av[]) {
    ::testing::InitGoogleTest(&ac, av);
    return RUN_ALL_TESTS();
}

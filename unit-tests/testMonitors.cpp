#include <blond/monitors/Monitors.h>
#include <blond/constants.h>
#include <blond/beams/Distributions.h>
#include <blond/input_parameters/GeneralParameters.h>
#include <blond/math_functions.h>
#include <blond/utilities.h>
#include <gtest/gtest.h>


class testMonitors : public ::testing::Test {

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
    const int dt_save = 1000;
    const int N_t = 5000;            // Number of turns to track
    const int N_p = 10000;            // Macro-particles
    const int N_slices = 50;
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


TEST_F(testMonitors, SlicesMonitor1)
{
    std::string params =
        TEST_FILES "/Monitors/SlicesMonitor1/";
    auto filename = "n_macroparticles.h5";

    longitudinal_bigaussian(tau_0 / 4, 0, -1, false);
    auto slice = Slices(N_slices);
    auto slicesMonitor = SlicesMonitor(filename, N_t / dt_save + 1, &slice);
    auto tracker = RingAndRfSection();

    for (int i = 0; i < N_t; i++) {
        tracker.track();
        slice.track();
        if (i % dt_save == 0)
            slicesMonitor.track();
    }
    slicesMonitor.track();

    hsize_t dimsCpp[2];
    int *real = (int *) read_2D(filename, "Slices/n_macroparticles",
                                "int", dimsCpp);
    hsize_t dimsPy[2];
    float *ref = (float *) read_2D(params + filename, "Slices/n_macroparticles",
                                   "float", dimsPy);

    for (int i = 0; i < 2; i++) ASSERT_EQ(dimsCpp[i], dimsPy[i]);

    for (uint i = 0; i < dimsCpp[0] * dimsCpp[1]; i++)
        ASSERT_EQ(ref[i], real[i])
                << "Testing of n_macroparticles failed on i " << i << '\n';
    free(real);
    free(ref);
    std::remove(filename);
}


TEST_F(testMonitors, BunchMonitor1)
{
    auto GP = Context::GP;
    auto RfP = Context::RfP;
    auto beam = Context::Beam;
    auto epsilon = 1e-6;
    std::string params =
        TEST_FILES "/Monitors/BunchMonitor1/";
    auto filename = "bunch.h5";

    longitudinal_bigaussian(tau_0 / 4, 0, -1, false);
    auto slice = Slices(N_slices);
    auto bunchmonitor = BunchMonitor(GP, RfP, beam, filename, 100);
    auto tracker = RingAndRfSection();

    for (int i = 0; i < N_t; i++) {
        tracker.track();
        slice.track();
        bunchmonitor.track();
    }
    bunchmonitor.track();
    bunchmonitor.close();

    hsize_t dimsPy[1];
    hsize_t dimsCpp[1];

    int *realInt = (int *) read_1D(filename,
                                   "Beam/n_macroparticles_alive",
                                   "int", dimsCpp);
    float *ref = (float *) read_1D(params + filename,
                                   "Beam/n_macroparticles_alive",
                                   "float", dimsPy);
    ASSERT_EQ(dimsCpp[0], dimsPy[0]);
    for (uint i = 0; i < dimsCpp[0]; i++)
        ASSERT_EQ(ref[i], realInt[i])
                << "Testing of n_macroparticles failed on i " << i << '\n';

    free(ref); free(realInt);


    double *realD = (double *) read_1D(filename,
                                       "Beam/mean_dt",
                                       "double", dimsCpp);
    ref = (float *) read_1D(params + filename,
                            "Beam/mean_dt",
                            "float", dimsPy);
    ASSERT_EQ(dimsCpp[0], dimsPy[0]);
    for (uint i = 0; i < dimsCpp[0]; ++i)
        ASSERT_NEAR(ref[i], realD[i], epsilon *
                    std::max((double)std::abs(ref[i]), std::abs(realD[i])))
                << "Testing of mean_dt failed on i " << i << '\n';
    free(realD); free(ref);


    realD = (double *) read_1D(filename,
                               "Beam/mean_dE",
                               "double", dimsCpp);
    ref = (float *) read_1D(params + filename,
                            "Beam/mean_dE",
                            "float", dimsPy);
    ASSERT_EQ(dimsCpp[0], dimsPy[0]);
    for (uint i = 0; i < dimsCpp[0]; ++i)
        ASSERT_NEAR(ref[i], realD[i], epsilon *
                    std::max(std::abs((double)ref[i]), std::abs(realD[i])))
                << "Testing of mean_dE failed on i " << i << '\n';
    free(realD); free(ref);


    realD = (double *) read_1D(filename,
                               "Beam/sigma_dt",
                               "double", dimsCpp);
    ref = (float *) read_1D(params + filename,
                            "Beam/sigma_dt",
                            "float", dimsPy);
    ASSERT_EQ(dimsCpp[0], dimsPy[0]);
    for (uint i = 0; i < dimsCpp[0]; ++i)
        ASSERT_NEAR(ref[i], realD[i], epsilon *
                    std::max((double)std::abs(ref[i]), std::abs(realD[i])))
                << "Testing of sigma_dt failed on i " << i << '\n';
    free(realD); free(ref);

    realD = (double *) read_1D(filename,
                               "Beam/sigma_dE",
                               "double", dimsCpp);
    ref = (float *) read_1D(params + filename,
                            "Beam/sigma_dE",
                            "float", dimsPy);
    ASSERT_EQ(dimsCpp[0], dimsPy[0]);
    for (uint i = 0; i < dimsCpp[0]; ++i)
        ASSERT_NEAR(ref[i], realD[i], epsilon *
                    std::max((double)std::abs(ref[i]), std::abs(realD[i])))
                << "Testing of sigma_dE failed on i " << i << '\n';
    free(realD); free(ref);

    realD = (double *) read_1D(filename,
                               "Beam/epsn_rms_l",
                               "double", dimsCpp);
    ref = (float *) read_1D(params + filename,
                            "Beam/epsn_rms_l",
                            "float", dimsPy);
    ASSERT_EQ(dimsCpp[0], dimsPy[0]);
    for (uint i = 0; i < dimsCpp[0]; ++i)
        ASSERT_NEAR(ref[i], realD[i], epsilon *
                    std::max((double)std::abs(ref[i]), std::abs(realD[i])))
                << "Testing of epsn_rms_l failed on i " << i << '\n';
    free(realD); free(ref);


    std::remove(filename);
}


int main(int ac, char *av[])
{
    ::testing::InitGoogleTest(&ac, av);
    return RUN_ALL_TESTS();
}

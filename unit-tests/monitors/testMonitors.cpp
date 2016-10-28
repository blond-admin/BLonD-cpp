#include <blond/monitors/Monitors.h>
#include <blond/constants.h>
#include <blond/beams/Distributions.h>
#include <blond/input_parameters/GeneralParameters.h>
#include <blond/math_functions.h>
#include <blond/utilities.h>
#include <gtest/gtest.h>

using namespace std;

class testMonitors : public ::testing::Test {

protected:
    const long long N_b = 1e9;      // Intensity
    const double tau_0 = 0.4e-9;     // Initial bunch length, 4 sigma [s]
    // Machine and RF parameters
    const double C = 26658.883;       // Machine circumference [m]
    const double p_i = 450e9;     // Synchronous momentum [eV/c]
    const double p_f = 460.005e9;     // Synchronous momentum, final
    const long long h = 35640;       // Harmonic number
    const double V = 6e6;             // RF voltage [V]
    const double dphi = 0;            // Phase modulation/offset
    const double gamma_t = 55.759505; // Transition gamma
    const double alpha =
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
        omp_set_num_threads(1);

        f_vector_2d_t momentumVec(n_sections, f_vector_t(N_t + 1, 0));
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
    auto RfP = Context::RfP;
    auto Beam = Context::Beam;
    auto GP = Context::GP;
    
    string params =
        TEST_FILES "/Monitors/SlicesMonitor1/";
    auto filename = "n_macroparticles.h5";
    remove(filename);

    longitudinal_bigaussian(GP, RfP, Beam, tau_0 / 4, 0, -1, false);
    auto slice = Slices(RfP, Beam, N_slices);
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
    remove(filename);
}


TEST_F(testMonitors, BunchMonitor1)
{
    auto GP = Context::GP;
    auto RfP = Context::RfP;
    auto Beam = Context::Beam;
    auto epsilon = 1e-7;
    string params =
        TEST_FILES "/Monitors/BunchMonitor1/";
    auto filename = "bunch.h5";
    remove(filename);

    longitudinal_bigaussian(GP, RfP, Beam, tau_0 / 4, 0, -1, false);
    auto slice = Slices(RfP, Beam, N_slices);
    auto bunchmonitor = BunchMonitor(GP, RfP, Beam, filename, 100);
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
        ASSERT_EQ((int)ref[i], realInt[i])
                << "Testing of n_macroparticles failed on i " << i << '\n';

    free(ref); free(realInt);


    float *realD = (float *) read_1D(filename,
                                     "Beam/mean_dt",
                                     "float", dimsCpp);
    ref = (float *) read_1D(params + filename,
                            "Beam/mean_dt",
                            "float", dimsPy);
    ASSERT_EQ(dimsCpp[0], dimsPy[0]);
    for (uint i = 0; i < dimsCpp[0]; ++i)
        ASSERT_NEAR(ref[i], realD[i], epsilon *
                    max(abs(ref[i]), abs(realD[i])))
                << "Testing of mean_dt failed on i " << i << '\n';
    free(realD); free(ref);


    realD = (float *) read_1D(filename,
                              "Beam/mean_dE",
                              "float", dimsCpp);
    ref = (float *) read_1D(params + filename,
                            "Beam/mean_dE",
                            "float", dimsPy);
    ASSERT_EQ(dimsCpp[0], dimsPy[0]);

    for (uint i = 0; i < dimsCpp[0]; ++i)
        ASSERT_NEAR(ref[i], realD[i], epsilon *
                    max(abs(ref[i]), abs(realD[i])))
                << "Testing of mean_dE failed on i " << i << '\n';
    free(realD); free(ref);


    realD = (float *) read_1D(filename,
                              "Beam/sigma_dt",
                              "float", dimsCpp);
    ref = (float *) read_1D(params + filename,
                            "Beam/sigma_dt",
                            "float", dimsPy);
    ASSERT_EQ(dimsCpp[0], dimsPy[0]);
    for (uint i = 0; i < dimsCpp[0]; ++i)
        ASSERT_NEAR(ref[i], realD[i], epsilon *
                    max(abs(ref[i]), abs(realD[i])))
                << "Testing of sigma_dt failed on i " << i << '\n';
    free(realD); free(ref);

    realD = (float *) read_1D(filename,
                              "Beam/sigma_dE",
                              "float", dimsCpp);
    ref = (float *) read_1D(params + filename,
                            "Beam/sigma_dE",
                            "float", dimsPy);
    ASSERT_EQ(dimsCpp[0], dimsPy[0]);
    for (uint i = 0; i < dimsCpp[0]; ++i)
        ASSERT_NEAR(ref[i], realD[i], epsilon *
                    max(abs(ref[i]), abs(realD[i])))
                << "Testing of sigma_dE failed on i " << i << '\n';
    free(realD); free(ref);

    realD = (float *) read_1D(filename,
                              "Beam/epsn_rms_l",
                              "float", dimsCpp);
    ref = (float *) read_1D(params + filename,
                            "Beam/epsn_rms_l",
                            "float", dimsPy);
    ASSERT_EQ(dimsCpp[0], dimsPy[0]);
    for (uint i = 0; i < dimsCpp[0]; ++i)
        ASSERT_NEAR(ref[i], realD[i], epsilon *
                    max(abs(ref[i]), abs(realD[i])))
                << "Testing of epsn_rms_l failed on i " << i << '\n';
    free(realD); free(ref);


    remove(filename);
}


TEST_F(testMonitors, BunchMonitor2)
{

    auto GP = Context::GP;
    auto RfP = Context::RfP;
    auto Beam = Context::Beam;
    auto epsilon = 1e-7;
    string params =
        TEST_FILES "/Monitors/BunchMonitor2/";
    auto filename = "bunch.h5";
    remove(filename);

    longitudinal_bigaussian(GP, RfP, Beam, tau_0 / 4, 0, -1, false);
    auto slice = Slices(RfP, Beam, N_slices);
    auto PL_gain = 1.0 / (5 * GP->t_rev[0]);
    auto SL_gain = PL_gain / 10.0;
    f_vector_t PL_gainVec(N_t + 1 , PL_gain);
    auto PL = LHC(PL_gainVec, SL_gain);
    auto bunchmonitor = BunchMonitor(GP, RfP, Beam, filename,
                                     100, &slice, &PL);
    auto tracker = RingAndRfSection();

    for (int i = 0; i < N_t; i++) {
        tracker.track();
        slice.track();
        // cout << "round " << i << "\n";
        bunchmonitor.track();
    }
    bunchmonitor.track();
    bunchmonitor.close();

    hsize_t dimsPy[1];
    hsize_t dimsCpp[1];

    double *realD = (double *) read_1D(filename,
                                       "Beam/PL_omegaRF",
                                       "double", dimsCpp);
    double *refD = (double *) read_1D(params + filename,
                                      "Beam/PL_omegaRF",
                                      "double", dimsPy);
    ASSERT_EQ(dimsCpp[0], dimsPy[0]);

    for (uint i = 0; i < dimsCpp[0]; i++)
        ASSERT_NEAR(refD[i], realD[i], epsilon *
                    max(abs(refD[i]), abs(realD[i])))
                << "Testing of PL_omegaRF failed on i " << i << '\n';

    free(refD); free(realD);

    float *real = (float *) read_1D(filename,
                                    "Beam/PL_phiRF",
                                    "float", dimsCpp);
    float *ref = (float *) read_1D(params + filename,
                                   "Beam/PL_phiRF",
                                   "float", dimsPy);
    ASSERT_EQ(dimsCpp[0], dimsPy[0]);
    for (uint i = 0; i < dimsCpp[0]; ++i)
        ASSERT_NEAR(ref[i], real[i], epsilon *
                    max(abs(ref[i]), abs(real[i])))
                << "Testing of PL_phiRF failed on i " << i << '\n';
    free(real); free(ref);


    real = (float *) read_1D(filename,
                             "Beam/PL_bunch_phase",
                             "float", dimsCpp);
    ref = (float *) read_1D(params + filename,
                            "Beam/PL_bunch_phase",
                            "float", dimsPy);
    ASSERT_EQ(dimsCpp[0], dimsPy[0]);
    for (uint i = 0; i < dimsCpp[0]; ++i)
        ASSERT_NEAR(ref[i], real[i], epsilon *
                    max(abs(ref[i]), abs(real[i])))
                << "Testing of PL_bunch_phase failed on i " << i << '\n';
    free(real); free(ref);


    real = (float *) read_1D(filename,
                             "Beam/PL_phase_corr",
                             "float", dimsCpp);
    ref = (float *) read_1D(params + filename,
                            "Beam/PL_phase_corr",
                            "float", dimsPy);
    ASSERT_EQ(dimsCpp[0], dimsPy[0]);
    for (uint i = 0; i < dimsCpp[0]; ++i)
        ASSERT_NEAR(ref[i], real[i], epsilon *
                    max(abs(ref[i]), abs(real[i])))
                << "Testing of PL_phase_corr failed on i " << i << '\n';
    free(real); free(ref);

    real = (float *) read_1D(filename,
                             "Beam/PL_omegaRF_corr",
                             "float", dimsCpp);
    ref = (float *) read_1D(params + filename,
                            "Beam/PL_omegaRF_corr",
                            "float", dimsPy);
    ASSERT_EQ(dimsCpp[0], dimsPy[0]);
    for (uint i = 0; i < dimsCpp[0]; ++i)
        ASSERT_NEAR(ref[i], real[i], epsilon *
                    max(abs(ref[i]), abs(real[i])))
                << "Testing of PL_omegaRF_corr failed on i " << i << '\n';
    free(real); free(ref);

    real = (float *) read_1D(filename,
                             "Beam/SL_dphiRF",
                             "float", dimsCpp);
    ref = (float *) read_1D(params + filename,
                            "Beam/SL_dphiRF",
                            "float", dimsPy);
    ASSERT_EQ(dimsCpp[0], dimsPy[0]);
    for (uint i = 0; i < dimsCpp[0]; ++i)
        ASSERT_NEAR(ref[i], real[i], epsilon *
                    max(abs(ref[i]), abs(real[i])))
                << "Testing of SL_dphiRF failed on i " << i << '\n';
    free(real); free(ref);


    real = (float *) read_1D(filename,
                             "Beam/RL_drho",
                             "float", dimsCpp);
    ref = (float *) read_1D(params + filename,
                            "Beam/RL_drho",
                            "float", dimsPy);
    ASSERT_EQ(dimsCpp[0], dimsPy[0]);
    for (uint i = 0; i < dimsCpp[0]; ++i)
        ASSERT_NEAR(ref[i], real[i], epsilon *
                    max(abs(ref[i]), abs(real[i])))
                << "Testing of RL_drho failed on i " << i << '\n';
    free(real); free(ref);

    remove(filename);
}


// TODO: test the LHC_noise_FB_bl_bbb dataset
TEST_F(testMonitors, BunchMonitor3)
{
    auto GP = Context::GP;
    auto RfP = Context::RfP;
    auto Beam = Context::Beam;
    auto epsilon = 1e-7;
    string params =
        TEST_FILES "/Monitors/BunchMonitor3/";
    auto filename = "bunch.h5";

    remove(filename);

    longitudinal_bigaussian(GP, RfP, Beam, tau_0 / 4, 0, -1, false);
    auto slice = Slices(RfP, Beam, N_slices);
    auto PL_gain = 1.0 / (5 * GP->t_rev[0]);
    auto SL_gain = PL_gain / 10.0;
    f_vector_t PL_gainVec(N_t + 1 , PL_gain);
    auto PL = LHC(PL_gainVec, SL_gain);

    auto noiseFB = LHCNoiseFB(1e-8, 0.1e9, 0.93, 100, true, {10, 20, 30});

    auto bunchmonitor = BunchMonitor(GP, RfP, Beam, filename,
                                     100, &slice, &PL, &noiseFB);
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

    float *real = (float *) read_1D(filename,
                                    "Beam/LHC_noise_FB_factor",
                                    "float", dimsCpp);
    float *ref = (float *) read_1D(params + filename,
                                   "Beam/LHC_noise_FB_factor",
                                   "float", dimsPy);
    ASSERT_EQ(dimsCpp[0], dimsPy[0]);
    for (uint i = 0; i < dimsCpp[0]; i++)
        ASSERT_NEAR(ref[i], real[i], epsilon *
                    max(abs(ref[i]), abs(real[i])))
                << "Testing of LHC_noise_FB_factor failed on i " << i << '\n';

    free(ref); free(real);


    real = (float *) read_1D(filename,
                             "Beam/LHC_noise_FB_bl",
                             "float", dimsCpp);
    ref = (float *) read_1D(params + filename,
                            "Beam/LHC_noise_FB_bl",
                            "float", dimsPy);
    ASSERT_EQ(dimsCpp[0], dimsPy[0]);
    for (uint i = 0; i < dimsCpp[0]; ++i)
        ASSERT_NEAR(ref[i], real[i], epsilon *
                    max(abs(ref[i]), abs(real[i])))
                << "Testing of LHC_noise_FB_bl failed on i " << i << '\n';
    free(real); free(ref);


    // hsize_t dimsCpp2[2];
    // hsize_t dimsPy2[2];
    // real = (float *) read_2D(filename, "Beam/LHC_noise_FB_bl_bbb",
    //                           "float", dimsCpp2);
    // ref = (float *) read_2D(params + filename, "Beam/LHC_noise_FB_bl_bbb",
    //                         "float", dimsPy2);

    // for (int i = 0; i < 2; i++) ASSERT_EQ(dimsCpp2[i], dimsPy2[i]);

    // for (uint i = 0; i < dimsCpp2[0] * dimsCpp2[1]; i++)
    //     ASSERT_NEAR(ref[i], real[i], epsilon *
    //                 max((float)abs(ref[i]), abs(real[i])))
    //             << "Testing of LHC_noise_FB_bl_bbb failed on i " << i << '\n';


    // free(real); free(ref);
    remove(filename);
}

int main(int ac, char *av[])
{
    ::testing::InitGoogleTest(&ac, av);
    return RUN_ALL_TESTS();
}

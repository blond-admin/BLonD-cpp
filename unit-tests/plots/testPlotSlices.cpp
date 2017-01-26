#include <blond/blond.h>

#include <gtest/gtest.h>

using namespace std;

class testPlotSlices : public ::testing::Test {

protected:
    // Machine and RF parameters
    const double C = 26658.883;       // Machine circumference [m]
    const double p_i = 450e9;     // Synchronous momentum [eV/c]
    const double p_f = 460.005e9;     // Synchronous momentum, final
    const long long h = 35640;       // Harmonic number
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
        omp_set_num_threads(1);


        f_vector_2d_t momentumVec(n_sections, f_vector_t(N_t + 1));
        for (auto &v : momentumVec)
            mymath::linspace(v.data(), p_i, p_f, N_t + 1);

        f_vector_2d_t alphaVec(n_sections, f_vector_t(alpha_order, alpha));

        f_vector_t CVec(n_sections, C);

        f_vector_2d_t hVec(n_sections, f_vector_t(N_t + 1, h));

        f_vector_2d_t voltageVec(n_sections, f_vector_t(N_t + 1, V));

        f_vector_2d_t dphiVec(n_sections, f_vector_t(N_t + 1, dphi));

        Context::GP = new GeneralParameters(N_t, CVec, alphaVec,
                                             momentumVec,
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


TEST_F(testPlotSlices, plot_beam_profile1)
{
    auto GP = Context::GP;
    auto RfP = Context::RfP;
    auto Beam = Context::Beam;
    longitudinal_bigaussian(GP, RfP, Beam, tau_0 / 4, 0, -1, false);

    auto slice = new Slices(RfP, Beam, N_slices);

    ASSERT_EQ(plot_beam_profile(slice, 0), 1);

    delete slice;
}


TEST_F(testPlotSlices, plot_beam_profile_derivative1)
{
    auto GP = Context::GP;
    auto RfP = Context::RfP;
    auto Beam = Context::Beam;
    longitudinal_bigaussian(GP, RfP, Beam, tau_0 / 4, 0, -1, false);

    auto slice = new Slices(RfP, Beam, N_slices);

    ASSERT_EQ(plot_beam_profile_derivative(slice, 0), 1);

    delete slice;
}


TEST_F(testPlotSlices, plot_beam_profile_derivative2)
{
    auto GP = Context::GP;
    auto RfP = Context::RfP;
    auto Beam = Context::Beam;
    longitudinal_bigaussian(GP, RfP, Beam, tau_0 / 4, 0, -1, false);
    auto slice = new Slices(RfP, Beam, N_slices);

    ASSERT_EQ(plot_beam_profile_derivative(slice, 0, "-", "fig",
    {"filter1d", "gradient", "diff"}), 1);

    delete slice;
}


TEST_F(testPlotSlices, plot_beam_spectrum1)
{


    auto GP = Context::GP;
    auto RfP = Context::RfP;
    auto Beam = Context::Beam;
    longitudinal_bigaussian(GP, RfP, Beam, tau_0 / 4, 0, -1, false);
    auto slice = new Slices(RfP, Beam, N_slices);

    ASSERT_EQ(plot_beam_spectrum(slice, 0), 1);

    delete slice;
}



int main(int ac, char *av[])
{
    python::initialize();
    ::testing::InitGoogleTest(&ac, av);
    auto ret = RUN_ALL_TESTS();
    python::finalize();
    return ret;
}

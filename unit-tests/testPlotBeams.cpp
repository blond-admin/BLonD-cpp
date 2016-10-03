#include <iostream>
#include <blond/beams/Distributions.h>
#include <blond/input_parameters/GeneralParameters.h>
#include <blond/math_functions.h>
#include <blond/utilities.h>
#include <blond/plots/plot_beams.h>
#include <gtest/gtest.h>

using namespace std;

string h5file = TEST_FILES "/Plots/beam.h5";

class testPlotBeams : public ::testing::Test {

protected:
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
    const long long N_b = 1e9;  // Intensity
    const ftype tau_0 = 0.4e-9; // Initial bunch length, 4 sigma [s]

    const int N_t = 2000; // Number of turns to track
    const int N_p = 10000;  // Macro-particles
    const int N_slices = 100;



    virtual void SetUp()
    {
        omp_set_num_threads(1);


        f_vector_2d_t momentumVec(n_sections, f_vector_t(N_t + 1));
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


TEST_F(testPlotBeams, plot_long_phase_space1)
{
    longitudinal_bigaussian(tau_0 / 4, 0, -1, false);
    auto GP = Context::GP;
    auto RfP = Context::RfP;
    auto Beam = Context::Beam;

    ASSERT_EQ(plot_long_phase_space(GP, RfP, Beam, 0., 1e-5, -1e7, 1e7), 1);

}


TEST_F(testPlotBeams, plot_long_phase_space2)
{
    longitudinal_bigaussian(tau_0 / 4, 0, -1, false);
    auto GP = Context::GP;
    auto RfP = Context::RfP;
    auto Beam = Context::Beam;

    ASSERT_EQ(plot_long_phase_space(GP, RfP, Beam, 0., 1e-5,
                                    -1e7, 1e7,  "s", 1, true), 1);

}

TEST_F(testPlotBeams, plot_bunch_length_evol1)
{
    auto GP = Context::GP;
    auto RfP = Context::RfP;
    auto Beam = Context::Beam;

    ASSERT_EQ(plot_bunch_length_evol(RfP, h5file), 1);

}


/*
TEST_F(testPlotBeams, plot_bunch_length_evol_gaussian1)
{
    auto GP = Context::GP;
    auto RfP = Context::RfP;
    auto Beam = Context::Beam;

    plot_bunch_length_evol_gaussian(RfP, h5file);

}
*/


TEST_F(testPlotBeams, plot_position_evol1)
{
    auto GP = Context::GP;
    auto RfP = Context::RfP;
    auto Beam = Context::Beam;

    ASSERT_EQ(plot_position_evol(RfP, h5file), 1);

}



TEST_F(testPlotBeams, plot_energy_evol1)
{
    auto GP = Context::GP;
    auto RfP = Context::RfP;
    auto Beam = Context::Beam;

    ASSERT_EQ(plot_energy_evol(RfP, h5file), 1);

}


TEST_F(testPlotBeams, plot_transmitted_particles1)
{
    auto GP = Context::GP;
    auto RfP = Context::RfP;
    auto Beam = Context::Beam;

    ASSERT_EQ(plot_transmitted_particles(RfP, h5file), 1);

}


int main(int ac, char *av[])
{
    python::initialize();
    ::testing::InitGoogleTest(&ac, av);
    auto ret = RUN_ALL_TESTS();
    python::finalize();
    return ret;
}

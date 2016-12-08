
#include <blond/globals.h>
#include <blond/input_parameters/GeneralParameters.h>
#include <blond/input_parameters/RfParameters.h>
#include <blond/impedances/Music.h>
#include <blond/utilities.h>
#include <testing_utilities.h>
#include <gtest/gtest.h>
#include <blond/beams/Distributions.h>
using namespace std;

class testMusic : public ::testing::Test {

protected:
    // Simulation parameters
    // --------------------------------------------------------
    // Bunch parameters
    const long int N_b = 1e12; // Intensity
    // Machine and RF parameters
    const double radius = 25;
    const double C = 2 * constant::pi * radius;   // Machine circumference [m]
    const double mass_rest = constant::m_p * constant::c * constant::c / constant::e;
    const double tot_energy = 13e9;
    // final
    const double h = 1;                     // Harmonic number
    const double V = 24e3;                        // RF voltage [V]
    const double dphi = 0;                         // Phase modulation/offset
    const double gamma_t = 4.076750841; // Transition gamma
    const double momentum = sqrt(tot_energy * tot_energy - mass_rest *mass_rest);
    const double alpha = 1.0 / gamma_t / gamma_t; // First order mom. comp. factor
    const int alpha_order = 1;
    const int n_rf_systems = 1;
    // Tracking details

    long N_t = 1;       // Number of turns to track
    long N_p = 10000; // Macro-particles
    long N_slices = 1000;

    virtual void SetUp()
    {

        omp_set_num_threads(1);

        f_vector_2d_t momentumVec(n_rf_systems, f_vector_t(N_t + 1, momentum));

        f_vector_2d_t alphaVec(n_rf_systems, f_vector_t(alpha_order, alpha));

        f_vector_t CVec(n_rf_systems, C);

        f_vector_2d_t hVec(n_rf_systems, f_vector_t(N_t + 1, h));

        f_vector_2d_t voltageVec(n_rf_systems, f_vector_t(N_t + 1, V));

        f_vector_2d_t dphiVec(n_rf_systems, f_vector_t(N_t + 1, dphi));

        Context::GP = new GeneralParameters(N_t, CVec, alphaVec,
                                             momentumVec,
                                            GeneralParameters::particle_t::proton);
        auto GP = Context::GP;

        auto Beam = Context::Beam = new Beams(GP, N_p, N_b);

        auto RfP = Context::RfP = new RfParameters(GP, n_rf_systems, hVec,
                voltageVec, dphiVec);



        longitudinal_bigaussian(GP, RfP, Beam, 1e-9, 0, -1);
        Context::Slice = new Slices(RfP, Beam, N_slices);
    }

    virtual void TearDown()
    {
        // Code here will be called immediately after each test
        // (right before the destructor).
        delete Context::GP;
        delete Context::Beam;
        delete Context::RfP;
        delete Context::Slice;
    }
};



TEST_F(testMusic, constructor1)
{
    double epsilon = 1e-8;
    auto slices = Context::Slice;
    auto beam = Context::Beam;
    auto params = std::string(TEST_FILES "/Impedances/") +
                  "Music/constructor1/";

    double R_S = 1e7;
    double frequency_R = 1e8;
    double Q = 1;
    auto music = Music(beam, {R_S, 2 * constant::pi * frequency_R, Q}, N_p, N_b);

    f_vector_t v;
    util::read_vector_from_file(v, params + "alpha.txt");
    ASSERT_EQ(v.size(), 1);
    EXPECT_NEAR(v[0], music.alpha, epsilon * abs(v[0]))
            << "Testing of alpha failed\n";

    util::read_vector_from_file(v, params + "omega_bar.txt");
    ASSERT_EQ(v.size(), 1);
    EXPECT_NEAR(v[0], music.omega_bar, epsilon * abs(v[0]))
            << "Testing of omega_bar failed\n";

    util::read_vector_from_file(v, params + "constant.txt");
    ASSERT_EQ(v.size(), 1);
    EXPECT_NEAR(v[0], music.constant, epsilon * abs(v[0]))
            << "Testing of constant failed\n";

    util::read_vector_from_file(v, params + "coeff1.txt");
    ASSERT_EQ(v.size(), 1);
    EXPECT_NEAR(v[0], music.coeff1, epsilon * abs(v[0]))
            << "Testing of coeff1 failed\n";

    util::read_vector_from_file(v, params + "coeff2.txt");
    ASSERT_EQ(v.size(), 1);
    EXPECT_NEAR(v[0], music.coeff2, epsilon * abs(v[0]))
            << "Testing of coeff2 failed\n";

    util::read_vector_from_file(v, params + "coeff3.txt");
    ASSERT_EQ(v.size(), 1);
    EXPECT_NEAR(v[0], music.coeff3, epsilon * abs(v[0]))
            << "Testing of coeff3 failed\n";

    util::read_vector_from_file(v, params + "coeff4.txt");
    ASSERT_EQ(v.size(), 1);
    EXPECT_NEAR(v[0], music.coeff4, epsilon * abs(v[0]))
            << "Testing of coeff4 failed\n";

}


TEST_F(testMusic, track1)
{
    double epsilon = 1e-8;
    auto slices = Context::Slice;
    auto beam = Context::Beam;
    f_vector_t v;
    string params = TEST_FILES "/Impedances/Music/track1/";

    slices->track();

    double R_S = 1e7;
    double frequency_R = 1e8;
    double Q = 1;
    auto music = Music(beam, {R_S, 2 * constant::pi * frequency_R, Q}, N_p, N_b);
    music.track();

    util::read_vector_from_file(v, params + "induced_voltage.txt");
    ASSERT_NEAR_LOOP(v, music.induced_voltage, "induced_voltage", epsilon);

    util::read_vector_from_file(v, params + "beam_dE.txt");
    ASSERT_NEAR_LOOP(v, music.Beam->dE, "beam_dE", epsilon);

}


int main(int ac, char *av[])
{
    ::testing::InitGoogleTest(&ac, av);
    return RUN_ALL_TESTS();
}

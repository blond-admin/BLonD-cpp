
#include <blond/globals.h>
#include <blond/input_parameters/GeneralParameters.h>
#include <blond/input_parameters/RfParameters.h>
#include <blond/impedances/Music.h>
#include <blond/utilities.h>
#include <testing_utilities.h>
#include <gtest/gtest.h>

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

        f_vector_2d_t alphaVec(n_rf_systems, f_vector_t(alpha_order + 1, alpha));

        f_vector_t CVec(n_rf_systems, C);

        f_vector_2d_t hVec(n_rf_systems, f_vector_t(N_t + 1, h));

        f_vector_2d_t voltageVec(n_rf_systems, f_vector_t(N_t + 1, V));

        f_vector_2d_t dphiVec(n_rf_systems, f_vector_t(N_t + 1, dphi));

        Context::GP = new GeneralParameters(N_t, CVec, alphaVec,
                                            alpha_order, momentumVec,
                                            GeneralParameters::particle_t::proton);
        auto GP = Context::GP;

        auto Beam = Context::Beam = new Beams(GP, N_p, N_b);

        auto RfP = Context::RfP = new RfParameters(GP, n_rf_systems, hVec,
                voltageVec, dphiVec);



        // longitudinal_bigaussian(GP, RfP, Beam, tau_0 / 4, 0, -1, false);

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

    // for (uint i = 0; i < v.size(); ++i) {
    //     double ref = v[i];
    //     double real = indVoltTime->fTimeArray[i];
    //     ASSERT_NEAR(ref, real, epsilon * std::max(std::abs(ref), std::abs(real)))
    //             << "Testing of indVoltTime->fTimeArray failed on i " << i
    //             << std::endl;
    // }

    // util::read_vector_from_file(v, params + "total_wake.txt");
    // ASSERT_EQ(v.size(), indVoltTime->fTotalWake.size());

    // for (uint i = 0; i < v.size(); ++i) {
    //     double ref = v[i];
    //     double real = indVoltTime->fTotalWake[i];
    //     ASSERT_NEAR(ref, real, epsilon * std::max(std::abs(ref), std::abs(real)))
    //             << "Testing of indVoltTime->fTotalWake failed on i " << i
    //             << std::endl;
    // }

    // util::read_vector_from_file(v, params + "cut.txt");
    // ASSERT_EQ(v.size(), 1);
    // ASSERT_EQ(v[0], indVoltTime->fCut) << "Testing of fCut failed\n";

    // util::read_vector_from_file(v, params + "fshape.txt");
    // ASSERT_EQ(v.size(), 1);
    // ASSERT_EQ(v[0], indVoltTime->fShape) << "Testing of fShape failed\n";

    // delete indVoltTime;
}




int main(int ac, char *av[])
{
    ::testing::InitGoogleTest(&ac, av);
    return RUN_ALL_TESTS();
}

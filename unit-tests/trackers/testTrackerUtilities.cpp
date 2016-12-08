#include <blond/beams/Distributions.h>
#include <blond/input_parameters/GeneralParameters.h>
#include <blond/llrf/PhaseLoop.h>
#include <blond/math_functions.h>
#include <blond/trackers/Tracker.h>
#include <blond/utilities.h>
#include <gtest/gtest.h>
#include <blond/trackers/utilities.h>
#include <testing_utilities.h>



using namespace std;

class testTrackerUtilities : public ::testing::Test {

protected:
    // Machine and RF parameters
    const double C = 26658.883;                  // Machine circumference [m]
    const int h = 35640;                        // Harmonic number
    const double dphi = 0.;                      // Phase modulation/offset
    const double gamma_t = 55.759505;            // Transition gamma
    const double alpha = 1. / gamma_t / gamma_t; // First order mom. comp. factor
    const double bl_target = 1.25e-9;            // 4 sigma r.m.s. target bunch length in [s]
    const double p_s = 450.e9;
    const double V = 6e6;
    const double tau_0 = 0.4e-9;
    // Tracking details
    const long long N_b = 1e9;                  // Intensity
    const int alpha_order = 1;
    const int n_sections = 1;
    const int N_t = 5000;                      // Number of turns to track; full ramp: 8700001
    const int N_p = 10000;                      // Macro-particles
    const int N_slices = 200;

    virtual void SetUp()
    {
        omp_set_num_threads(1);


        f_vector_2d_t momentumVec(n_sections, f_vector_t(N_t + 1));
        for (auto &v : momentumVec)
            mymath::linspace(v.data(), p_s, 1.01 * p_s, N_t + 1);

        f_vector_2d_t alphaVec(n_sections, f_vector_t(alpha_order, alpha));

        f_vector_t CVec(n_sections, C);

        f_vector_2d_t hVec(n_sections, f_vector_t(N_t + 1, h));

        f_vector_2d_t voltageVec(n_sections, f_vector_t(N_t + 1, V));

        f_vector_2d_t dphiVec(n_sections, f_vector_t(N_t + 1, dphi));

        Context::GP = new GeneralParameters(N_t, CVec, alphaVec, momentumVec,
                                            GeneralParameters::particle_t::proton);

        // Context::Beam = new Beams(N_p, N_b);

        auto GP = Context::GP;
        auto Beam = Context::Beam = new Beams(GP, N_p, N_b);

        auto RfP = Context::RfP = new RfParameters(GP, n_sections, hVec,
                voltageVec, dphiVec);


        longitudinal_bigaussian(GP, RfP, Beam, tau_0 / 4, 0, -1, false);
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

TEST(phaseModuloBelow, test1)
{
    auto epsilon = 1e-8;
    auto real = phase_modulo_below_transition(constant::pi);
    ASSERT_NEAR(real, -3.141592653, epsilon);

    real = phase_modulo_below_transition(constant::pi / 2);
    ASSERT_NEAR(real, 1.570796326, epsilon);

    real = phase_modulo_below_transition(constant::pi / 3);
    ASSERT_NEAR(real, 1.0471975511, epsilon);

    real = phase_modulo_below_transition(constant::pi / 4);
    ASSERT_NEAR(real, 0.78539816339, epsilon);

    real = phase_modulo_below_transition(3 * constant::pi / 2);
    ASSERT_NEAR(real, -1.570796326, epsilon);

}

TEST(phaseModuloBelow, test2)
{
    auto epsilon = 1e-8;
    auto real = phase_modulo_below_transition(constant::pi / 13);
    ASSERT_NEAR(real, 0.241660973353, epsilon);

    real = phase_modulo_below_transition(11 * constant::pi / 72);
    ASSERT_NEAR(real, 0.479965544298, epsilon);

    real = phase_modulo_below_transition(39 * constant::pi / 6);
    ASSERT_NEAR(real, 1.570796326794 , epsilon);

    real = phase_modulo_below_transition(constant::pi * 1.324);
    ASSERT_NEAR(real, -2.12371663382 , epsilon);
}

TEST(phaseModuloBelow, test3)
{
    auto epsilon = 1e-8;
    auto real = phase_modulo_below_transition(-constant::pi);
    ASSERT_NEAR(real, -3.1415926535897, epsilon);

    real = phase_modulo_below_transition(- constant::pi / 6);
    ASSERT_NEAR(real, -0.523598775598, epsilon);

    real = phase_modulo_below_transition(- 5 * constant::pi / 4);
    ASSERT_NEAR(real, 2.3561944901923 , epsilon);

    real = phase_modulo_below_transition(- 2 * constant::pi);
    ASSERT_NEAR(real, 0.0, epsilon);
}


TEST(phaseModuloAbove, test1)
{
    auto epsilon = 1e-8;
    auto real = phase_modulo_above_transition(constant::pi);
    ASSERT_NEAR(real, 3.141592653, epsilon);

    real = phase_modulo_above_transition(constant::pi / 2);
    ASSERT_NEAR(real, 1.570796326, epsilon);

    real = phase_modulo_above_transition(constant::pi / 3);
    ASSERT_NEAR(real, 1.0471975511, epsilon);

    real = phase_modulo_above_transition(constant::pi / 4);
    ASSERT_NEAR(real, 0.78539816339, epsilon);

    real = phase_modulo_above_transition(3 * constant::pi / 2);
    ASSERT_NEAR(real, 4.7123889803846, epsilon);

}

TEST(phaseModuloAbove, test2)
{
    auto epsilon = 1e-8;
    auto real = phase_modulo_above_transition(constant::pi / 13);
    ASSERT_NEAR(real, 0.241660973353, epsilon);

    real = phase_modulo_above_transition(11 * constant::pi / 72);
    ASSERT_NEAR(real, 0.479965544298, epsilon);

    real = phase_modulo_above_transition(39 * constant::pi / 6);
    ASSERT_NEAR(real, 1.570796326794 , epsilon);

    real = phase_modulo_above_transition(constant::pi * 1.324);
    ASSERT_NEAR(real, 4.1594686733528 , epsilon);
}

TEST(phaseModuloAbove, test3)
{
    auto epsilon = 1e-8;
    auto real = phase_modulo_above_transition(-constant::pi);
    ASSERT_NEAR(real, 3.1415926535897, epsilon);

    real = phase_modulo_above_transition(- constant::pi / 6);
    ASSERT_NEAR(real, 5.75958653158128, epsilon);

    real = phase_modulo_above_transition(- 5 * constant::pi / 4);
    ASSERT_NEAR(real, 2.35619449019234 , epsilon);

    real = phase_modulo_above_transition(- 2 * constant::pi);
    ASSERT_NEAR(real, 0.0, epsilon);
}


TEST_F(testTrackerUtilities, hamiltonian1)
{
    auto epsilon = 1e-8;
    auto GP = Context::GP;
    auto Beam = Context::Beam;
    auto RfP = Context::RfP;
    auto real = hamiltonian(GP, RfP, Beam, 0.0, 0.0);
    ASSERT_NEAR(real, 2915556.62852, 2915556.62852 * epsilon);
}

TEST_F(testTrackerUtilities, hamiltonian2)
{
    auto epsilon = 1e-8;
    auto GP = Context::GP;
    auto Beam = Context::Beam;
    auto RfP = Context::RfP;
    auto real = hamiltonian(GP, RfP, Beam, 1e-7, 0.0);
    ASSERT_NEAR(real, 2828529.22602, 2828529.22602 * epsilon);
}

TEST_F(testTrackerUtilities, hamiltonian3)
{
    auto epsilon = 1e-8;
    auto GP = Context::GP;
    auto Beam = Context::Beam;
    auto RfP = Context::RfP;
    auto real = hamiltonian(GP, RfP, Beam, 1e-7, 1e7);
    ASSERT_NEAR(real, 2831020.19352, 2831020.19352* epsilon);
}

TEST_F(testTrackerUtilities, hamiltonian4)
{
    auto epsilon = 1e-8;
    auto GP = Context::GP;
    auto Beam = Context::Beam;
    auto RfP = Context::RfP;
    auto real = hamiltonian(GP, RfP, Beam, Beam->dt[100], Beam->dE[200]);
    ASSERT_NEAR(real, 4767.48080202, 4767.48080202 * epsilon);
}


TEST_F(testTrackerUtilities, hamiltonian5)
{
    auto params = std::string(TEST_FILES) +
                  "/TrackerUtilities/hamiltonian5/";
    auto epsilon = 1e-8;
    auto GP = Context::GP;
    auto Beam = Context::Beam;
    auto RfP = Context::RfP;
    auto dt = f_vector_t(Beam->dt.begin() + 100, Beam->dt.begin() + 200);
    auto dE = f_vector_t(Beam->dE.begin() + 100, Beam->dE.begin() + 200);

    auto ham = hamiltonian(GP, RfP, Beam, dt.data(), dE.data(), dt.size());

    f_vector_t v;

    util::read_vector_from_file(v, params + "hamiltonian.txt");
    ASSERT_EQ(v.size(), ham.size());
    for (uint i = 0; i < v.size(); ++i) {
        auto ref = v[i];
        auto real = ham[i];
        ASSERT_NEAR(ref, real, epsilon * min(fabs(ref), fabs(real)))
                << "Testing of hamiltonian failed on i " << i << "\n";
    }

}


TEST_F(testTrackerUtilities, hamiltonian6)
{
    auto params = std::string(TEST_FILES) +
                  "/TrackerUtilities/hamiltonian6/";
    auto epsilon = 1e-6;
    auto GP = Context::GP;
    auto Beam = Context::Beam;
    auto RfP = Context::RfP;

    auto ham = hamiltonian(GP, RfP, Beam, Beam->dt.data(),
                           Beam->dE.data(), Beam->dt.size());

    f_vector_t v;

    util::read_vector_from_file(v, params + "hamiltonian.txt");
    ASSERT_EQ(v.size(), ham.size());
    for (uint i = 0; i < v.size(); ++i) {
        auto ref = v[i];
        auto real = ham[i];
        ASSERT_NEAR(ref, real, epsilon * min(fabs(ref), fabs(real)))
                << "Testing of hamiltonian failed on i " << i << "\n";
    }

}


TEST_F(testTrackerUtilities, is_in_separatrix2)
{
    auto params = std::string(TEST_FILES) +
                  "/TrackerUtilities/is_in_separatrix2/";
    auto GP = Context::GP;
    auto Beam = Context::Beam;
    auto RfP = Context::RfP;

    int i = 0;
    for (auto &t : Beam->dt)
        t += (1.0 * (i++) / N_p) * t;

    i = 0;
    for (auto &t : Beam->dE)
        t += (1.0 * (i++) / N_p) * t;


    auto dt = f_vector_t(Beam->dt.begin() + 5000, Beam->dt.begin() + 6000);
    auto dE = f_vector_t(Beam->dE.begin() + 5000, Beam->dE.begin() + 6000);


    auto sep = is_in_separatrix(GP, RfP, Beam, dt, dE);

    f_vector_t v;

    util::read_vector_from_file(v, params + "separatrix.txt");
    ASSERT_EQ(v.size(), sep.size());
    for (uint i = 0; i < v.size(); ++i) {
        bool ref = v[i];
        bool real = sep[i];
        ASSERT_EQ(ref, real)
                << "Testing of separatrix failed on i " << i << "\n";
    }

}


TEST_F(testTrackerUtilities, is_in_separatrix1)
{
    auto params = std::string(TEST_FILES) +
                  "/TrackerUtilities/is_in_separatrix1/";
    auto GP = Context::GP;
    auto Beam = Context::Beam;
    auto RfP = Context::RfP;

    int i = 0;
    for (auto &t : Beam->dt)
        t += (1.0 * (i++) / N_p) * t;

    auto sep = is_in_separatrix(GP, RfP, Beam, Beam->dt, Beam->dE);

    f_vector_t v;

    util::read_vector_from_file(v, params + "separatrix.txt");
    ASSERT_EQ(v.size(), sep.size());
    for (uint i = 0; i < v.size(); ++i) {
        bool ref = v[i];
        bool real = sep[i];
        ASSERT_EQ(ref, real)
                << "Testing of separatrix failed on i " << i << "\n";
    }

}



TEST(testMinMax, minmax_location1)
{
    auto epsilon = 1e-8;
    auto params = string(TEST_FILES "/Distributions/minmax_location1/");
    int size = 100;
    f_vector_t x = mymath::linspace(0.0, size, size);
    f_vector_t f(size);

    for (int i = 0; i < size; i++) f[i] = cos(i);
    f_vector_t min_x_position, max_x_position;
    f_vector_t min_values, max_values;

    minmax_location(x, f, min_x_position, max_x_position,
                    min_values, max_values);
    f_vector_t v;

    util::read_vector_from_file(v, params + "min_x_position.txt");
    ASSERT_NEAR_LOOP(v, min_x_position, "min_x_position", epsilon);

    util::read_vector_from_file(v, params + "max_x_position.txt");
    ASSERT_NEAR_LOOP(v, max_x_position, "max_x_position", epsilon);

    util::read_vector_from_file(v, params + "min_values.txt");
    ASSERT_NEAR_LOOP(v, min_values, "min_values", epsilon);

    util::read_vector_from_file(v, params + "max_values.txt");
    ASSERT_NEAR_LOOP(v, max_values, "max_values", epsilon);

}


TEST(testMinMax, minmax_location2)
{
    auto epsilon = 1e-8;
    auto params = string(TEST_FILES "/Distributions/minmax_location2/");
    int size = 200;
    f_vector_t x = mymath::linspace(0.0, size, size);
    f_vector_t f(size);

    for (int i = 0; i < size; i++) f[i] = (i % 2) ? cos(i) : sin(i);

    f_vector_t min_x_position, max_x_position;
    f_vector_t min_values, max_values;

    minmax_location(x, f, min_x_position, max_x_position,
                    min_values, max_values);
    f_vector_t v;

    util::read_vector_from_file(v, params + "min_x_position.txt");
    ASSERT_NEAR_LOOP(v, min_x_position, "min_x_position", epsilon);

    util::read_vector_from_file(v, params + "max_x_position.txt");
    ASSERT_NEAR_LOOP(v, max_x_position, "max_x_position", epsilon);

    util::read_vector_from_file(v, params + "min_values.txt");
    ASSERT_NEAR_LOOP(v, min_values, "min_values", epsilon);

    util::read_vector_from_file(v, params + "max_values.txt");
    ASSERT_NEAR_LOOP(v, max_values, "max_values", epsilon);

}



TEST(testPotentialWell, potential_well_cut1)
{
    auto epsilon = 1e-8;
    auto params = string(TEST_FILES "/Distributions/potential_well_cut1/");
    int size = 100;
    f_vector_t theta_coord_array = mymath::linspace(0.0, size, size);
    f_vector_t potential_array(size);

    for (int i = 0; i < size; i++)
        potential_array[i] = (i % 2) ? cos(i) : sin(i);

    f_vector_t theta_coord_sep, potential_well_sep;
    f_vector_t min_values, max_values;

    potential_well_cut(theta_coord_array, potential_array,
                       theta_coord_sep, potential_well_sep);

    f_vector_t v;

    util::read_vector_from_file(v, params + "theta_coord_sep.txt");
    ASSERT_NEAR_LOOP(v, theta_coord_sep, "theta_coord_sep", epsilon);

    util::read_vector_from_file(v, params + "potential_well_sep.txt");
    ASSERT_NEAR_LOOP(v, potential_well_sep, "potential_well_sep", epsilon);

}


TEST(testPotentialWell, potential_well_cut_deathtest1)
{
    int size = 5;
    f_vector_t theta_coord_array = mymath::linspace(0.0, size, size);
    f_vector_t potential_array(size);

    for (int i = 0; i < size; i++)
        potential_array[i] = (i % 2) ? cos(i) : sin(i);

    f_vector_t theta_coord_sep, potential_well_sep;
    f_vector_t min_values, max_values;

    ASSERT_DEATH(potential_well_cut(theta_coord_array, potential_array,
                                    theta_coord_sep, potential_well_sep)
                 , "[potential_well_cut]\\s+");

}


TEST(testPotentialWell, potential_well_cut_deathtest2)
{
    int size = 10;
    f_vector_t theta_coord_array = mymath::linspace(0.0, size, size);
    f_vector_t potential_array(size);

    for (int i = 0; i < size; i++)
        potential_array[i] = cos(i);

    f_vector_t theta_coord_sep, potential_well_sep;
    f_vector_t min_values, max_values;

    ASSERT_DEATH(potential_well_cut(theta_coord_array, potential_array,
                                    theta_coord_sep, potential_well_sep)
                 , "[potential_well_cut]\\s+");

}


TEST(testPotentialWell, potential_well_cut2)
{
    auto epsilon = 1e-8;
    auto params = string(TEST_FILES "/Distributions/potential_well_cut2/");
    int size = 10;
    f_vector_t theta_coord_array = mymath::linspace(0.0, size, size);
    f_vector_t potential_array(size);

    for (int i = 0; i < size; i++)
        potential_array[i] = sin(i);

    f_vector_t theta_coord_sep, potential_well_sep;
    f_vector_t min_values, max_values;

    potential_well_cut(theta_coord_array, potential_array,
                       theta_coord_sep, potential_well_sep);

    f_vector_t v;

    util::read_vector_from_file(v, params + "theta_coord_sep.txt");
    ASSERT_NEAR_LOOP(v, theta_coord_sep, "theta_coord_sep", epsilon);

    util::read_vector_from_file(v, params + "potential_well_sep.txt");
    ASSERT_NEAR_LOOP(v, potential_well_sep, "potential_well_sep", epsilon);

}



TEST(testPotentialWell, potential_well_cut_deathtest3)
{
    auto epsilon = 1e-8;
    int size = 10;
    f_vector_t theta_coord_array = mymath::linspace(0.0, size, size);
    f_vector_t potential_array(size);

    for (int i = 0; i < size; i++)
        potential_array[i] = exp(i);

    f_vector_t theta_coord_sep, potential_well_sep;
    f_vector_t min_values, max_values;

    ASSERT_DEATH(potential_well_cut(theta_coord_array, potential_array,
                                    theta_coord_sep, potential_well_sep)
                 , "[potential_well_cut]\\s+");
}






int main(int ac, char *av[])
{
    ::testing::InitGoogleTest(&ac, av);
    return RUN_ALL_TESTS();
}

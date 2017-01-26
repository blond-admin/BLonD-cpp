#include <blond/blond.h>
#include <testing_utilities.h>
#include <gtest/gtest.h>

using namespace std;

class testGP : public ::testing::Test {

protected:
    virtual void SetUp()
    {
        omp_set_num_threads(1);

        f_vector_2d_t momentumVec(n_sections, f_vector_t(N_t + 1));
        for (auto &v : momentumVec)
            mymath::linspace(v.data(), p_i, p_f, N_t + 1);

        f_vector_2d_t alphaVec(n_sections, f_vector_t(alpha_order, alpha));

        f_vector_t CVec(n_sections, C);

        Context::GP = new GeneralParameters(N_t, CVec, alphaVec, momentumVec,
                                            GeneralParameters::particle_t::proton);
    }

    virtual void TearDown()
    {
        // Code here will be called immediately after each test
        // (right before the destructor).
        delete Context::GP;
    }

private:
    // Machine and RF parameters
    const double C = 26658.883;       // Machine circumference [m]
    const double p_i = 450e9;     // Synchronous momentum [eV/c]
    const double p_f = 460.005e9;     // Synchronous momentum, final
    const double gamma_t = 55.759505; // Transition gamma
    const double alpha =
        1.0 / gamma_t / gamma_t; // First order mom. comp. factor
    const int alpha_order = 1;
    const int n_sections = 1;
    // Tracking details

    const int N_t = 2000; // Number of turns to track


};

TEST_F(testGP, charge)
{
    const double epsilon = 1e-8;

    string GP_params = TEST_FILES "/GP/GP_params/";
    vector<double> v;
    util::read_vector_from_file(v, GP_params + "charge");
    double ref = v[0];
    double real = Context::GP->charge;
    ASSERT_NEAR(ref, real, epsilon * max(abs(ref), abs(real)));
}

TEST_F(testGP, mass)
{
    const double epsilon = 1e-8;
    string GP_params = TEST_FILES "/GP/GP_params/";
    vector<double> v;
    util::read_vector_from_file(v, GP_params + "mass");
    double ref = v[0];
    double real = Context::GP->mass;
    ASSERT_NEAR(ref, real, epsilon * max(abs(ref), abs(real)));
}

TEST_F(testGP, ring_radius)
{
    const double epsilon = 1e-8;
    string GP_params = TEST_FILES "/GP/GP_params/";
    vector<double> v;
    util::read_vector_from_file(v, GP_params + "ring_radius");
    double ref = v[0];
    double real = Context::GP->ring_radius;
    ASSERT_NEAR(ref, real, epsilon * max(abs(ref), abs(real)));
}

TEST_F(testGP, t_rev)
{
    const double epsilon = 1e-8;
    string GP_params = TEST_FILES "/GP/GP_params/";
    vector<double> v;
    util::read_vector_from_file(v, GP_params + "t_rev.txt");
    ASSERT_NEAR_LOOP(v, Context::GP->t_rev, "t_rev", epsilon);
}

TEST_F(testGP, cycle_time)
{
    const double epsilon = 1e-8;
    string GP_params = TEST_FILES "/GP/GP_params/";
    vector<double> v;
    util::read_vector_from_file(v, GP_params + "cycle_time");
    ASSERT_NEAR_LOOP(v, Context::GP->cycle_time, "cycle_time", epsilon);
}

TEST_F(testGP, omega_rev)
{
    const double epsilon = 1e-8;
    string GP_params = TEST_FILES "/GP/GP_params/";
    vector<double> v;
    util::read_vector_from_file(v, GP_params + "omega_rev.txt");
    ASSERT_NEAR_LOOP(v, Context::GP->omega_rev, "omega_rev", epsilon);
}

TEST_F(testGP, eta_0)
{
    const double epsilon = 1e-8;
    string GP_params = TEST_FILES "/GP/GP_params/";
    vector<double> v;
    util::read_vector_from_file(v, GP_params + "eta_0[0].txt");
    ASSERT_NEAR_LOOP(v, Context::GP->eta_0[0], "eta_0[0]", epsilon);
}


TEST(higher_alpha_order, eta_1_1)
{
    const double C = 26658.883;       // Machine circumference [m]
    const double p_i = 100e8;     // Synchronous momentum [eV/c]
    const double p_f = 110e8;     // Synchronous momentum, final
    const double gamma_t = 55.759505; // Transition gamma
    const double alpha =
        1.0 / gamma_t / gamma_t; // First order mom. comp. factor
    const int alpha_order = 2;
    const int n_sections = 1;
    // Tracking details

    const int N_t = 100; // Number of turns to track

    omp_set_num_threads(1);

    f_vector_2d_t momentumVec(n_sections);
    for (auto &v : momentumVec)
        v = mymath::linspace(p_i, p_f, N_t + 1);

    f_vector_2d_t alphaVec(n_sections, {alpha, alpha / gamma_t});

    f_vector_t CVec(n_sections, C);

    auto GP = GeneralParameters(N_t, CVec, alphaVec,
                                momentumVec,
                                GeneralParameters::particle_t::electron);

    f_vector_t v;
    string params = TEST_FILES "/GP/higher_alpha_order/eta_1_1/";

    util::read_vector_from_file(v, params + "eta_0[0].txt");
    ASSERT_NEAR_LOOP(v, GP.eta_0[0], "eta_0[0]");

    util::read_vector_from_file(v, params + "eta_1[0].txt");
    ASSERT_NEAR_LOOP(v, GP.eta_1[0], "eta_1[0]");

}

TEST(higher_alpha_order, eta_1_2)
{
    const double C = 26658.883;       // Machine circumference [m]
    const double p_i = 50e5;     // Synchronous momentum [eV/c]
    const double p_f = 80e5;     // Synchronous momentum, final
    const double gamma_t = 15.759505; // Transition gamma
    const double alpha =
        1.0 / gamma_t / gamma_t; // First order mom. comp. factor
    const int alpha_order = 2;
    const int n_sections = 1;
    // Tracking details

    const int N_t = 50; // Number of turns to track

    omp_set_num_threads(1);

    f_vector_2d_t momentumVec(n_sections);
    for (auto &v : momentumVec)
        v = mymath::linspace(p_i, p_f, N_t + 1);

    f_vector_2d_t alphaVec(n_sections, {alpha, alpha / gamma_t});

    f_vector_t CVec(n_sections, C);

    auto GP = GeneralParameters(N_t, CVec, alphaVec,
                                momentumVec,
                                GeneralParameters::particle_t::electron);

    f_vector_t v;
    string params = TEST_FILES "/GP/higher_alpha_order/eta_1_2/";

    util::read_vector_from_file(v, params + "eta_0[0].txt");
    ASSERT_NEAR_LOOP(v, GP.eta_0[0], "eta_0[0]");

    util::read_vector_from_file(v, params + "eta_1[0].txt");
    ASSERT_NEAR_LOOP(v, GP.eta_1[0], "eta_1[0]");

}


TEST(higher_alpha_order, eta_2_1)
{
    const double C = 26658.883;       // Machine circumference [m]
    const double p_i = 50e5;     // Synchronous momentum [eV/c]
    const double p_f = 80e5;     // Synchronous momentum, final
    const double gamma_t = 15.759505; // Transition gamma
    const double alpha =
        1.0 / gamma_t / gamma_t; // First order mom. comp. factor
    const int alpha_order = 3;
    const int n_sections = 1;
    // Tracking details

    const int N_t = 50; // Number of turns to track

    omp_set_num_threads(1);

    f_vector_2d_t momentumVec(n_sections);
    for (auto &v : momentumVec)
        v = mymath::linspace(p_i, p_f, N_t + 1);

    f_vector_2d_t alphaVec(n_sections, {alpha, alpha / gamma_t, 2 * alpha / gamma_t});

    f_vector_t CVec(n_sections, C);

    auto GP = GeneralParameters(N_t, CVec, alphaVec, momentumVec,
                                GeneralParameters::particle_t::electron);

    f_vector_t v;
    string params = TEST_FILES "/GP/higher_alpha_order/eta_2_1/";

    util::read_vector_from_file(v, params + "eta_0[0].txt");
    ASSERT_NEAR_LOOP(v, GP.eta_0[0], "eta_0[0]");

    util::read_vector_from_file(v, params + "eta_1[0].txt");
    ASSERT_NEAR_LOOP(v, GP.eta_1[0], "eta_1[0]");

    util::read_vector_from_file(v, params + "eta_2[0].txt");
    ASSERT_NEAR_LOOP(v, GP.eta_2[0], "eta_2[0]");

}


TEST(particle_types, electron_electron1)
{
    const double C = 26658.883;       // Machine circumference [m]
    const double p_i = 50e5;     // Synchronous momentum [eV/c]
    const double p_f = 80e5;     // Synchronous momentum, final
    const double gamma_t = 15.759505; // Transition gamma
    const double alpha =
        1.0 / gamma_t / gamma_t; // First order mom. comp. factor
    const int alpha_order = 1;
    const int n_sections = 1;
    // Tracking details

    const int N_t = 50; // Number of turns to track

    omp_set_num_threads(1);

    f_vector_2d_t momentumVec({mymath::linspace(p_i, p_f, N_t + 1)});

    f_vector_2d_t alphaVec({{alpha}});

    f_vector_t CVec{C};

    auto GP = GeneralParameters(N_t, CVec, alphaVec,
                                momentumVec,
                                GeneralParameters::particle_t::electron,
                                0, 0, GeneralParameters::electron, 0, 0,
                                n_sections);

    f_vector_t v;
    string params = TEST_FILES "/GP/particle_types/electron_electron1/";

    util::read_vector_from_file(v, params + "mass.txt");
    ASSERT_NEAR_LOOP(v, {GP.mass}, "mass");
    util::read_vector_from_file(v, params + "charge.txt");
    ASSERT_NEAR_LOOP(v, {GP.charge}, "charge");

    util::read_vector_from_file(v, params + "mass2.txt");
    ASSERT_NEAR_LOOP(v, {GP.mass2}, "mass2");
    util::read_vector_from_file(v, params + "charge2.txt");
    ASSERT_NEAR_LOOP(v, {GP.charge2}, "charge2");

}


TEST(particle_types, input_input1)
{
    const double C = 26658.883;       // Machine circumference [m]
    const double p_i = 50e5;     // Synchronous momentum [eV/c]
    const double p_f = 80e5;     // Synchronous momentum, final
    const double gamma_t = 15.759505; // Transition gamma
    const double alpha =
        1.0 / gamma_t / gamma_t; // First order mom. comp. factor
    const int alpha_order = 1;
    const int n_sections = 1;
    // Tracking details

    const int N_t = 50; // Number of turns to track

    omp_set_num_threads(1);

    f_vector_2d_t momentumVec({mymath::linspace(p_i, p_f, N_t + 1)});

    f_vector_2d_t alphaVec({{alpha}});

    f_vector_t CVec{C};

    auto GP = GeneralParameters(N_t, CVec, alphaVec,
                                momentumVec,
                                GeneralParameters::user_input, 1.0, 2.0,
                                GeneralParameters::user_input, 3.0, 4.0,
                                n_sections);

    auto epsilon = 1e-8;
    ASSERT_NEAR(GP.mass, 1.0,  epsilon * min(1.0, GP.mass));
    ASSERT_NEAR(GP.charge, 2.0,  epsilon * min(1.0, GP.charge));

    ASSERT_NEAR(GP.mass2, 3.0,  epsilon * min(1.0, GP.mass2));
    ASSERT_NEAR(GP.charge2, 4.0,  epsilon * min(1.0, GP.charge2));
}

TEST(particle_types, death_test1)
{
    const double C = 26658.883;       // Machine circumference [m]
    const double p_i = 50e5;     // Synchronous momentum [eV/c]
    const double p_f = 80e5;     // Synchronous momentum, final
    const double gamma_t = 15.759505; // Transition gamma
    const double alpha =
        1.0 / gamma_t / gamma_t; // First order mom. comp. factor
    const int alpha_order = 1;
    const int n_sections = 1;
    // Tracking details

    const int N_t = 50; // Number of turns to track

    omp_set_num_threads(1);

    f_vector_2d_t momentumVec({mymath::linspace(p_i, p_f, N_t + 1)});

    f_vector_2d_t alphaVec({{alpha}});

    f_vector_t CVec{C};

    auto GP = GeneralParameters(N_t, CVec, alphaVec,
                                momentumVec,
                                GeneralParameters::user_input, 1.0, 2.0,
                                GeneralParameters::user_input, 3.0, 4.0,
                                n_sections);

    auto epsilon = 1e-8;
    ASSERT_NEAR(GP.mass, 1.0,  epsilon * min(1.0, GP.mass));
    ASSERT_NEAR(GP.charge, 2.0,  epsilon * min(1.0, GP.charge));

    ASSERT_NEAR(GP.mass2, 3.0,  epsilon * min(1.0, GP.mass2));
    ASSERT_NEAR(GP.charge2, 4.0,  epsilon * min(1.0, GP.charge2));
}

int main(int ac, char *av[])
{
    ::testing::InitGoogleTest(&ac, av);
    return RUN_ALL_TESTS();
}

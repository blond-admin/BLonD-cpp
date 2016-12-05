#include <list>
#include <blond/globals.h>
#include <blond/input_parameters/GeneralParameters.h>
#include <blond/math_functions.h>
#include <blond/utilities.h>
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

        f_vector_2d_t alphaVec(n_sections, f_vector_t(alpha_order + 1, alpha));

        f_vector_t CVec(n_sections, C);

        Context::GP = new GeneralParameters(N_t, CVec, alphaVec,
                                            alpha_order, momentumVec,
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

TEST_F(testGP, test_charge)
{
    const double epsilon = 1e-8;

    string GP_params = TEST_FILES "/GP/GP_params/";
    vector<double> v;
    util::read_vector_from_file(v, GP_params + "charge");
    // cout << v[0];
    double ref = v[0];
    double real = Context::GP->charge;
    ASSERT_NEAR(ref, real, epsilon * max(std::abs(ref), std::abs(real)));
}

TEST_F(testGP, test_mass)
{
    const double epsilon = 1e-8;
    string GP_params = TEST_FILES "/GP/GP_params/";
    vector<double> v;
    util::read_vector_from_file(v, GP_params + "mass");
    // cout << v[0];
    double ref = v[0];
    double real = Context::GP->mass;
    ASSERT_NEAR(ref, real, epsilon * max(std::abs(ref), std::abs(real)));
}

TEST_F(testGP, test_ring_radius)
{
    const double epsilon = 1e-8;
    string GP_params = TEST_FILES "/GP/GP_params/";
    vector<double> v;
    util::read_vector_from_file(v, GP_params + "ring_radius");
    // cout << v[0];
    double ref = v[0];
    double real = Context::GP->ring_radius;
    ASSERT_NEAR(ref, real, epsilon * max(std::abs(ref), std::abs(real)));
}

TEST_F(testGP, test_t_rev)
{
    const double epsilon = 1e-8;
    string GP_params = TEST_FILES "/GP/GP_params/";
    vector<double> v;
    util::read_vector_from_file(v, GP_params + "t_rev");
    // cout << v.size() << endl;
    // ASSERT_EQ(v.size(), GP->n_turns+1 );
    for (unsigned int i = 0; i < v.size(); ++i) {
        double ref = v[i];
        double real = Context::GP->t_rev[i];
        ASSERT_NEAR(ref, real, epsilon * max(std::abs(ref), std::abs(real)));
    }
}

TEST_F(testGP, test_cycle_time)
{
    const double epsilon = 1e-8;
    string GP_params = TEST_FILES "/GP/GP_params/";
    vector<double> v;
    util::read_vector_from_file(v, GP_params + "cycle_time");
    // ASSERT_EQ(v.size(), GP->n_turns+1);
    for (unsigned int i = 0; i < v.size(); ++i) {
        double ref = v[i];
        double real = Context::GP->cycle_time[i];
        ASSERT_NEAR(ref, real, epsilon * max(std::abs(ref), std::abs(real)));
    }
}

TEST_F(testGP, test_omega_rev)
{
    const double epsilon = 1e-8;
    string GP_params = TEST_FILES "/GP/GP_params/";
    vector<double> v;
    util::read_vector_from_file(v, GP_params + "omega_rev");
    // cout << v.size() << endl;
    // ASSERT_EQ(v.size(), GP->n_turns+1);
    for (unsigned int i = 0; i < v.size(); ++i) {
        double ref = v[i];
        double real = Context::GP->omega_rev[i];
        ASSERT_NEAR(ref, real, epsilon * max(std::abs(ref), std::abs(real)));
    }
}

TEST_F(testGP, test_eta_0)
{
    const double epsilon = 1e-8;
    string GP_params = TEST_FILES "/GP/GP_params/";
    vector<double> v;
    util::read_vector_from_file(v, GP_params + "eta_0[0]");
    // cout << v.size() << endl;
    // ASSERT_EQ(v.size(), GP->n_turns +1);
    for (unsigned int i = 0; i < v.size(); ++i) {
        double ref = v[i];
        double real = Context::GP->eta_0[0][i];
        ASSERT_NEAR(ref, real, epsilon * max(std::abs(ref), std::abs(real)));
    }
}


TEST(testGP_higher_alpha_order, eta_1_1)
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
                                alpha_order, momentumVec,
                                GeneralParameters::particle_t::electron);

    f_vector_t v;
    string params = TEST_FILES "/GP/higher_alpha_order/eta_1_1/";

    util::read_vector_from_file(v, params + "eta_0[0].txt");
    ASSERT_NEAR_LOOP(v, GP.eta_0[0], "eta_0[0]");

    util::read_vector_from_file(v, params + "eta_1[0].txt");
    ASSERT_NEAR_LOOP(v, GP.eta_1[0], "eta_1[0]");

}

TEST(testGP_higher_alpha_order, eta_1_2)
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
                                alpha_order, momentumVec,
                                GeneralParameters::particle_t::electron);

    f_vector_t v;
    string params = TEST_FILES "/GP/higher_alpha_order/eta_1_2/";

    util::read_vector_from_file(v, params + "eta_0[0].txt");
    ASSERT_NEAR_LOOP(v, GP.eta_0[0], "eta_0[0]");

    util::read_vector_from_file(v, params + "eta_1[0].txt");
    ASSERT_NEAR_LOOP(v, GP.eta_1[0], "eta_1[0]");

}


TEST(testGP_higher_alpha_order, eta_2_1)
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

    auto GP = GeneralParameters(N_t, CVec, alphaVec,
                                alpha_order, momentumVec,
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


TEST(testGP_higher_alpha_order, eta_2_2)
{
    const double C = 26658.883;       // Machine circumference [m]
    const double p_i = 50e5;     // Synchronous momentum [eV/c]
    const double p_f = 80e5;     // Synchronous momentum, final
    const double gamma_t = 15.759505; // Transition gamma
    const double alpha =
        1.0 / gamma_t / gamma_t; // First order mom. comp. factor
    const int alpha_order = 3;
    const int n_sections = 2;
    // Tracking details

    const int N_t = 50; // Number of turns to track

    omp_set_num_threads(1);

    f_vector_2d_t momentumVec({mymath::linspace(p_i, p_f, N_t + 1),
                               mymath::linspace(1.1 * p_i, 1.2 * p_f, N_t + 1)
                              });
    // for (auto &v : momentumVec)
    //     v = mymath::linspace(p_i, p_f, N_t + 1);

    f_vector_2d_t alphaVec({{alpha, alpha / gamma_t, 2 * alpha / gamma_t},
        {1.1 * alpha, 1.1 * alpha / gamma_t, 1.1 * 2 * alpha / gamma_t}
    });

    f_vector_t CVec{C, 1.1 * C};

    auto GP = GeneralParameters(N_t, CVec, alphaVec,
                                alpha_order, momentumVec,
                                GeneralParameters::particle_t::electron,
                                0, 0, GeneralParameters::proton, 0, 0,
                                n_sections);

    f_vector_t v;
    string params = TEST_FILES "/GP/higher_alpha_order/eta_2_2/";

    util::read_vector_from_file(v, params + "eta_0[1].txt");
    ASSERT_NEAR_LOOP(v, GP.eta_0[1], "eta_0[1]");

    util::read_vector_from_file(v, params + "eta_1[1].txt");
    ASSERT_NEAR_LOOP(v, GP.eta_1[1], "eta_1[1]");

    util::read_vector_from_file(v, params + "eta_2[1].txt");
    ASSERT_NEAR_LOOP(v, GP.eta_2[1], "eta_2[1]");


}

int main(int ac, char *av[])
{
    ::testing::InitGoogleTest(&ac, av);
    return RUN_ALL_TESTS();
}

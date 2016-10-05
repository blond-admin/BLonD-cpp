#include <list>
#include <blond/globals.h>
#include <blond/input_parameters/GeneralParameters.h>
#include <blond/math_functions.h>
#include <blond/utilities.h>
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

        Context::GP = new GeneralParameters(N_t, CVec, alphaVec, alpha_order,
                                            momentumVec, proton);
    }

    virtual void TearDown()
    {
        // Code here will be called immediately after each test
        // (right before the destructor).
        delete Context::GP;
    }

private:
    // Machine and RF parameters
    const ftype C = 26658.883;       // Machine circumference [m]
    const long long p_i = 450e9;     // Synchronous momentum [eV/c]
    const ftype p_f = 460.005e9;     // Synchronous momentum, final
    const ftype gamma_t = 55.759505; // Transition gamma
    const ftype alpha =
        1.0 / gamma_t / gamma_t; // First order mom. comp. factor
    const int alpha_order = 1;
    const int n_sections = 1;
    // Tracking details

    const int N_t = 2000; // Number of turns to track


};

TEST_F(testGP, test_charge)
{
    const ftype epsilon = 1e-8;

    string GP_params = TEST_FILES "/GP/GP_params/";
    vector<ftype> v;
    util::read_vector_from_file(v, GP_params + "charge");
    // cout << v[0];
    ftype ref = v[0];
    ftype real = Context::GP->charge;
    ASSERT_NEAR(ref, real, epsilon * max(std::abs(ref), std::abs(real)));
}

TEST_F(testGP, test_mass)
{
    const ftype epsilon = 1e-8;
    string GP_params = TEST_FILES "/GP/GP_params/";
    vector<ftype> v;
    util::read_vector_from_file(v, GP_params + "mass");
    // cout << v[0];
    ftype ref = v[0];
    ftype real = Context::GP->mass;
    ASSERT_NEAR(ref, real, epsilon * max(std::abs(ref), std::abs(real)));
}

TEST_F(testGP, test_ring_radius)
{
    const ftype epsilon = 1e-8;
    string GP_params = TEST_FILES "/GP/GP_params/";
    vector<ftype> v;
    util::read_vector_from_file(v, GP_params + "ring_radius");
    // cout << v[0];
    ftype ref = v[0];
    ftype real = Context::GP->ring_radius;
    ASSERT_NEAR(ref, real, epsilon * max(std::abs(ref), std::abs(real)));
}

TEST_F(testGP, test_t_rev)
{
    const ftype epsilon = 1e-8;
    string GP_params = TEST_FILES "/GP/GP_params/";
    vector<ftype> v;
    util::read_vector_from_file(v, GP_params + "t_rev");
    // cout << v.size() << endl;
    // ASSERT_EQ(v.size(), GP->n_turns+1 );
    for (unsigned int i = 0; i < v.size(); ++i) {
        ftype ref = v[i];
        ftype real = Context::GP->t_rev[i];
        ASSERT_NEAR(ref, real, epsilon * max(std::abs(ref), std::abs(real)));
    }
}

TEST_F(testGP, test_cycle_time)
{
    const ftype epsilon = 1e-8;
    string GP_params = TEST_FILES "/GP/GP_params/";
    vector<ftype> v;
    util::read_vector_from_file(v, GP_params + "cycle_time");
    // ASSERT_EQ(v.size(), GP->n_turns+1);
    for (unsigned int i = 0; i < v.size(); ++i) {
        ftype ref = v[i];
        ftype real = Context::GP->cycle_time[i];
        ASSERT_NEAR(ref, real, epsilon * max(std::abs(ref), std::abs(real)));
    }
}

TEST_F(testGP, test_omega_rev)
{
    const ftype epsilon = 1e-8;
    string GP_params = TEST_FILES "/GP/GP_params/";
    vector<ftype> v;
    util::read_vector_from_file(v, GP_params + "omega_rev");
    // cout << v.size() << endl;
    // ASSERT_EQ(v.size(), GP->n_turns+1);
    for (unsigned int i = 0; i < v.size(); ++i) {
        ftype ref = v[i];
        ftype real = Context::GP->omega_rev[i];
        ASSERT_NEAR(ref, real, epsilon * max(std::abs(ref), std::abs(real)));
    }
}

TEST_F(testGP, test_eta_0)
{
    const ftype epsilon = 1e-8;
    string GP_params = TEST_FILES "/GP/GP_params/";
    vector<ftype> v;
    util::read_vector_from_file(v, GP_params + "eta_0[0]");
    // cout << v.size() << endl;
    // ASSERT_EQ(v.size(), GP->n_turns +1);
    for (unsigned int i = 0; i < v.size(); ++i) {
        ftype ref = v[i];
        ftype real = Context::GP->eta_0[0][i];
        ASSERT_NEAR(ref, real, epsilon * max(std::abs(ref), std::abs(real)));
    }
}

int main(int ac, char *av[])
{
    ::testing::InitGoogleTest(&ac, av);
    return RUN_ALL_TESTS();
}

#include <iostream>

#include <blond/globals.h>
#include <blond/input_parameters/GeneralParameters.h>
#include <blond/math_functions.h>
#include <blond/utilities.h>
#include <gtest/gtest.h>

using namespace std;

class testRFP : public ::testing::Test {

    const long long N_b = 1e9;  // Intensity
    // const double tau_0 = 0.4e-9; // Initial bunch length, 4 sigma [s]

    // Machine and RF parameters
    const double C = 26658.883;       // Machine circumference [m]
    const double p_i = 450e9;     // Synchronous momentum [eV/c]
    const double p_f = 460.005e9;     // Synchronous momentum, final
    const double h = 35640;       // Harmonic number
    const double V = 6e6;             // RF voltage [V]
    const double dphi = 0;            // Phase modulation/offset
    const double gamma_t = 55.759505; // Transition gamma
    const double alpha =
        1. / gamma_t / gamma_t; // First order mom. comp. factor
    const int alpha_order = 1;
    const int n_sections = 1;
    // Tracking details

    const int N_t = 2000; // Number of turns to track
    const int N_p = 100;  // Macro-particles
    // const int N_slices = 10;

    virtual void SetUp()
    {
        omp_set_num_threads(1);

        f_vector_2d_t momentumVec(n_sections, f_vector_t(N_t + 1));
        for (auto &v : momentumVec)
            mymath::linspace(v.data(), p_i, p_f, N_t + 1);

        // f_vector_2d_t alphaVec(alpha_order + 1, f_vector_t(n_sections,
        // alpha));
        f_vector_2d_t alphaVec(n_sections, f_vector_t(alpha_order, alpha));

        f_vector_t CVec(n_sections, C);

        f_vector_2d_t hVec(n_sections, f_vector_t(N_t + 1, h));

        f_vector_2d_t voltageVec(n_sections, f_vector_t(N_t + 1, V));

        f_vector_2d_t dphiVec(n_sections, f_vector_t(N_t + 1, dphi));
        // util::dump(CVec, "C\n");
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

TEST_F(testRFP, test_length_ratio)
{
    auto rfp = Context::RfP;
    double epsilon = 1e-8;
    string params = TEST_FILES "/RFP/RFP_params/";

    f_vector_t v;
    util::read_vector_from_file(v, params + "length_ratio");
    // cout << v[0];
    double ref = v[0];
    double real = rfp->length_ratio;
    ASSERT_NEAR(ref, real, epsilon * max(abs(ref), abs(real)))
            << "Testing of length_ratio has failed\n";
}

TEST_F(testRFP, test_E_increment)
{
    auto rfp = Context::RfP;
    double epsilon = 1e-8;
    string params = TEST_FILES "/RFP/RFP_params/";

    f_vector_t v;
    util::read_vector_from_file(v, params + "E_increment");
    // cout << v.size() << endl;
    for (uint i = 0; i < v.size(); ++i) {
        double ref = v[i];
        double real = rfp->E_increment[i];
        ASSERT_NEAR(ref, real, epsilon * max(abs(ref), abs(real)))
                << "Testing of E_increment has failed on i " << i << "\n";
    }
}

TEST_F(testRFP, test_phi_s)
{
    auto rfp = Context::RfP;
    double epsilon = 1e-8;
    string params = TEST_FILES "/RFP/RFP_params/";

    f_vector_t v;
    util::read_vector_from_file(v, params + "phi_s");
    // cout << v.size() << endl;
    for (uint i = 0; i < v.size(); ++i) {
        double ref = v[i];
        double real = rfp->phi_s[i];
        ASSERT_NEAR(ref, real, epsilon * max(abs(ref), abs(real)))
                << "Testing of phi_s has failed on i " << i << "\n";
    }
}

TEST_F(testRFP, test_Qs)
{
    auto rfp = Context::RfP;
    double epsilon = 1e-8;
    string params = TEST_FILES "/RFP/RFP_params/";

    f_vector_t v;
    util::read_vector_from_file(v, params + "Qs");
    // cout << v.size() << endl;
    for (uint i = 0; i < v.size(); ++i) {
        double ref = v[i];
        double real = rfp->Qs[i];
        ASSERT_NEAR(ref, real, epsilon * max(abs(ref), abs(real)))
                << "Testing of Qs has failed on i " << i << "\n";
    }
}

TEST_F(testRFP, test_omega_s0)
{
    auto rfp = Context::RfP;
    double epsilon = 1e-8;
    string params = TEST_FILES "/RFP/RFP_params/";

    f_vector_t v;
    util::read_vector_from_file(v, params + "omega_s0");
    // cout << v.size() << endl;
    for (uint i = 0; i < v.size(); ++i) {
        double ref = v[i];
        double real = rfp->omega_s0[i];
        ASSERT_NEAR(ref, real, epsilon * max(abs(ref), abs(real)))
                << "Testing of omega_s0 has failed on i " << i << "\n";
    }
}

TEST_F(testRFP, test_omega_RF_d)
{
    auto rfp = Context::RfP;
    double epsilon = 1e-8;
    string params = TEST_FILES "/RFP/RFP_params/";

    f_vector_t v;
    util::read_vector_from_file(v, params + "omega_RF_d[0]");
    // cout << v.size() << endl;
    for (uint i = 0; i < v.size(); ++i) {
        double ref = v[i];
        double real = rfp->omega_rf_d[0][i];
        ASSERT_NEAR(ref, real, epsilon * max(abs(ref), abs(real)))
                << "Testing of omega_rf_d[0] has failed on i " << i << "\n";
    }
}

TEST_F(testRFP, test_omega_RF)
{

    auto rfp = Context::RfP;
    double epsilon = 1e-8;
    string params = TEST_FILES "/RFP/RFP_params/";
    f_vector_t v;
    util::read_vector_from_file(v, params + "omega_RF[0]");
    // cout << v.size() << endl;
    for (uint i = 0; i < v.size(); ++i) {
        double ref = v[i];
        double real = rfp->omega_rf[0][i];
        ASSERT_NEAR(ref, real, epsilon * max(abs(ref), abs(real)))
                << "Testing of omega_rf[0] has failed on i " << i << "\n";
    }
}

TEST_F(testRFP, test_t_RF)
{
    auto rfp = Context::RfP;
    double epsilon = 1e-8;
    string params = TEST_FILES "/RFP/RFP_params/";

    f_vector_t v;
    util::read_vector_from_file(v, params + "t_RF");
    // cout << v.size() << endl;
    for (uint i = 0; i < v.size(); ++i) {
        double ref = v[i];
        double real = rfp->t_rf[i];
        ASSERT_NEAR(ref, real, epsilon * max(abs(ref), abs(real)))
                << "Testing of t_RF has failed on i " << i << "\n";
    }
}

int main(int ac, char *av[])
{
    ::testing::InitGoogleTest(&ac, av);
    return RUN_ALL_TESTS();
}

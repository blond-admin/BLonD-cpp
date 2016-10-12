#include <iostream>

#include <blond/globals.h>
#include <blond/input_parameters/GeneralParameters.h>
#include <blond/math_functions.h>
#include <blond/utilities.h>
#include <gtest/gtest.h>

const ftype epsilon = 1e-8;
const std::string params = TEST_FILES "/RFP/RFP_params/";

class testRFP : public ::testing::Test {

protected:
    virtual void SetUp()
    {
        omp_set_num_threads(1);

        f_vector_2d_t momentumVec(n_sections, f_vector_t(N_t + 1));
        for (auto &v : momentumVec)
            mymath::linspace(v.data(), p_i, p_f, N_t + 1);

        // f_vector_2d_t alphaVec(alpha_order + 1, f_vector_t(n_sections,
        // alpha));
        f_vector_2d_t alphaVec(n_sections, f_vector_t(alpha_order + 1, alpha));

        f_vector_t CVec(n_sections, C);

        f_vector_2d_t hVec(n_sections, f_vector_t(N_t + 1, h));

        f_vector_2d_t voltageVec(n_sections, f_vector_t(N_t + 1, V));

        f_vector_2d_t dphiVec(n_sections, f_vector_t(N_t + 1, dphi));

        Context::GP = new GeneralParameters(N_t, CVec, alphaVec,
                                            alpha_order, momentumVec,
                                            GeneralParameters::particle_t::proton);

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

private:
    const long long N_b = 1e9;  // Intensity
    // const ftype tau_0 = 0.4e-9; // Initial bunch length, 4 sigma [s]

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

    const int N_t = 2000; // Number of turns to track
    const int N_p = 100;  // Macro-particles
    // const int N_slices = 10;
};

TEST_F(testRFP, test_length_ratio)
{
    std::vector<ftype> v;
    util::read_vector_from_file(v, params + "length_ratio");
    // std::cout << v[0];
    ftype ref = v[0];
    ftype real = Context::RfP->length_ratio;
    ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)));
}

TEST_F(testRFP, test_E_increment)
{
    std::vector<ftype> v;
    util::read_vector_from_file(v, params + "E_increment");
    // std::cout << v.size() << std::endl;
    for (unsigned int i = 0; i < v.size(); ++i) {
        ftype ref = v[i];
        ftype real = Context::RfP->E_increment[i];
        ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)));
    }
}

TEST_F(testRFP, test_phi_s)
{
    std::vector<ftype> v;
    util::read_vector_from_file(v, params + "phi_s");
    // std::cout << v.size() << std::endl;
    for (unsigned int i = 0; i < v.size(); ++i) {
        ftype ref = v[i];
        ftype real = Context::RfP->phi_s[i];
        ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)));
    }
}

TEST_F(testRFP, test_Qs)
{
    std::vector<ftype> v;
    util::read_vector_from_file(v, params + "Qs");
    // std::cout << v.size() << std::endl;
    for (unsigned int i = 0; i < v.size(); ++i) {
        ftype ref = v[i];
        ftype real = Context::RfP->Qs[i];
        ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)));
    }
}

TEST_F(testRFP, test_omega_s0)
{
    std::vector<ftype> v;
    util::read_vector_from_file(v, params + "omega_s0");
    // std::cout << v.size() << std::endl;
    for (unsigned int i = 0; i < v.size(); ++i) {
        ftype ref = v[i];
        ftype real = Context::RfP->omega_s0[i];
        ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)));
    }
}

TEST_F(testRFP, test_omega_RF_d)
{
    std::vector<ftype> v;
    util::read_vector_from_file(v, params + "omega_RF_d[0]");
    // std::cout << v.size() << std::endl;
    for (unsigned int i = 0; i < v.size(); ++i) {
        ftype ref = v[i];
        ftype real = Context::RfP->omega_RF_d[0][i];
        ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)));
    }
}

TEST_F(testRFP, test_omega_RF)
{
    std::vector<ftype> v;
    util::read_vector_from_file(v, params + "omega_RF[0]");
    // std::cout << v.size() << std::endl;
    for (unsigned int i = 0; i < v.size(); ++i) {
        ftype ref = v[i];
        ftype real = Context::RfP->omega_RF[0][i];
        ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)));
    }
}

TEST_F(testRFP, test_t_RF)
{
    std::vector<ftype> v;
    util::read_vector_from_file(v, params + "t_RF");
    // std::cout << v.size() << std::endl;
    for (unsigned int i = 0; i < v.size(); ++i) {
        ftype ref = v[i];
        ftype real = Context::RfP->t_RF[i];
        ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)));
    }
}

int main(int ac, char *av[])
{
    ::testing::InitGoogleTest(&ac, av);
    return RUN_ALL_TESTS();
}

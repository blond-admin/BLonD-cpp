#include <iostream>
#include <blond/beams/Distributions.h>
#include <blond/input_parameters/GeneralParameters.h>
#include <blond/math_functions.h>
#include <blond/utilities.h>
#include <gtest/gtest.h>


class testSlices : public ::testing::Test {

protected:
    const long long N_b = 1e8;  // Intensity
    const ftype tau_0 = 1e-8; // Initial bunch length, 4 sigma [s]
// Machine and RF parameters
    const ftype C = 26658.883;       // Machine circumference [m]
    const long long p_i = 450e9;     // Synchronous momentum [eV/c]
    const ftype p_f = 455e9;     // Synchronous momentum, final
    const long long h = 35640;       // Harmonic number
    const ftype V = 7e6;             // RF voltage [V]
    const ftype dphi = 0;            // Phase modulation/offset
    const ftype gamma_t = 55.759505; // Transition gamma
    const ftype alpha =
        1.0 / gamma_t / gamma_t; // First order mom. comp. factor
    const int alpha_order = 1;
    const int n_sections = 1;
    // Tracking details

    const int N_t = 2000; // Number of turns to track
    const int N_p = 10000;  // Macro-particles
    const int N_slices = 100;
    virtual void SetUp()
    {
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

        longitudinal_bigaussian(tau_0 / 4, 0, -1, false);

    }

    virtual void TearDown()
    {
        // Code here will be called immediately after each test
        // (right before the destructor).
        delete Context::GP;
        delete Context::Beam;
        delete Context::RfP;
        // delete Context::Slice;
    }

};


TEST_F(testSlices, set_cuts1)
{
    auto Slice = new Slices(N_slices);

    std::string params = TEST_FILES "/Slices/set_cuts1/";
    auto epsilon = 1e-8;
    f_vector_t v;

    util::read_vector_from_file(v, params + "cut_left.txt");
    auto ref = v[0];
    auto real = Slice->cut_left;
    ASSERT_NEAR(ref, real, epsilon * std::max(std::abs(ref), std::abs(real)))
            << "Testing of cut_left failed\n";

    util::read_vector_from_file(v, params + "cut_right.txt");
    ref = v[0];
    real = Slice->cut_right;
    ASSERT_NEAR(ref, real, epsilon * std::max(std::abs(ref), std::abs(real)))
            << "Testing of cut_right failed\n";

    util::read_vector_from_file(v, params + "bin_centers.txt");
    for (uint i = 0; i < v.size(); ++i) {
        ftype ref = v[i];
        ftype real = Slice->bin_centers[i];
        ASSERT_NEAR(ref, real, epsilon * std::max(std::abs(ref), std::abs(real)))
                << "Testing of bin_centers failed on i " << i << '\n';
    }
    delete Slice;
}


TEST_F(testSlices, set_cuts2)
{
    auto Slice = new Slices(N_slices, 10);

    std::string params = TEST_FILES "/Slices/set_cuts2/";
    auto epsilon = 1e-8;
    f_vector_t v;

    util::read_vector_from_file(v, params + "cut_left.txt");
    auto ref = v[0];
    auto real = Slice->cut_left;
    ASSERT_NEAR(ref, real, epsilon * std::max(std::abs(ref), std::abs(real)))
            << "Testing of cut_left failed\n";

    util::read_vector_from_file(v, params + "cut_right.txt");
    ref = v[0];
    real = Slice->cut_right;
    ASSERT_NEAR(ref, real, epsilon * std::max(std::abs(ref), std::abs(real)))
            << "Testing of cut_right failed\n";

    util::read_vector_from_file(v, params + "bin_centers.txt");
    for (uint i = 0; i < v.size(); ++i) {
        ftype ref = v[i];
        ftype real = Slice->bin_centers[i];
        ASSERT_NEAR(ref, real, epsilon * std::max(std::abs(ref), std::abs(real)))
                << "Testing of bin_centers failed on i " << i << '\n';
    }
    delete Slice;

}


TEST_F(testSlices, set_cuts3)
{
    auto Slice = new Slices(N_slices, 0, -1e-8, 1e8);

    std::string params = TEST_FILES "/Slices/set_cuts3/";
    auto epsilon = 1e-8;
    f_vector_t v;

    util::read_vector_from_file(v, params + "cut_left.txt");
    auto ref = v[0];
    auto real = Slice->cut_left;
    ASSERT_NEAR(ref, real, epsilon * std::max(std::abs(ref), std::abs(real)))
            << "Testing of cut_left failed\n";

    util::read_vector_from_file(v, params + "cut_right.txt");
    ref = v[0];
    real = Slice->cut_right;
    ASSERT_NEAR(ref, real, epsilon * std::max(std::abs(ref), std::abs(real)))
            << "Testing of cut_right failed\n";

    util::read_vector_from_file(v, params + "bin_centers.txt");
    for (uint i = 0; i < v.size(); ++i) {
        ftype ref = v[i];
        ftype real = Slice->bin_centers[i];
        ASSERT_NEAR(ref, real, epsilon * std::max(std::abs(ref), std::abs(real)))
                << "Testing of bin_centers failed on i " << i << '\n';
    }
    delete Slice;

}


TEST_F(testSlices, sort_particles1)
{
    auto Slice = new Slices(N_slices);

    auto Beam = Context::Beam;
    std::string params = TEST_FILES "/Slices/sort_particles1/";
    auto epsilon = 1e-8;

    // util::dump(Beam->dE, "dE\n");
    // util::dump(Beam->dt, "dt\n");
    // util::dump(Beam->id, "id\n");
    f_vector_t v;
    util::read_vector_from_file(v, params + "dE.txt");
    for (uint i = 0; i < v.size(); ++i) {
        auto ref = v[i];
        auto real = Beam->dE[i];
        ASSERT_NEAR(ref, real, epsilon * std::max(std::abs(ref), std::abs(real)))
                << "Testing of dE failed on i " << i << '\n';
    }

    util::read_vector_from_file(v, params + "dt.txt");
    for (uint i = 0; i < v.size(); ++i) {
        auto ref = v[i];
        auto real = Beam->dt[i];
        ASSERT_NEAR(ref, real, epsilon * std::max(std::abs(ref), std::abs(real)))
                << "Testing of dt failed on i " << i << '\n';
    }

    /*
    // #NOTE : removed because python's argsort is not
    // a stable sort -> undefined behaviour on duplicates
    util::read_vector_from_file(v, params + "id.txt");
    for (uint i = 0; i < v.size(); ++i) {
        auto ref = 1.0 * v[i];
        auto real = 1.0 * Beam->id[i];
        EXPECT_NEAR(ref, real, epsilon * std::max(std::abs(ref), std::abs(real)))
                << "Testing of id failed on i " << i << '\n';
    }
    */
    delete Slice;

}

TEST_F(testSlices, track1)
{
    auto Slice = new Slices(N_slices);

    std::string params = TEST_FILES "/Slices/track1/";
    auto epsilon = 1e-8;
    f_vector_t v;

    for (int i = 0; i < 10; i++) Slice->track();

    util::read_vector_from_file(v, params + "n_macroparticles.txt");
    for (unsigned int i = 0; i < v.size(); ++i) {
        auto ref = 1.0 * v[i];
        auto real = 1.0 * Slice->n_macroparticles[i];
        ASSERT_NEAR(ref, real, epsilon * std::max(std::abs(ref), std::abs(real)))
                << "Testing of n_macroparticles failed on i " << i << '\n';
    }
    delete Slice;

}


TEST_F(testSlices, beam_profile_derivative1)
{
    auto Slice = new Slices(N_slices);
    std::string params = TEST_FILES "/Slices/beam_profile_derivative1/";
    auto epsilon = 1e-8;
    f_vector_t v;

    Slice->track();
    f_vector_t x, derivative;
    Slice->beam_profile_derivative(x, derivative, "gradient");

    util::read_vector_from_file(v, params + "x.txt");
    ASSERT_EQ(x.size(), v.size());
    for (unsigned int i = 0; i < v.size(); ++i) {
        auto ref = v[i];
        auto real = x[i];
        ASSERT_NEAR(ref, real, epsilon * std::max(std::abs(ref), std::abs(real)))
                << "Testing of x failed on i " << i << '\n';
    }

    util::read_vector_from_file(v, params + "derivative.txt");
    ASSERT_EQ(derivative.size(), v.size());
    for (unsigned int i = 0; i < v.size(); ++i) {
        auto ref = v[i];
        auto real = derivative[i];
        ASSERT_NEAR(ref, real, epsilon * std::max(std::abs(ref), std::abs(real)))
                << "Testing of derivative failed on i " << i << '\n';
    }
    delete Slice;

}


TEST_F(testSlices, beam_profile_derivative2)
{
    auto Slice = new Slices(N_slices);
    std::string params = TEST_FILES "/Slices/beam_profile_derivative2/";
    auto epsilon = 1e-8;
    f_vector_t v;

    Slice->track();
    f_vector_t x, derivative;
    Slice->beam_profile_derivative(x, derivative, "filter1d");

    util::read_vector_from_file(v, params + "x.txt");
    ASSERT_EQ(x.size(), v.size());
    for (unsigned int i = 0; i < v.size(); ++i) {
        auto ref = v[i];
        auto real = x[i];
        ASSERT_NEAR(ref, real, epsilon * std::max(std::abs(ref), std::abs(real)))
                << "Testing of x failed on i " << i << '\n';
    }

    util::read_vector_from_file(v, params + "derivative.txt");

    ASSERT_EQ(derivative.size(), v.size());
    for (unsigned int i = 0; i < v.size(); ++i) {
        auto ref = v[i];
        auto real = derivative[i];
        ASSERT_NEAR(ref, real, epsilon * std::max(std::abs(ref), std::abs(real)))
                << "Testing of derivative failed on i " << i << '\n';
    }
    delete Slice;

}


TEST_F(testSlices, beam_profile_derivative3)
{
    auto Slice = new Slices(N_slices);
    std::string params = TEST_FILES "/Slices/beam_profile_derivative3/";
    auto epsilon = 1e-8;
    f_vector_t v;

    Slice->track();
    f_vector_t x, derivative;
    Slice->beam_profile_derivative(x, derivative, "diff");
    util::read_vector_from_file(v, params + "x.txt");
    ASSERT_EQ(x.size(), v.size());
    for (unsigned int i = 0; i < v.size(); ++i) {
        auto ref = v[i];
        auto real = x[i];
        ASSERT_NEAR(ref, real, epsilon * std::max(std::abs(ref), std::abs(real)))
                << "Testing of x failed on i " << i << '\n';
    }

    util::read_vector_from_file(v, params + "derivative.txt");
    ftype max = *max_element(v.begin(), v.end(), [](ftype i, ftype j) {
        return fabs(i) < fabs(j);
    });
    max = std::abs(max);
    ASSERT_EQ(derivative.size(), v.size());
    for (unsigned int i = 0; i < v.size(); ++i) {
        auto ref = v[i];
        auto real = derivative[i];
        ASSERT_NEAR(ref, real, epsilon * max)
                << "Testing of derivative failed on i " << i << '\n';
    }
    delete Slice;

}



int main(int ac, char *av[])
{
    python::initialize();
    ::testing::InitGoogleTest(&ac, av);
    auto ret = RUN_ALL_TESTS();
    python::finalize();

    return ret;
}

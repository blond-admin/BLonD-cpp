#include <iostream>
#include <blond/beams/Distributions.h>
#include <blond/input_parameters/GeneralParameters.h>
#include <blond/math_functions.h>
#include <blond/utilities.h>
#include <gtest/gtest.h>

using namespace std;

class testSlices : public ::testing::Test {

protected:
    const long long N_b = 1e8;  // Intensity
    const double tau_0 = 1e-8; // Initial bunch length, 4 sigma [s]
// Machine and RF parameters
    const double C = 26658.883;       // Machine circumference [m]
    const double p_i = 450e9;     // Synchronous momentum [eV/c]
    const double p_f = 455e9;     // Synchronous momentum, final
    const long long h = 35640;       // Harmonic number
    const double V = 7e6;             // RF voltage [V]
    const double dphi = 0;            // Phase modulation/offset
    const double gamma_t = 55.759505; // Transition gamma
    const double alpha =
        1.0 / gamma_t / gamma_t; // First order mom. comp. factor
    const int alpha_order = 1;
    const int n_sections = 1;
    // Tracking details

    const int N_t = 2000; // Number of turns to track
    const int N_p = 10000;  // Macro-particles
    const int N_slices = 100;
    virtual void SetUp()
    {
        omp_set_num_threads(1);

        f_vector_2d_t momentumVec(n_sections, f_vector_t(N_t + 1));
        for (auto &v : momentumVec)
            mymath::linspace(v.data(), p_i, p_f, N_t + 1);

        f_vector_2d_t alphaVec(n_sections, f_vector_t(alpha_order, alpha));

        f_vector_t CVec(n_sections, C);

        f_vector_2d_t hVec(n_sections, f_vector_t(N_t + 1, h));

        f_vector_2d_t voltageVec(n_sections, f_vector_t(N_t + 1, V));

        f_vector_2d_t dphiVec(n_sections, f_vector_t(N_t + 1, dphi));

        Context::GP = new GeneralParameters(N_t, CVec, alphaVec,
                                             momentumVec,
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
        // delete Context::Slice;
    }

};



TEST_F(testSlices, set_cuts1)
{
    auto RfP = Context::RfP;
    auto Beam = Context::Beam;

    auto Slice = new Slices(RfP, Beam, N_slices);

    string params = TEST_FILES "/Slices/set_cuts1/";
    auto epsilon = 1e-8;
    f_vector_t v;

    util::read_vector_from_file(v, params + "cut_left.txt");
    auto ref = v[0];
    auto real = Slice->cut_left;
    ASSERT_NEAR(ref, real, epsilon * max(abs(ref), abs(real)))
            << "Testing of cut_left failed\n";

    util::read_vector_from_file(v, params + "cut_right.txt");
    ref = v[0];
    real = Slice->cut_right;
    ASSERT_NEAR(ref, real, epsilon * max(abs(ref), abs(real)))
            << "Testing of cut_right failed\n";

    util::read_vector_from_file(v, params + "bin_centers.txt");
    for (uint i = 0; i < v.size(); ++i) {
        double ref = v[i];
        double real = Slice->bin_centers[i];
        ASSERT_NEAR(ref, real, epsilon * max(abs(ref), abs(real)))
                << "Testing of bin_centers failed on i " << i << '\n';
    }
    delete Slice;
}


TEST_F(testSlices, set_cuts2)
{
    auto RfP = Context::RfP;
    auto Beam = Context::Beam;

    auto Slice = new Slices(RfP, Beam, N_slices, 10);

    string params = TEST_FILES "/Slices/set_cuts2/";
    auto epsilon = 1e-8;
    f_vector_t v;

    util::read_vector_from_file(v, params + "cut_left.txt");
    auto ref = v[0];
    auto real = Slice->cut_left;
    ASSERT_NEAR(ref, real, epsilon * max(abs(ref), abs(real)))
            << "Testing of cut_left failed\n";

    util::read_vector_from_file(v, params + "cut_right.txt");
    ref = v[0];
    real = Slice->cut_right;
    ASSERT_NEAR(ref, real, epsilon * max(abs(ref), abs(real)))
            << "Testing of cut_right failed\n";

    util::read_vector_from_file(v, params + "bin_centers.txt");
    for (uint i = 0; i < v.size(); ++i) {
        double ref = v[i];
        double real = Slice->bin_centers[i];
        ASSERT_NEAR(ref, real, epsilon * max(abs(ref), abs(real)))
                << "Testing of bin_centers failed on i " << i << '\n';
    }
    delete Slice;

}


TEST_F(testSlices, set_cuts3)
{
    auto RfP = Context::RfP;
    auto Beam = Context::Beam;

    auto Slice = new Slices(RfP, Beam, N_slices, 0, -1e-8, 1e8);

    string params = TEST_FILES "/Slices/set_cuts3/";
    auto epsilon = 1e-8;
    f_vector_t v;

    util::read_vector_from_file(v, params + "cut_left.txt");
    auto ref = v[0];
    auto real = Slice->cut_left;
    ASSERT_NEAR(ref, real, epsilon * max(abs(ref), abs(real)))
            << "Testing of cut_left failed\n";

    util::read_vector_from_file(v, params + "cut_right.txt");
    ref = v[0];
    real = Slice->cut_right;
    ASSERT_NEAR(ref, real, epsilon * max(abs(ref), abs(real)))
            << "Testing of cut_right failed\n";

    util::read_vector_from_file(v, params + "bin_centers.txt");
    for (uint i = 0; i < v.size(); ++i) {
        double ref = v[i];
        double real = Slice->bin_centers[i];
        ASSERT_NEAR(ref, real, epsilon * max(abs(ref), abs(real)))
                << "Testing of bin_centers failed on i " << i << '\n';
    }
    delete Slice;

}


TEST_F(testSlices, sort_particles1)
{
    auto RfP = Context::RfP;
    auto Beam = Context::Beam;

    auto Slice = new Slices(RfP, Beam, N_slices);


    string params = TEST_FILES "/Slices/sort_particles1/";
    auto epsilon = 1e-7;

    // util::dump(Beam->dE, "dE\n");
    // util::dump(Beam->dt, "dt\n");
    // util::dump(Beam->id, "id\n");
    f_vector_t v;
    util::read_vector_from_file(v, params + "dE.txt");
    for (uint i = 0; i < v.size(); ++i) {
        auto ref = v[i];
        auto real = Beam->dE[i];
        ASSERT_NEAR(ref, real, epsilon * max(abs(ref), abs(real)))
                << "Testing of dE failed on i " << i << '\n';
    }

    util::read_vector_from_file(v, params + "dt.txt");
    for (uint i = 0; i < v.size(); ++i) {
        auto ref = v[i];
        auto real = Beam->dt[i];
        ASSERT_NEAR(ref, real, epsilon * max(abs(ref), abs(real)))
                << "Testing of dt failed on i " << i << '\n';
    }

    /*
    // #NOTE : removed because python's argsort is not
    // a stable sort -> undefined behaviour on duplicates
    util::read_vector_from_file(v, params + "id.txt");
    for (uint i = 0; i < v.size(); ++i) {
        auto ref = 1.0 * v[i];
        auto real = 1.0 * Beam->id[i];
        EXPECT_NEAR(ref, real, epsilon * max(abs(ref), abs(real)))
                << "Testing of id failed on i " << i << '\n';
    }
    */
    delete Slice;

}

TEST_F(testSlices, convert_coordinates1)
{
    auto RfP = Context::RfP;
    auto Beam = Context::Beam;

    auto slice = Slices(RfP, Beam, N_slices);
    auto ret = slice.convert_coordinates(1.0, Slices::cuts_unit_t::s);
    ASSERT_DOUBLE_EQ(1.0, ret);
}


TEST_F(testSlices, convert_coordinates2)
{
    auto RfP = Context::RfP;
    auto Beam = Context::Beam;
    auto slice = Slices(RfP, Beam, N_slices);
    auto a = slice.convert_coordinates(1.0, Slices::cuts_unit_t::rad);
    auto b = slice.convert_coordinates(2.0, Slices::cuts_unit_t::rad);
    ASSERT_DOUBLE_EQ(1.0 / 2.0, a / b);
}


TEST_F(testSlices, smooth_histogram1)
{
    auto RfP = Context::RfP;
    auto Beam = Context::Beam;

    auto epsilon = 1e-8;
    string params = TEST_FILES "/Slices/smooth_histogram1/";
    f_vector_t v;

    auto slice = Slices(RfP, Beam, N_slices);
    slice.slice_constant_space_histogram_smooth();

    util::read_vector_from_file(v, params + "n_macroparticles.txt");

    ASSERT_EQ(v.size(), slice.n_macroparticles.size());

    for (uint i = 0; i < v.size(); ++i) {
        double ref = v[i];
        double real = slice.n_macroparticles[i];
        ASSERT_NEAR(ref, real, epsilon * max(abs(ref), abs(real)))
                << "Testing of n_macroparticles failed on i " << i << '\n';
    }
}



TEST_F(testSlices, smooth_histogram2)
{
    auto epsilon = 1e-8;
    auto RfP = Context::RfP;
    auto Beam = Context::Beam;
    string params = TEST_FILES "/Slices/smooth_histogram2/";
    f_vector_t v;

    auto cut_left = Beam->dt.front();
    auto cut_right = Beam->dt.back();
    if (cut_left > cut_right) swap(cut_left, cut_right);

    auto slice = Slices(RfP, Beam, N_slices, 0, cut_left, cut_right);
    slice.slice_constant_space_histogram_smooth();

    util::read_vector_from_file(v, params + "n_macroparticles.txt");

    ASSERT_EQ(v.size(), slice.n_macroparticles.size());

    for (uint i = 0; i < v.size(); ++i) {
        double ref = v[i];
        double real = slice.n_macroparticles[i];
        ASSERT_NEAR(ref, real, epsilon * max(abs(ref), abs(real)))
                << "Testing of n_macroparticles failed on i " << i << '\n';
    }
}


TEST_F(testSlices, track1)
{
    auto RfP = Context::RfP;
    auto Beam = Context::Beam;

    auto Slice = new Slices(RfP, Beam, N_slices);

    string params = TEST_FILES "/Slices/track1/";
    auto epsilon = 1e-8;
    f_vector_t v;

    for (int i = 0; i < 10; i++) Slice->track();

    util::read_vector_from_file(v, params + "n_macroparticles.txt");
    for (unsigned int i = 0; i < v.size(); ++i) {
        int ref = v[i];
        int real = Slice->n_macroparticles[i];
        ASSERT_EQ(ref, real)
                << "Testing of n_macroparticles failed on i " << i << '\n';
    }
    delete Slice;

}


TEST_F(testSlices, beam_profile_derivative1)
{
    auto RfP = Context::RfP;
    auto Beam = Context::Beam;

    auto Slice = new Slices(RfP, Beam, N_slices);
    string params = TEST_FILES "/Slices/beam_profile_derivative1/";
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
        ASSERT_NEAR(ref, real, epsilon * max(abs(ref), abs(real)))
                << "Testing of x failed on i " << i << '\n';
    }

    util::read_vector_from_file(v, params + "derivative.txt");
    ASSERT_EQ(derivative.size(), v.size());
    for (unsigned int i = 0; i < v.size(); ++i) {
        auto ref = v[i];
        auto real = derivative[i];
        ASSERT_NEAR(ref, real, epsilon * max(abs(ref), abs(real)))
                << "Testing of derivative failed on i " << i << '\n';
    }
    delete Slice;

}


TEST_F(testSlices, beam_profile_derivative2)
{
    auto RfP = Context::RfP;
    auto Beam = Context::Beam;

    auto Slice = new Slices(RfP, Beam, N_slices);
    string params = TEST_FILES "/Slices/beam_profile_derivative2/";
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
        ASSERT_NEAR(ref, real, epsilon * max(abs(ref), abs(real)))
                << "Testing of x failed on i " << i << '\n';
    }

    util::read_vector_from_file(v, params + "derivative.txt");

    ASSERT_EQ(derivative.size(), v.size());
    for (unsigned int i = 0; i < v.size(); ++i) {
        auto ref = v[i];
        auto real = derivative[i];
        ASSERT_NEAR(ref, real, epsilon * max(abs(ref), abs(real)))
                << "Testing of derivative failed on i " << i << '\n';
    }
    delete Slice;

}


TEST_F(testSlices, beam_profile_derivative3)
{
    auto RfP = Context::RfP;
    auto Beam = Context::Beam;

    auto Slice = new Slices(RfP, Beam, N_slices);
    string params = TEST_FILES "/Slices/beam_profile_derivative3/";
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
        ASSERT_NEAR(ref, real, epsilon * max(abs(ref), abs(real)))
                << "Testing of x failed on i " << i << '\n';
    }

    util::read_vector_from_file(v, params + "derivative.txt");
    double max = *max_element(v.begin(), v.end(), [](double i, double j) {
        return fabs(i) < fabs(j);
    });
    max = abs(max);
    ASSERT_EQ(derivative.size(), v.size());
    for (unsigned int i = 0; i < v.size(); ++i) {
        auto ref = v[i];
        auto real = derivative[i];
        ASSERT_NEAR(ref, real, epsilon * max)
                << "Testing of derivative failed on i " << i << '\n';
    }
    delete Slice;

}

TEST_F(testSlices, beam_profile_derivative_deathtest1)
{
    auto RfP = Context::RfP;
    auto Beam = Context::Beam;

    auto Slice = new Slices(RfP, Beam, N_slices);
    auto epsilon = 1e-8;
    f_vector_t x, derivative;
    ASSERT_DEATH(Slice->beam_profile_derivative(x, derivative, "blabla"),
                 "Option for derivative is not recognized.\n");
    delete Slice;

}



TEST_F(testSlices, rms1)
{
    auto RfP = Context::RfP;
    auto Beam = Context::Beam;

    auto slice = Slices(RfP, Beam, N_slices);
    string params = TEST_FILES "/Slices/rms1/";
    auto epsilon = 1e-8;
    f_vector_t v;

    slice.track();
    slice.rms();

    util::read_vector_from_file(v, params + "bp_rms.txt");
    ASSERT_EQ(v.size(), 1);
    EXPECT_NEAR(v[0], slice.bp_rms,
                epsilon * max(abs(v[0]), abs(slice.bp_rms)))
            << "Testing failed on slice.bp_rms\n";


    util::read_vector_from_file(v, params + "bl_rms.txt");
    ASSERT_EQ(v.size(), 1);
    EXPECT_NEAR(v[0], slice.bl_rms,
                epsilon * max(abs(v[0]), abs(slice.bl_rms)))
            << "Testing failed on slice.bl_rms\n";

}


TEST_F(testSlices, rms2)
{
    auto RfP = Context::RfP;
    auto Beam = Context::Beam;
    string params = TEST_FILES "/Slices/rms2/";
    auto epsilon = 1e-8;
    f_vector_t v;

    auto cut_left = Beam->dt.front();
    auto cut_right = Beam->dt.back();
    if (cut_left > cut_right) swap(cut_left, cut_right);

    auto slice = Slices(RfP, Beam, N_slices, 0, cut_left, cut_right);

    slice.track();
    slice.rms();

    util::read_vector_from_file(v, params + "bp_rms.txt");
    ASSERT_EQ(v.size(), 1);
    EXPECT_NEAR(v[0], slice.bp_rms,
                epsilon * max(abs(v[0]), abs(slice.bp_rms)))
            << "Testing failed on slice.bp_rms\n";


    util::read_vector_from_file(v, params + "bl_rms.txt");
    ASSERT_EQ(v.size(), 1);
    EXPECT_NEAR(v[0], slice.bl_rms,
                epsilon * max(abs(v[0]), abs(slice.bl_rms)))
            << "Testing failed on slice.bl_rms\n";

}


TEST_F(testSlices, fwhm1)
{
    auto RfP = Context::RfP;
    auto Beam = Context::Beam;

    auto slice = Slices(RfP, Beam, N_slices);
    string params = TEST_FILES "/Slices/fwhm1/";
    auto epsilon = 1e-8;
    f_vector_t v;

    slice.track();
    slice.fwhm();

    util::read_vector_from_file(v, params + "bp_fwhm.txt");
    ASSERT_EQ(v.size(), 1);
    EXPECT_NEAR(v[0], slice.bp_fwhm,
                epsilon * max(abs(v[0]), abs(slice.bp_fwhm)))
            << "Testing failed on slice.bp_fwhm\n";


    util::read_vector_from_file(v, params + "bl_fwhm.txt");
    ASSERT_EQ(v.size(), 1);
    EXPECT_NEAR(v[0], slice.bl_fwhm,
                epsilon * max(abs(v[0]), abs(slice.bl_fwhm)))
            << "Testing failed on slice.bl_fwhm\n";

}


TEST_F(testSlices, fwhm2)
{
    auto RfP = Context::RfP;
    auto Beam = Context::Beam;
    string params = TEST_FILES "/Slices/fwhm2/";
    auto epsilon = 1e-8;
    f_vector_t v;

    auto cut_left = Beam->dt.front();
    auto cut_right = Beam->dt.back();
    if (cut_left > cut_right) swap(cut_left, cut_right);

    auto slice = Slices(RfP, Beam, N_slices, 0, cut_left, cut_right);

    slice.track();
    slice.fwhm(10.);
    util::read_vector_from_file(v, params + "bp_fwhm.txt");
    ASSERT_EQ(v.size(), 1);
    EXPECT_NEAR(v[0], slice.bp_fwhm,
                epsilon * max(abs(v[0]), abs(slice.bp_fwhm)))
            << "Testing failed on slice.bp_fwhm\n";


    util::read_vector_from_file(v, params + "bl_fwhm.txt");
    ASSERT_EQ(v.size(), 1);
    EXPECT_NEAR(v[0], slice.bl_fwhm,
                epsilon * max(abs(v[0]), abs(slice.bl_fwhm)))
            << "Testing failed on slice.bl_fwhm\n";

}


TEST_F(testSlices, fwhm3)
{
    auto RfP = Context::RfP;
    auto Beam = Context::Beam;

    string params = TEST_FILES "/Slices/fwhm3/";
    auto epsilon = 1e-8;
    f_vector_t v;

    auto cut_left = Beam->dt.front();
    auto cut_right = Beam->dt.back();
    if (cut_left > cut_right) swap(cut_left, cut_right);

    auto slice = Slices(RfP, Beam, N_slices, 0, cut_left, cut_right);

    slice.track();
    slice.fwhm(0.1);
    util::read_vector_from_file(v, params + "bp_fwhm.txt");
    ASSERT_EQ(v.size(), 1);
    EXPECT_NEAR(v[0], slice.bp_fwhm,
                epsilon * max(abs(v[0]), abs(slice.bp_fwhm)))
            << "Testing failed on slice.bp_fwhm\n";


    util::read_vector_from_file(v, params + "bl_fwhm.txt");
    ASSERT_EQ(v.size(), 1);
    EXPECT_NEAR(v[0], slice.bl_fwhm,
                epsilon * max(abs(v[0]), abs(slice.bl_fwhm)))
            << "Testing failed on slice.bl_fwhm\n";

}


TEST_F(testSlices, gaussian_fit1)
{
    auto RfP = Context::RfP;
    auto Beam = Context::Beam;

    auto slice = Slices(RfP, Beam, N_slices, 0, 0, 0, Slices::cuts_unit_t::s,
                        Slices::fit_t::gaussian);
    string params = TEST_FILES "/Slices/gaussian_fit1/";
    auto epsilon = 1e-8;
    f_vector_t v;

    slice.track();
    slice.gaussian_fit();
    // std::cout << "bl_gauss: " << slice.bl_gauss << "\n";
    // std::cout << "bp_gauss: " << slice.bp_gauss << "\n";
    util::read_vector_from_file(v, params + "bp_gauss.txt");
    ASSERT_EQ(v.size(), 1);
    EXPECT_NEAR(v[0], slice.bp_gauss,
                epsilon * max(abs(v[0]), abs(slice.bp_gauss)))
            << "Testing failed on slice.bp_gauss\n";


    util::read_vector_from_file(v, params + "bl_gauss.txt");
    ASSERT_EQ(v.size(), 1);
    EXPECT_NEAR(v[0], slice.bl_gauss,
                epsilon * max(abs(v[0]), abs(slice.bl_gauss)))
            << "Testing failed on slice.bl_fwhm\n";

}


TEST_F(testSlices, gaussian_fit2)
{
    auto RfP = Context::RfP;
    auto Beam = Context::Beam;

    auto slice = Slices(RfP, Beam, N_slices, 0, 0, 0, Slices::cuts_unit_t::s,
                        Slices::fit_t::gaussian);
    string params = TEST_FILES "/Slices/gaussian_fit2/";
    auto epsilon = 1e-8;
    f_vector_t v;

    slice.track();
    slice.gaussian_fit();
    slice.gaussian_fit();
    // std::cout << "bl_gauss: " << slice.bl_gauss << "\n";
    // std::cout << "bp_gauss: " << slice.bp_gauss << "\n";
    util::read_vector_from_file(v, params + "bp_gauss.txt");
    ASSERT_EQ(v.size(), 1);
    EXPECT_NEAR(v[0], slice.bp_gauss,
                epsilon * max(abs(v[0]), abs(slice.bp_gauss)))
            << "Testing failed on slice.bp_gauss\n";


    util::read_vector_from_file(v, params + "bl_gauss.txt");
    ASSERT_EQ(v.size(), 1);
    EXPECT_NEAR(v[0], slice.bl_gauss,
                epsilon * max(abs(v[0]), abs(slice.bl_gauss)))
            << "Testing failed on slice.bl_fwhm\n";

}

TEST_F(testSlices, gaussian_fit_deathtest1)
{
    auto RfP = Context::RfP;
    auto Beam = Context::Beam;

    auto slice = Slices(RfP, Beam, N_slices, 0, 0, 0, Slices::cuts_unit_t::s,
                        Slices::fit_t::gaussian);
    slice.n_macroparticles.clear();
    ASSERT_DEATH(slice.gaussian_fit(), "[gaussian_fit]\\s*");
}


TEST_F(testSlices, chebyshev1)
{
    auto RfP = Context::RfP;
    auto Beam = Context::Beam;

    auto slice = Slices(RfP, Beam, N_slices);
    string params = TEST_FILES "/Slices/chebyshev1/";
    auto epsilon = 1e-8;
    f_vector_t v;

    slice.track();
    map<string, string> filter_option = {
        {"type", "chebyshev"},
        {"pass_frequency", "0.2"},
        {"stop_frequency", "0.3"},
        {"gain_pass", "3"},
        {"gain_stop", "60"}
    };
    int nCoeff;
    f_vector_t a, b;
    slice.beam_profile_filter_chebyshev(filter_option, nCoeff, b, a);

    util::read_vector_from_file(v, params + "n_coefficients.txt");
    ASSERT_EQ(v.size(), 1);
    ASSERT_EQ(v[0], nCoeff) << "Testing failed on nCoeff\n";


    util::read_vector_from_file(v, params + "a.txt");
    ASSERT_EQ(v.size(), a.size());
    for (uint i = 0; i < v.size(); ++i) {
        double ref = v[i];
        double real = a[i];
        ASSERT_NEAR(ref, real, epsilon * max(abs(ref), abs(real)))
                << "Testing of a failed on i " << i << '\n';
    }

    util::read_vector_from_file(v, params + "b.txt");
    ASSERT_EQ(v.size(), b.size());
    for (uint i = 0; i < v.size(); ++i) {
        double ref = v[i];
        double real = b[i];
        ASSERT_NEAR(ref, real, epsilon * max(abs(ref), abs(real)))
                << "Testing of b failed on i " << i << '\n';
    }

}


TEST_F(testSlices, chebyshev2)
{
    auto RfP = Context::RfP;
    auto Beam = Context::Beam;

    auto slice = Slices(RfP, Beam, N_slices);
    string params = TEST_FILES "/Slices/chebyshev2/";
    auto epsilon = 1e-8;
    f_vector_t v;

    slice.track();
    map<string, string> filter_option = {
        {"type", "chebyshev"},
        {"pass_frequency", "0.6"},
        {"stop_frequency", "0.2"},
        {"gain_pass", "2"},
        {"gain_stop", "40"}
    };
    int nCoeff;
    f_vector_t a, b;
    slice.beam_profile_filter_chebyshev(filter_option, nCoeff, b, a);

    util::read_vector_from_file(v, params + "n_coefficients.txt");
    ASSERT_EQ(v.size(), 1);
    ASSERT_EQ(v[0], nCoeff) << "Testing failed on nCoeff\n";


    util::read_vector_from_file(v, params + "a.txt");
    ASSERT_EQ(v.size(), a.size());
    for (uint i = 0; i < v.size(); ++i) {
        double ref = v[i];
        double real = a[i];
        ASSERT_NEAR(ref, real, epsilon * max(abs(ref), abs(real)))
                << "Testing of a failed on i " << i << '\n';
    }

    util::read_vector_from_file(v, params + "b.txt");
    ASSERT_EQ(v.size(), b.size());
    for (uint i = 0; i < v.size(); ++i) {
        double ref = v[i];
        double real = b[i];
        ASSERT_NEAR(ref, real, epsilon * max(abs(ref), abs(real)))
                << "Testing of b failed on i " << i << '\n';
    }

}


TEST_F(testSlices, chebyshev3)
{
    auto RfP = Context::RfP;
    auto Beam = Context::Beam;

    auto slice = Slices(RfP, Beam, N_slices);
    string params = TEST_FILES "/Slices/chebyshev3/";
    auto epsilon = 1e-8;
    f_vector_t v;

    slice.track();
    map<string, string> filter_option = {
        {"type", "chebyshev"},
        {"pass_frequency", "0.1"},
        {"stop_frequency", "0.4"},
        {"gain_pass", "3"},
        {"gain_stop", "50"},
        {"transfer_function_plot", "true"}
    };
    int nCoeff;
    f_vector_t transferFreq;
    complex_vector_t transferGain;
    slice.beam_profile_filter_chebyshev(filter_option,
                                        nCoeff, transferFreq,
                                        transferGain);

    util::read_vector_from_file(v, params + "n_coefficients.txt");
    ASSERT_EQ(v.size(), 1);
    ASSERT_EQ(v[0], nCoeff) << "Testing failed on nCoeff\n";


    util::read_vector_from_file(v, params + "transferFreq.txt");
    ASSERT_EQ(v.size(), transferFreq.size());
    for (uint i = 0; i < v.size(); ++i) {
        double ref = v[i];
        double real = transferFreq[i];
        ASSERT_NEAR(ref, real, epsilon * max(abs(ref), abs(real)))
                << "Testing of transferFreq failed on i " << i << '\n';
    }

    util::read_vector_from_file(v, params + "transferGain.txt");
    ASSERT_EQ(v.size(), transferGain.size());
    for (uint i = 0; i < v.size(); ++i) {
        double ref = v[i];
        double real = abs(transferGain[i]);
        ASSERT_NEAR(ref, real, epsilon * max(abs(ref), abs(real)))
                << "Testing of transferGain failed on i " << i << '\n';
    }

}

TEST_F(testSlices, chebyshev4)
{
    auto RfP = Context::RfP;
    auto Beam = Context::Beam;

    auto slice = Slices(RfP, Beam, N_slices);
    string params = TEST_FILES "/Slices/chebyshev4/";
    auto epsilon = 1e-8;
    f_vector_t v;

    slice.track();
    map<string, string> filter_option = {
        {"type", "chebyshev"},
        {"pass_frequency", "0.1"},
        {"stop_frequency", "1.1"},
        {"gain_pass", "30"},
        {"gain_stop", "35"},
        {"transfer_function_plot", "true"}
    };
    int nCoeff;
    f_vector_t transferFreq;
    complex_vector_t transferGain;
    slice.beam_profile_filter_chebyshev(filter_option,
                                        nCoeff, transferFreq,
                                        transferGain);

    util::read_vector_from_file(v, params + "n_coefficients.txt");
    ASSERT_EQ(v.size(), 1);
    ASSERT_EQ(v[0], nCoeff) << "Testing failed on nCoeff\n";


    util::read_vector_from_file(v, params + "transferFreq.txt");
    ASSERT_EQ(v.size(), transferFreq.size());
    for (uint i = 0; i < v.size(); ++i) {
        double ref = v[i];
        double real = transferFreq[i];
        ASSERT_NEAR(ref, real, epsilon * max(abs(ref), abs(real)))
                << "Testing of transferFreq failed on i " << i << '\n';
    }

    util::read_vector_from_file(v, params + "transferGain.txt");
    ASSERT_EQ(v.size(), transferGain.size());
    for (uint i = 0; i < v.size(); ++i) {
        double ref = v[i];
        double real = abs(transferGain[i]);
        ASSERT_NEAR(ref, real, epsilon * max(abs(ref), abs(real)))
                << "Testing of transferGain failed on i " << i << '\n';
    }

}

TEST_F(testSlices, chebyshev_deathtest1)
{
    auto RfP = Context::RfP;
    auto Beam = Context::Beam;

    auto slice = Slices(RfP, Beam, N_slices);
    auto epsilon = 1e-8;
    f_vector_t v;

    slice.track();
    map<string, string> filter_option = {
        {"type", "chebyshev"},
        {"pass_frequency", "0.1"},
        {"stop_frequency", "1.1"},
        {"gain_pass", "30"},
        {"gain_stop", "35"},
        {"transfer_function_plot", "false"}
    };
    int nCoeff;
    f_vector_t transferFreq;
    complex_vector_t transferGain;

    ASSERT_DEATH(slice.beam_profile_filter_chebyshev(filter_option,
                 nCoeff, transferFreq, transferGain),
                 "[beam_profile_filter_chebyshev] \\s*");

}


TEST_F(testSlices, chebyshev_deathtest2)
{
    auto RfP = Context::RfP;
    auto Beam = Context::Beam;

    auto slice = Slices(RfP, Beam, N_slices);
    auto epsilon = 1e-8;
    f_vector_t v;

    slice.track();
    map<string, string> filter_option = {
        {"type", "chebyshev"},
        {"pass_frequency", "0.1"},
        {"stop_frequency", "1.1"},
        {"gain_pass", "30"},
        {"gain_stop", "35"},
        {"transfer_function_plot", "true"}
    };
    int nCoeff;
    f_vector_t a;
    f_vector_t b;

    ASSERT_DEATH(slice.beam_profile_filter_chebyshev(filter_option,
                 nCoeff, a, b),
                 "[beam_profile_filter_chebyshev] \\s*");

}


TEST_F(testSlices, chebyshev_deathtest3)
{
    auto RfP = Context::RfP;
    auto Beam = Context::Beam;

    auto slice = Slices(RfP, Beam, N_slices);
    auto epsilon = 1e-8;
    f_vector_t v;

    slice.track();
    map<string, string> filter_option = {
        {"type", "chebyshev"},
    };
    int nCoeff;
    f_vector_t a;
    f_vector_t b;

    ASSERT_DEATH(slice.beam_profile_filter_chebyshev(filter_option,
                 nCoeff, a, b),
                 "[beam_profile_filter_chebyshev] \\s*");

}

TEST_F(testSlices, chebyshev_deathtest4)
{
    auto RfP = Context::RfP;
    auto Beam = Context::Beam;
    auto slice = Slices(RfP, Beam, N_slices);
    auto epsilon = 1e-8;
    f_vector_t v;

    slice.track();
    map<string, string> filter_option = {
        {"type", "chebyshev"},
    };
    int nCoeff;
    f_vector_t a;
    complex_vector_t b;

    ASSERT_DEATH(slice.beam_profile_filter_chebyshev(filter_option,
                 nCoeff, a, b),
                 "[beam_profile_filter_chebyshev] \\s*");

}

// TEST_F(testSlices, gaussian_fit3)
// {
//     auto slice = Slices(RfP, Beam, N_slices, 0, 0, 0, Slices::cuts_unit_t::s,
//                         Slices::fit_t::gaussian);
//     string params = TEST_FILES "/Slices/gaussian_fit3/";
//     auto epsilon = 1e-6;
//     f_vector_t v;

//     slice.bp_gauss = 100e-9;
//     slice.bl_gauss = 150e-8;
//     slice.track();
//     slice.gaussian_fit();
//     // std::cout << "bl_gauss: " << slice.bl_gauss << "\n";
//     // std::cout << "bp_gauss: " << slice.bp_gauss << "\n";
//     util::read_vector_from_file(v, params + "bp_gauss.txt");
//     ASSERT_EQ(v.size(), 1);
//     EXPECT_NEAR(v[0], slice.bp_gauss,
//                 epsilon * max(abs(v[0]), abs(slice.bp_gauss)))
//             << "Testing failed on slice.bp_gauss\n";


//     util::read_vector_from_file(v, params + "bl_gauss.txt");
//     ASSERT_EQ(v.size(), 1);
//     EXPECT_NEAR(v[0], slice.bl_gauss,
//                 epsilon * max(abs(v[0]), abs(slice.bl_gauss)))
//             << "Testing failed on slice.bl_fwhm\n";

// }


int main(int ac, char *av[])
{
    python::initialize();
    ::testing::InitGoogleTest(&ac, av);
    auto ret = RUN_ALL_TESTS();
    python::finalize();

    return ret;
}

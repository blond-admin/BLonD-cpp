#include <iostream>

#include <blond/beams/Distributions.h>
#include <blond/input_parameters/GeneralParameters.h>
#include <blond/math_functions.h>
#include <blond/trackers/Tracker.h>
#include <blond/utilities.h>
#include <gtest/gtest.h>


class testTracker : public ::testing::Test {

protected:
    const long long N_b = 1e9;  // Intensity
    const ftype tau_0 = 0.4e-9; // Initial bunch length, 4 sigma [s]
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
    const int N_p = 1000;  // Macro-particles
    const int N_slices = 10;

    virtual void SetUp()
    {
        omp_set_num_threads(1);

        f_vector_2d_t momentumVec(n_sections, f_vector_t(N_t + 1));
        for (auto &v : momentumVec)
            mymath::linspace(v.data(), p_i, p_f, N_t + 1);

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

        longitudinal_bigaussian(tau_0 / 4, 0, -1, false);

        Context::Slice = new Slices(N_slices);
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



TEST_F(testTracker, kick1)
{
    auto Beam = Context::Beam;
    auto epsilon = 1e-8;
    std::string params = TEST_FILES "/Tracker/kick1/";

    auto long_tracker = new RingAndRfSection();

    for (int i = 0; i < 10; i++) long_tracker->kick(Beam->dt, Beam->dE, i);

    f_vector_t v;
    util::read_vector_from_file(v, params + "dE.txt");
    ASSERT_EQ(v.size(), Beam->dE.size());
    for (uint i = 0; i < v.size(); ++i) {
        auto ref = v[i];
        auto real = Beam->dE[i];
        ASSERT_NEAR(ref, real, epsilon * std::max(std::abs(ref), std::abs(real)))
                << "Testing of dE failed on i " << i << '\n';

    }

    delete long_tracker;
}


TEST_F(testTracker, drift1)
{
    auto Beam = Context::Beam;
    auto epsilon = 1e-8;
    std::string params = TEST_FILES "/Tracker/drift1/";

    auto long_tracker = new RingAndRfSection();

    for (int i = 0; i < 10; i++) long_tracker->drift(Beam->dt, Beam->dE, i);

    f_vector_t v;
    util::read_vector_from_file(v, params + "dt.txt");
    ASSERT_EQ(v.size(), Beam->dt.size());
    for (uint i = 0; i < v.size(); ++i) {
        auto ref = v[i];
        auto real = Beam->dt[i];
        ASSERT_NEAR(ref, real, epsilon * std::max(std::abs(ref), std::abs(real)))
                << "Testing of dt failed on i " << i << '\n';

    }

    delete long_tracker;
}


TEST_F(testTracker, track1)
{

    auto epsilon = 1e-8;
    std::string params = TEST_FILES "/Tracker/track1/";

    auto long_tracker = new RingAndRfSection();
    long_tracker->track();

    f_vector_t v;
    util::read_vector_from_file(v, params + "dE.txt");
    ASSERT_EQ(v.size(), Context::Beam->dE.size());
    for (uint i = 0; i < v.size(); ++i) {
        auto ref = v[i];
        auto real = Context::Beam->dE[i];
        ASSERT_NEAR(ref, real, epsilon * std::max(std::abs(ref), std::abs(real)))
                << "Testing of dE failed on i " << i << '\n';

    }

    util::read_vector_from_file(v, params + "dt.txt");
    ASSERT_EQ(v.size(), Context::Beam->dt.size());
    for (uint i = 0; i < v.size(); ++i) {
        auto ref = v[i];
        auto real = Context::Beam->dt[i];
        ASSERT_NEAR(ref, real, epsilon * std::max(std::abs(ref), std::abs(real)))
                << "Testing of dt failed on i " << i << '\n';

    }

    delete long_tracker;
}

TEST_F(testTracker, track2)
{

    auto epsilon = 1e-8;
    std::string params = TEST_FILES "/Tracker/track2/";

    auto long_tracker = new RingAndRfSection();

    for (int i = 0; i < 10; i++) long_tracker->track();

    f_vector_t v;
    util::read_vector_from_file(v, params + "dE.txt");
    ASSERT_EQ(v.size(), Context::Beam->dE.size());
    for (uint i = 0; i < v.size(); ++i) {
        auto ref = v[i];
        auto real = Context::Beam->dE[i];
        ASSERT_NEAR(ref, real, epsilon * std::max(std::abs(ref), std::abs(real)))
                << "Testing of dE failed on i " << i << '\n';

    }

    util::read_vector_from_file(v, params + "dt.txt");
    ASSERT_EQ(v.size(), Context::Beam->dt.size());
    for (uint i = 0; i < v.size(); ++i) {
        auto ref = v[i];
        auto real = Context::Beam->dt[i];
        ASSERT_NEAR(ref, real, epsilon * std::max(std::abs(ref), std::abs(real)))
                << "Testing of dt failed on i " << i << '\n';

    }

    delete long_tracker;
}


TEST_F(testTracker, rf_voltage_calculation1)
{
    auto epsilon = 1e-8;
    std::string params = TEST_FILES "/Tracker/rf_voltage_calculation1/";
    auto Slice = Context::Slice;
    auto long_tracker = new RingAndRfSection();

    long_tracker->rf_voltage_calculation(0, Slice);

    f_vector_t v;
    util::read_vector_from_file(v, params + "rf_voltage.txt");
    ASSERT_EQ(v.size(), long_tracker->fRfVoltage.size());
    for (uint i = 0; i < v.size(); ++i) {
        auto ref = v[i];
        auto real = long_tracker->fRfVoltage[i];
        ASSERT_NEAR(ref, real, epsilon * std::max(std::abs(ref), std::abs(real)))
                << "Testing of rf_voltage failed on i " << i << '\n';

    }

    delete long_tracker;
}


TEST_F(testTracker, rf_voltage_calculation2)
{
    auto epsilon = 1e-8;
    std::string params = TEST_FILES "/Tracker/rf_voltage_calculation2/";
    auto Slice = Context::Slice;
    auto Beam = Context::Beam;
    auto long_tracker = new RingAndRfSection();

    for (int i = 0; i < 10; i++) {
        Slice->track();
        long_tracker->track();
        Beam->statistics();
        Slice->track_cuts();
        long_tracker->rf_voltage_calculation(i, Slice);
    }

    f_vector_t v;
    util::read_vector_from_file(v, params + "rf_voltage.txt");
    ASSERT_EQ(v.size(), long_tracker->fRfVoltage.size());
    for (uint i = 0; i < v.size(); ++i) {
        auto ref = v[i];
        auto real = long_tracker->fRfVoltage[i];
        ASSERT_NEAR(ref, real, epsilon * std::max(std::abs(ref), std::abs(real)))
                << "Testing of rf_voltage failed on i " << i << '\n';

    }

    delete long_tracker;
}


TEST_F(testTracker, rf_voltage_calculation3)
{
    auto epsilon = 1e-8;
    std::string params = TEST_FILES "/Tracker/rf_voltage_calculation3/";
    auto Slice = Context::Slice;
    auto Beam = Context::Beam;
    auto long_tracker = new RingAndRfSection();

    for (int i = 0; i < N_t; i++) {
        Slice->track();
        long_tracker->track();
        Beam->statistics();
        Slice->track_cuts();
        long_tracker->rf_voltage_calculation(i, Slice);
    }

    f_vector_t v;
    util::read_vector_from_file(v, params + "rf_voltage.txt");
    ASSERT_EQ(v.size(), long_tracker->fRfVoltage.size());
    for (uint i = 0; i < v.size(); ++i) {
        auto ref = v[i];
        auto real = long_tracker->fRfVoltage[i];
        ASSERT_NEAR(ref, real, epsilon * std::max(std::abs(ref), std::abs(real)))
                << "Testing of rf_voltage failed on i " << i << '\n';

    }

    delete long_tracker;
}



class testTrackerMultiRf : public ::testing::Test {

protected:
    const long long N_b = 1e9;  // Intensity
    const ftype tau_0 = 0.4e-9; // Initial bunch length, 4 sigma [s]
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
    const int n_rf = 3;
    // Tracking details

    const int N_t = 2000; // Number of turns to track
    const int N_p = 1000;  // Macro-particles
    const int N_slices = 100;

    virtual void SetUp()
    {
        omp_set_num_threads(1);

        f_vector_2d_t momentumVec(n_sections, f_vector_t(N_t + 1));
        for (auto &v : momentumVec)
            mymath::linspace(v.data(), p_i, p_f, N_t + 1);

        f_vector_2d_t alphaVec(n_sections, f_vector_t(alpha_order + 1, alpha));

        f_vector_t CVec(n_sections, C);

        f_vector_2d_t hVec(n_rf);
        for (int i = 0; i < n_rf; i++)
            hVec[i] = f_vector_t(N_t + 1, (i + 1) * h);

        f_vector_2d_t voltageVec(n_rf);

        for (int i = 0; i < n_rf; i++)
            voltageVec[i] = f_vector_t(N_t + 1, (i + 1) * V);

        f_vector_2d_t dphiVec(n_rf, f_vector_t(N_t + 1, dphi));

        Context::GP = new GeneralParameters(N_t, CVec, alphaVec,
                                            alpha_order, momentumVec,
                                            GeneralParameters::particle_t::proton);

        Context::Beam = new Beams(N_p, N_b);

        Context::RfP = new RfParameters(n_rf, hVec, voltageVec, dphiVec);

        longitudinal_bigaussian(tau_0 / 4, 0, -1, false);

        Context::Slice = new Slices(N_slices);
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


TEST_F(testTrackerMultiRf, rf_voltage_calculation4)
{
    auto epsilon = 1e-8;
    std::string params = TEST_FILES "/Tracker/rf_voltage_calculation4/";
    auto Slice = Context::Slice;
    auto Beam = Context::Beam;
    auto long_tracker = new RingAndRfSection();

    for (int i = 0; i < N_t; i++) {
        Slice->track();
        long_tracker->track();
        Beam->statistics();
        Slice->track_cuts();
        long_tracker->rf_voltage_calculation(i, Slice);
    }

    f_vector_t v;
    util::read_vector_from_file(v, params + "rf_voltage.txt");
    ASSERT_EQ(v.size(), long_tracker->fRfVoltage.size());
    for (uint i = 0; i < v.size(); ++i) {
        auto ref = v[i];
        auto real = long_tracker->fRfVoltage[i];
        ASSERT_NEAR(ref, real, epsilon * std::max(std::abs(ref), std::abs(real)))
                << "Testing of rf_voltage failed on i " << i << '\n';

    }

    delete long_tracker;
}


class testTrackerPeriodicity : public ::testing::Test {

protected:

    const long long N_b = 1e9;  // Intensity
    const ftype tau_0 = 1e-9; // Initial bunch length, 4 sigma [s]
    // Machine and RF parameters
    const ftype C = 26658.883;       // Machine circumference [m]
    const long long p_i = 455e9;     // Synchronous momentum [eV/c]
    const ftype p_f = 460.005e9;     // Synchronous momentum, final
    const long long h = 35640;       // Harmonic number
    const ftype V = 10e6;             // RF voltage [V]
    const ftype dphi = 0;            // Phase modulation/offset
    const ftype gamma_t = 55.759505; // Transition gamma
    const ftype alpha =
        1.0 / gamma_t / gamma_t; // First order mom. comp. factor
    const int alpha_order = 1;
    const int n_sections = 1;
    // Tracking details

    const int N_t = 1111; // Number of turns to track
    const int N_p = 7777;  // Macro-particles
    const int N_slices = 100;


    virtual void SetUp()
    {
        omp_set_num_threads(1);

        f_vector_2d_t momentumVec(n_sections, f_vector_t(N_t + 1));
        for (auto &v : momentumVec)
            mymath::linspace(v.data(), p_i, p_f, N_t + 1);

        f_vector_2d_t alphaVec(alpha_order + 1, f_vector_t(n_sections, alpha));

        f_vector_t CVec(n_sections, C);

        f_vector_2d_t hVec(n_sections, f_vector_t(N_t + 1, h));

        f_vector_2d_t voltageVec(n_sections, f_vector_t(N_t + 1, V));

        f_vector_2d_t dphiVec(n_sections, f_vector_t(N_t + 1, dphi));

        Context::GP = new GeneralParameters(N_t, CVec, alphaVec,
                                            alpha_order, momentumVec,
                                            GeneralParameters::particle_t::proton);

        Context::Beam = new Beams(N_p, N_b);

        Context::RfP = new RfParameters(n_sections, hVec, voltageVec, dphiVec);

        longitudinal_bigaussian(tau_0 / 4, 0, -1, false);

        // Context::Slice = new Slices(N_slices);
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



TEST_F(testTrackerPeriodicity, set_periodicity1)
{
    auto Beam = Context::Beam;
    // auto epsilon = 1e-8;
    auto params = std::string(TEST_FILES "/Tracker/periodicity/set_periodicity1/");

    auto long_tracker = new RingAndRfSection(Context::RfP, RingAndRfSection::simple,
            NULL, NULL, true, 0.0);

    auto mean = mymath::mean(Beam->dt.data(), Beam->dt.size());
    Context::GP->t_rev[Context::RfP->counter + 1] = mean;

    long_tracker->set_periodicity();

    f_vector_t v;
    util::read_vector_from_file(v, params + "indices_right_outside.txt");
    auto res = long_tracker->indices_right_outside;
    ASSERT_EQ(v.size(), res.size());
    for (uint i = 0; i < v.size(); ++i) {
        int ref = v[i];
        int real = res[i];
        ASSERT_EQ(ref, real)
                << "Testing of Beam->indices_right_outside failed on i " << i
                << std::endl;
    }

    util::read_vector_from_file(v, params + "indices_inside_frame.txt");
    res = long_tracker->indices_inside_frame;
    ASSERT_EQ(v.size(), res.size());
    for (uint i = 0; i < v.size(); ++i) {
        int ref = v[i];
        int real = res[i];
        ASSERT_EQ(ref, real)
                << "Testing of Beam->indices_inside_frame failed on i " << i
                << std::endl;
    }

    delete long_tracker;
}


TEST_F(testTrackerPeriodicity, track1)
{
    auto Beam = Context::Beam;
    auto epsilon = 1e-8;
    auto params = std::string(TEST_FILES "/Tracker/periodicity/track1/");

    auto mean = mymath::mean(Beam->dt.data(), Beam->dt.size());
    Context::GP->t_rev[Context::RfP->counter + 1] = mean;

    auto long_tracker = new RingAndRfSection(Context::RfP, RingAndRfSection::simple,
            NULL, NULL, true, 0.0);


    long_tracker->track();

    f_vector_t v;
    util::read_vector_from_file(v, params + "dE.txt");
    ASSERT_EQ(v.size(), Context::Beam->dE.size());
    for (uint i = 0; i < v.size(); ++i) {
        auto ref = v[i];
        auto real = Context::Beam->dE[i];
        ASSERT_NEAR(ref, real, epsilon * std::max(std::abs(ref), std::abs(real)))
                << "Testing of dE failed on i " << i << '\n';

    }

    util::read_vector_from_file(v, params + "dt.txt");
    ASSERT_EQ(v.size(), Context::Beam->dt.size());
    for (uint i = 0; i < v.size(); ++i) {
        auto ref = v[i];
        auto real = Context::Beam->dt[i];
        ASSERT_NEAR(ref, real, epsilon * std::max(std::abs(ref), std::abs(real)))
                << "Testing of dt failed on i " << i << '\n';

    }

    delete long_tracker;
}


TEST_F(testTrackerPeriodicity, track2)
{
    auto Beam = Context::Beam;
    auto GP = Context::GP;
    auto epsilon = 1e-8;
    auto params = std::string(TEST_FILES "/Tracker/periodicity/track2/");

    auto mean = mymath::mean(Beam->dt.data(), Beam->dt.size());
    int size = GP->t_rev.size();
    GP->t_rev = f_vector_t(size, mean);

    auto long_tracker = new RingAndRfSection(Context::RfP, RingAndRfSection::simple,
            NULL, NULL, true, 0.0);


    for (int i = 0; i < N_t; i++) long_tracker->track();

    f_vector_t v;
    util::read_vector_from_file(v, params + "dE.txt");
    ASSERT_EQ(v.size(), Beam->dE.size());
    for (uint i = 0; i < v.size(); ++i) {
        auto ref = v[i];
        auto real = Beam->dE[i];
        ASSERT_NEAR(ref, real, epsilon * std::max(std::abs(ref), std::abs(real)))
                << "Testing of dE failed on i " << i << '\n';

    }

    util::read_vector_from_file(v, params + "dt.txt");
    ASSERT_EQ(v.size(), Beam->dt.size());
    for (uint i = 0; i < v.size(); ++i) {
        auto ref = v[i];
        auto real = Beam->dt[i];
        ASSERT_NEAR(ref, real, epsilon * std::max(std::abs(ref), std::abs(real)))
                << "Testing of dt failed on i " << i << '\n';

    }

    delete long_tracker;
}


int main(int ac, char *av[])
{
    ::testing::InitGoogleTest(&ac, av);
    return RUN_ALL_TESTS();
}

#include <blond/blond.h>

#include <testing_utilities.h>
#include <gtest/gtest.h>

using namespace std;

class testTracker : public ::testing::Test {

protected:
    const long long N_b = 1e9;  // Intensity
    const double tau_0 = 0.4e-9; // Initial bunch length, 4 sigma [s]
    // Machine and RF parameters
    const double C = 26658.883;       // Machine circumference [m]
    const double p_i = 450e9;     // Synchronous momentum [eV/c]
    const double p_f = 460.005e9;     // Synchronous momentum, final
    const long long h = 35640;       // Harmonic number
    const double V = 6e6;             // RF voltage [V]
    const double dphi = 0;            // Phase modulation/offset
    const double gamma_t = 55.759505; // Transition gamma
    const double alpha =
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

        f_vector_2d_t alphaVec(n_sections, f_vector_t(alpha_order, alpha));

        f_vector_t CVec(n_sections, C);

        f_vector_2d_t hVec(n_sections, f_vector_t(N_t + 1, h));

        f_vector_2d_t voltageVec(n_sections, f_vector_t(N_t + 1, V));

        f_vector_2d_t dphiVec(n_sections, f_vector_t(N_t + 1, dphi));

        Context::GP = new GeneralParameters(N_t, CVec, alphaVec,
                                             momentumVec,
                                            GeneralParameters::particle_t::proton);

        // auto Beam = Context::Beam = new Beams(N_p, N_b);

        auto GP = Context::GP;
        auto Beam = Context::Beam = new Beams(GP, N_p, N_b);

        auto RfP = Context::RfP = new RfParameters(GP, n_sections, hVec,
                voltageVec, dphiVec);


        longitudinal_bigaussian(GP, RfP, Beam, tau_0 / 4, 0, -1, false);

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



TEST_F(testTracker, kick1)
{
    auto Beam = Context::Beam;
    auto epsilon = 1e-8;
    std::string params = TEST_FILES "/Tracker/kick1/";

    auto long_tracker = new RingAndRfSection();

    for (int i = 0; i < 10; i++) long_tracker->kick(Beam->dt, Beam->dE, i);

    f_vector_t v;
    util::read_vector_from_file(v, params + "dE.txt");
    ASSERT_NEAR_LOOP(v, Beam->dE, "dE", epsilon);

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
    ASSERT_NEAR_LOOP(v, Beam->dt, "dt", epsilon);

    delete long_tracker;
}


TEST_F(testTracker, track1)
{

    auto epsilon = 1e-8;
    auto Beam = Context::Beam;
    std::string params = TEST_FILES "/Tracker/track1/";

    auto long_tracker = new RingAndRfSection();
    long_tracker->track();

    f_vector_t v;
    util::read_vector_from_file(v, params + "dE.txt");
    ASSERT_NEAR_LOOP(v, Beam->dE, "dE", epsilon);

    util::read_vector_from_file(v, params + "dt.txt");
    ASSERT_NEAR_LOOP(v, Beam->dt, "dt", epsilon);


    delete long_tracker;
}

TEST_F(testTracker, track2)
{

    auto epsilon = 1e-8;
    std::string params = TEST_FILES "/Tracker/track2/";
    auto Beam = Context::Beam;

    auto long_tracker = new RingAndRfSection();

    for (int i = 0; i < 10; i++) long_tracker->track();

    f_vector_t v;
    util::read_vector_from_file(v, params + "dE.txt");
    ASSERT_NEAR_LOOP(v, Beam->dE, "dE", epsilon);

    util::read_vector_from_file(v, params + "dt.txt");
    ASSERT_NEAR_LOOP(v, Beam->dt, "dt", epsilon);

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
    ASSERT_NEAR_LOOP(v, long_tracker->fRfVoltage, "rf_voltage", epsilon);

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
    ASSERT_NEAR_LOOP(v, long_tracker->fRfVoltage, "rf_voltage", epsilon);

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
    ASSERT_NEAR_LOOP(v, long_tracker->fRfVoltage, "rf_voltage", epsilon);

    delete long_tracker;
}


class testTracker2 : public ::testing::Test {

protected:
    const long long N_b = 1e9;  // Intensity
    const double tau_0 = 0.4e-9; // Initial bunch length, 4 sigma [s]
    // Machine and RF parameters
    const double C = 26658.883;       // Machine circumference [m]
    const double p_i = 460e9;     // Synchronous momentum [eV/c]
    const double p_f = 460.005e9;     // Synchronous momentum, final
    const long long h = 35640;       // Harmonic number
    const double V = 6e6;             // RF voltage [V]
    const double dphi = 0;            // Phase modulation/offset
    const double gamma_t = 55.759505; // Transition gamma
    const double alpha =
        1.0 / gamma_t / gamma_t; // First order mom. comp. factor
    // const int alpha_order = 3;
    const int n_sections = 1;
    // Tracking details

    const int N_t = 100; // Number of turns to track
    const int N_p = 100;  // Macro-particles
    const int N_slices = 10;

    virtual void SetUp()
    {
        omp_set_num_threads(1);
    }

    virtual void TearDown()
    {
    }
};


TEST_F(testTracker2, full_solver1)
{
    auto epsilon = 1e-8;
    string params = TEST_FILES "/Tracker/full_solver1/";

    int alpha_order = 3;

    f_vector_2d_t momentumVec(n_sections);
    for (auto &v : momentumVec)
        v = mymath::linspace(p_i, p_f, N_t + 1);

    f_vector_2d_t alphaVec({{alpha, alpha / gamma_t, 2 * alpha / gamma_t}});

    f_vector_t CVec(n_sections, C);

    f_vector_2d_t hVec(n_sections, f_vector_t(N_t + 1, h));

    f_vector_2d_t voltageVec(n_sections, f_vector_t(N_t + 1, V));

    f_vector_2d_t dphiVec(n_sections, f_vector_t(N_t + 1, dphi));

    auto GP = GeneralParameters(N_t, CVec, alphaVec,
                                 momentumVec,
                                GeneralParameters::particle_t::proton);


    auto Beam = Beams(&GP, N_p, N_b);

    auto RfP = RfParameters(&GP, n_sections, hVec,
                            voltageVec, dphiVec);


    longitudinal_bigaussian(&GP, &RfP, &Beam, tau_0 / 4, 0, -1, false);
    auto long_tracker = RingAndRfSection(&RfP, &Beam, RingAndRfSection::full);
    f_vector_t v;

    for (int i = 0; i < 10; i++)
        long_tracker.track();

    util::read_vector_from_file(v, params + "dE.txt");
    ASSERT_NEAR_LOOP(v, Beam.dE, "dE", epsilon);

    util::read_vector_from_file(v, params + "dt.txt");
    ASSERT_NEAR_LOOP(v, Beam.dt, "dt", epsilon);

}


TEST_F(testTracker2, full_solver2)
{
    auto epsilon = 1e-8;
    string params = TEST_FILES "/Tracker/full_solver2/";
    int alpha_order = 3;

    f_vector_2d_t momentumVec(n_sections);
    for (auto &v : momentumVec)
        v = mymath::linspace(p_i, p_f, N_t + 1);

    f_vector_2d_t alphaVec({{alpha, alpha / gamma_t, 2 * alpha / gamma_t}});

    f_vector_t CVec(n_sections, C);

    f_vector_2d_t hVec(n_sections, f_vector_t(N_t + 1, h));

    f_vector_2d_t voltageVec(n_sections, f_vector_t(N_t + 1, V));

    f_vector_2d_t dphiVec(n_sections, f_vector_t(N_t + 1, dphi));

    auto GP = GeneralParameters(N_t, CVec, alphaVec,
                                 momentumVec,
                                GeneralParameters::particle_t::proton);


    auto Beam = Beams(&GP, N_p, N_b);

    auto RfP = RfParameters(&GP, n_sections, hVec,
                            voltageVec, dphiVec);


    longitudinal_bigaussian(&GP, &RfP, &Beam, tau_0 / 4, 0, -1, false);
    auto long_tracker = RingAndRfSection(&RfP, &Beam, RingAndRfSection::full);

    f_vector_t v;

    for (int i = 0; i < N_t; i++)
        long_tracker.track();

    util::read_vector_from_file(v, params + "dE.txt");
    ASSERT_NEAR_LOOP(v, Beam.dE, "dE", epsilon);

    util::read_vector_from_file(v, params + "dt.txt");
    ASSERT_NEAR_LOOP(v, Beam.dt, "dt", epsilon);

}

TEST_F(testTracker2, full_solver3)
{
    auto epsilon = 1e-8;
    string params = TEST_FILES "/Tracker/full_solver3/";
    int alpha_order = 1;

    f_vector_2d_t momentumVec(n_sections);
    for (auto &v : momentumVec)
        v = mymath::linspace(p_i, p_f, N_t + 1);

    f_vector_2d_t alphaVec({{alpha}});

    f_vector_t CVec(n_sections, C);

    f_vector_2d_t hVec(n_sections, f_vector_t(N_t + 1, h));

    f_vector_2d_t voltageVec(n_sections, f_vector_t(N_t + 1, V));

    f_vector_2d_t dphiVec(n_sections, f_vector_t(N_t + 1, dphi));

    auto GP = GeneralParameters(N_t, CVec, alphaVec,
                                 momentumVec,
                                GeneralParameters::particle_t::proton);


    auto Beam = Beams(&GP, N_p, N_b);

    auto RfP = RfParameters(&GP, n_sections, hVec,
                            voltageVec, dphiVec);


    longitudinal_bigaussian(&GP, &RfP, &Beam, tau_0 / 4, 0, -1, false);
    auto long_tracker = RingAndRfSection(&RfP, &Beam, RingAndRfSection::full);

    f_vector_t v;

    for (int i = 0; i < N_t; i++)
        long_tracker.track();

    util::read_vector_from_file(v, params + "dE.txt");
    ASSERT_NEAR_LOOP(v, Beam.dE, "dE", epsilon);

    util::read_vector_from_file(v, params + "dt.txt");
    ASSERT_NEAR_LOOP(v, Beam.dt, "dt", epsilon);

}

TEST_F(testTracker2, full_solver4)
{
    auto epsilon = 1e-8;
    string params = TEST_FILES "/Tracker/full_solver4/";
    int alpha_order = 2;

    f_vector_2d_t momentumVec(n_sections);
    for (auto &v : momentumVec)
        v = mymath::linspace(p_i, p_f, N_t + 1);

    f_vector_2d_t alphaVec({{alpha, alpha / gamma_t}});

    f_vector_t CVec(n_sections, C);

    f_vector_2d_t hVec(n_sections, f_vector_t(N_t + 1, h));

    f_vector_2d_t voltageVec(n_sections, f_vector_t(N_t + 1, V));

    f_vector_2d_t dphiVec(n_sections, f_vector_t(N_t + 1, dphi));

    auto GP = GeneralParameters(N_t, CVec, alphaVec,
                                 momentumVec,
                                GeneralParameters::particle_t::proton);


    auto Beam = Beams(&GP, N_p, N_b);

    auto RfP = RfParameters(&GP, n_sections, hVec,
                            voltageVec, dphiVec);


    longitudinal_bigaussian(&GP, &RfP, &Beam, tau_0 / 4, 0, -1, false);
    auto long_tracker = RingAndRfSection(&RfP, &Beam, RingAndRfSection::full);

    f_vector_t v;

    for (int i = 0; i < N_t; i++)
        long_tracker.track();

    util::read_vector_from_file(v, params + "dE.txt");
    ASSERT_NEAR_LOOP(v, Beam.dE, "dE", epsilon);

    util::read_vector_from_file(v, params + "dt.txt");
    ASSERT_NEAR_LOOP(v, Beam.dt, "dt", epsilon);

}



class testTrackerMultiRf : public ::testing::Test {

protected:
    const long long N_b = 1e9;  // Intensity
    const double tau_0 = 0.4e-9; // Initial bunch length, 4 sigma [s]
    // Machine and RF parameters
    const double C = 26658.883;       // Machine circumference [m]
    const double p_i = 450e9;     // Synchronous momentum [eV/c]
    const double p_f = 460.005e9;     // Synchronous momentum, final
    const long long h = 35640;       // Harmonic number
    const double V = 6e6;             // RF voltage [V]
    const double dphi = 0;            // Phase modulation/offset
    const double gamma_t = 55.759505; // Transition gamma
    const double alpha =
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

        f_vector_2d_t alphaVec(n_sections, f_vector_t(alpha_order, alpha));

        f_vector_t CVec(n_sections, C);

        f_vector_2d_t hVec(n_rf);
        for (int i = 0; i < n_rf; i++)
            hVec[i] = f_vector_t(N_t + 1, (i + 1) * h);

        f_vector_2d_t voltageVec(n_rf);

        for (int i = 0; i < n_rf; i++)
            voltageVec[i] = f_vector_t(N_t + 1, (i + 1) * V);

        f_vector_2d_t dphiVec(n_rf, f_vector_t(N_t + 1, dphi));

        Context::GP = new GeneralParameters(N_t, CVec, alphaVec,
                                             momentumVec,
                                            GeneralParameters::particle_t::proton);
        auto GP = Context::GP;

        Context::Beam = new Beams(GP, N_p, N_b);
        auto Beam = Context::Beam;

        Context::RfP = new RfParameters(GP, n_rf, hVec, voltageVec, dphiVec);
        auto RfP = Context::RfP;

        longitudinal_bigaussian(GP, RfP, Beam, tau_0 / 4, 0, -1, false);

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
    ASSERT_NEAR_LOOP(v, long_tracker->fRfVoltage, "rf_voltage", epsilon);

    delete long_tracker;
}


class testTrackerPeriodicity : public ::testing::Test {

protected:

    const long long N_b = 1e9;  // Intensity
    const double tau_0 = 1e-9; // Initial bunch length, 4 sigma [s]
    // Machine and RF parameters
    const double C = 26658.883;       // Machine circumference [m]
    const double p_i = 455e9;     // Synchronous momentum [eV/c]
    const double p_f = 460.005e9;     // Synchronous momentum, final
    const long long h = 35640;       // Harmonic number
    const double V = 10e6;             // RF voltage [V]
    const double dphi = 0;            // Phase modulation/offset
    const double gamma_t = 55.759505; // Transition gamma
    const double alpha =
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
                                             momentumVec,
                                            GeneralParameters::particle_t::proton);

        auto GP = Context::GP;

        auto Beam = Context::Beam = new Beams(GP, N_p, N_b);

        auto RfP = Context::RfP = new RfParameters(GP, n_sections, hVec,
                voltageVec, dphiVec);


        longitudinal_bigaussian(GP, RfP, Beam, tau_0 / 4, 0, -1, false);

        // Context::Slice = new Slices(RfP, Beam, N_slices);
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

    auto long_tracker = new RingAndRfSection(Context::RfP, Beam,
            RingAndRfSection::simple, NULL, NULL, true, 0.0);

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

    auto long_tracker = new RingAndRfSection(Context::RfP, Beam,
            RingAndRfSection::simple, NULL, NULL, true, 0.0);


    long_tracker->track();

    f_vector_t v;
    util::read_vector_from_file(v, params + "dE.txt");
    ASSERT_NEAR_LOOP(v, Beam->dE, "dE", epsilon);

    util::read_vector_from_file(v, params + "dt.txt");
    ASSERT_NEAR_LOOP(v, Beam->dt, "dt", epsilon);

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

    auto long_tracker = new RingAndRfSection(Context::RfP, Beam,
            RingAndRfSection::simple, NULL, NULL, true, 0.0);


    for (int i = 0; i < N_t; i++) long_tracker->track();

    f_vector_t v;
    util::read_vector_from_file(v, params + "dE.txt");
    ASSERT_NEAR_LOOP(v, Beam->dE, "dE", epsilon);

    util::read_vector_from_file(v, params + "dt.txt");
    ASSERT_NEAR_LOOP(v, Beam->dt, "dt", epsilon);

    delete long_tracker;
}


int main(int ac, char *av[])
{
    ::testing::InitGoogleTest(&ac, av);
    return RUN_ALL_TESTS();
}

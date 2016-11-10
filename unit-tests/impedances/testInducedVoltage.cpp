
#include <blond/beams/Beams.h>
#include <blond/beams/Distributions.h>
#include <blond/beams/Slices.h>
#include <blond/globals.h>
#include <blond/impedances/InducedVoltage.h>
#include <blond/input_parameters/GeneralParameters.h>
#include <blond/input_parameters/RfParameters.h>
#include <blond/math_functions.h>
#include <blond/trackers/Tracker.h>
#include <blond/utilities.h>
#include <gtest/gtest.h>
#include <stdio.h>

using namespace std;

class testInducedVoltage : public ::testing::Test {

protected:
    // Simulation parameters
    // --------------------------------------------------------
    // Bunch parameters
    const long long N_b = 1e10; // Intensity
    const double tau_0 = 2e-9; // Initial bunch length, 4 sigma [s]
    // const particle_type particle = proton;
    // Machine and RF parameters
    const double C = 6911.56;   // Machine circumference [m]
    const double p_i = 25.92e9; // Synchronous momentum [eV/c]
    // const double p_f = 460.005e9;                  // Synchronous momentum,
    // final
    const double h = 4620;                     // Harmonic number
    const double V = 0.9e6;                        // RF voltage [V]
    const double dphi = 0;                         // Phase modulation/offset
    const double gamma_t = 1 / std::sqrt(0.00192); // Transition gamma
    const double alpha =
        1.0 / gamma_t / gamma_t; // First order mom. comp. factor
    const int alpha_order = 1;
    const int n_sections = 1;
    // Tracking details

    long N_t = 1000;       // Number of turns to track
    long N_p = 5000; // Macro-particles

    long N_slices = 1 << 8; // = (2^8)
    const string datafiles = DEMO_FILES "/TC5_Wake_impedance/";
    Resonators *resonator;

    virtual void SetUp()
    {

        omp_set_num_threads(1);

        f_vector_2d_t momentumVec(n_sections, f_vector_t(N_t + 1, p_i));

        f_vector_2d_t alphaVec(n_sections, f_vector_t(alpha_order + 1, alpha));

        f_vector_t CVec(n_sections, C);

        f_vector_2d_t hVec(n_sections, f_vector_t(N_t + 1, h));

        f_vector_2d_t voltageVec(n_sections, f_vector_t(N_t + 1, V));

        f_vector_2d_t dphiVec(n_sections, f_vector_t(N_t + 1, dphi));

        Context::GP = new GeneralParameters(N_t, CVec, alphaVec,
                                            alpha_order, momentumVec,
                                            GeneralParameters::particle_t::proton);
        // Context::Beam = new Beams(N_p, N_b);

        auto GP = Context::GP;
        auto Beam = Context::Beam = new Beams(GP, N_p, N_b);

        auto RfP = Context::RfP = new RfParameters(GP, n_sections, hVec,
                voltageVec, dphiVec);



        longitudinal_bigaussian(GP, RfP, Beam, tau_0 / 4, 0, -1, false);

        Context::Slice = new Slices(RfP, Beam, N_slices, 0, 0, 2 * constant::pi,
                                    Slices::cuts_unit_t::rad);

        f_vector_t v;
        util::read_vector_from_file(v, datafiles + "TC5_new_HQ_table.dat");
        assert(v.size() % 3 == 0);

        f_vector_t R_shunt, f_res, Q_factor;

        R_shunt.reserve(v.size() / 3);
        f_res.reserve(v.size() / 3);
        Q_factor.reserve(v.size() / 3);

        for (uint i = 0; i < v.size(); i += 3) {
            f_res.push_back(v[i] * 1e9);
            Q_factor.push_back(v[i + 1]);
            R_shunt.push_back(v[i + 2] * 1e6);
        }

        resonator = new Resonators(R_shunt, f_res, Q_factor);
    }

    virtual void TearDown()
    {
        // Code here will be called immediately after each test
        // (right before the destructor).
        delete Context::GP;
        delete Context::Beam;
        delete Context::RfP;
        delete Context::Slice;
        delete resonator;
    }
};


class testTotalInducedVoltage : public testInducedVoltage {};


TEST_F(testInducedVoltage, constructor1)
{
    double epsilon = 1e-8;
    auto slices = Context::Slice;
    auto indVoltTime = new InducedVoltageTime(slices, {resonator});

    auto params = std::string(TEST_FILES "/Impedances/") +
                  "InducedVoltage/InducedVoltageTime/constructor1/";

    f_vector_t v;
    util::read_vector_from_file(v, params + "time_array.txt");

    ASSERT_EQ(v.size(), indVoltTime->fTimeArray.size());

    for (uint i = 0; i < v.size(); ++i) {
        double ref = v[i];
        double real = indVoltTime->fTimeArray[i];
        ASSERT_NEAR(ref, real, epsilon * std::max(std::abs(ref), std::abs(real)))
                << "Testing of indVoltTime->fTimeArray failed on i " << i
                << std::endl;
    }

    util::read_vector_from_file(v, params + "total_wake.txt");
    ASSERT_EQ(v.size(), indVoltTime->fTotalWake.size());

    for (uint i = 0; i < v.size(); ++i) {
        double ref = v[i];
        double real = indVoltTime->fTotalWake[i];
        ASSERT_NEAR(ref, real, epsilon * std::max(std::abs(ref), std::abs(real)))
                << "Testing of indVoltTime->fTotalWake failed on i " << i
                << std::endl;
    }

    util::read_vector_from_file(v, params + "cut.txt");
    ASSERT_EQ(v.size(), 1);
    ASSERT_EQ(v[0], indVoltTime->fCut) << "Testing of fCut failed\n";

    util::read_vector_from_file(v, params + "fshape.txt");
    ASSERT_EQ(v.size(), 1);
    ASSERT_EQ(v[0], indVoltTime->fShape) << "Testing of fShape failed\n";

    delete indVoltTime;
}


TEST_F(testInducedVoltage, reprocess1)
{
    auto slices = Context::Slice;
    auto beam = Context::Beam;
    double epsilon = 1e-8;

    auto indVoltTime = new InducedVoltageTime(slices, {resonator});
    slices->track();

    for (int i = 0; i < slices->n_slices; i++)
        slices->bin_centers[i] *= 1.1;

    indVoltTime->reprocess(slices);

    std::string params = std::string(TEST_FILES "/Impedances/") +
                         "InducedVoltage/InducedVoltageTime/reprocess1/";

    f_vector_t v;
    util::read_vector_from_file(v, params + "time_array.txt");

    ASSERT_EQ(v.size(), indVoltTime->fTimeArray.size());


    for (uint i = 0; i < v.size(); ++i) {
        double ref = v[i];
        double real = indVoltTime->fTimeArray[i];
        ASSERT_NEAR(ref, real, epsilon * std::max(std::abs(ref), std::abs(real)))
                << "Testing of indVoltTime->fTimeArray failed on i " << i
                << std::endl;
    }

    util::read_vector_from_file(v, params + "total_wake.txt");
    ASSERT_EQ(v.size(), indVoltTime->fTotalWake.size());

    epsilon = 1e-8;
    for (uint i = 0; i < v.size(); ++i) {
        double ref = v[i];
        double real = indVoltTime->fTotalWake[i];
        ASSERT_NEAR(ref, real, epsilon * std::max(std::abs(ref), std::abs(real)))
                << "Testing of indVoltTime->fTotalWake failed on i " << i
                << std::endl;
    }

    util::read_vector_from_file(v, params + "cut.txt");

    for (uint i = 0; i < v.size(); ++i) {
        uint ref = v[i];
        uint real = indVoltTime->fCut;
        ASSERT_EQ(ref, real) << "Testing of indVoltTime->fCut failed on i " << i
                             << std::endl;
    }

    util::read_vector_from_file(v, params + "fshape.txt");

    for (uint i = 0; i < v.size(); ++i) {
        uint ref = v[i];
        uint real = indVoltTime->fShape;
        ASSERT_EQ(ref, real) << "Testing of fShape failed on i " << i
                             << std::endl;
    }

    delete indVoltTime;
}

TEST_F(testInducedVoltage, generation1)
{
    auto slices = Context::Slice;
    auto beam = Context::Beam;

    slices->track();
    auto epsilon = 1e-8;

    auto indVoltTime = new InducedVoltageTime(slices, {resonator});
    indVoltTime->induced_voltage_generation(beam);
    auto res = indVoltTime->fInducedVoltage;
    std::string params = std::string(TEST_FILES "/Impedances/") +
                         "InducedVoltage/InducedVoltageTime/generation1/";

    f_vector_t v;
    util::read_vector_from_file(v, params + "induced_voltage.txt");

    ASSERT_EQ(v.size(), res.size());

    double max = *max_element(res.begin(), res.end(), [](double i, double j) {
        return std::abs(i) < std::abs(j);
    });
    max = std::abs(max);
    for (uint i = 0; i < v.size(); ++i) {
        double ref = v[i];
        double real = res[i];

        ASSERT_NEAR(ref, real, epsilon * max)
                << "Testing of indVoltTime->fInducedVoltage failed on i " << i
                << std::endl;
    }

    delete indVoltTime;
}


TEST_F(testInducedVoltage, generation2)
{
    auto slices = Context::Slice;
    auto beam = Context::Beam;

    slices->track();
    auto epsilon = 1e-7;

    auto indVoltTime = new InducedVoltageTime(slices, {resonator});
    f_vector_t res = indVoltTime->induced_voltage_generation(beam, 100);

    std::string params = std::string(TEST_FILES "/Impedances/") +
                         "InducedVoltage/InducedVoltageTime/generation2/";

    f_vector_t v;
    util::read_vector_from_file(v, params + "induced_voltage.txt");

    ASSERT_EQ(v.size(), res.size());

    double max = *max_element(res.begin(), res.end(), [](double i, double j) {
        return std::abs(i) < std::abs(j);
    });
    max = std::abs(max);

    for (uint i = 0; i < v.size(); ++i) {
        double ref = v[i];
        double real = res[i];
        ASSERT_NEAR(ref, real, epsilon * max)
                << "Testing of indVoltTime->fInducedVoltage failed on i " << i
                << std::endl;
    }

    delete indVoltTime;
}


TEST_F(testInducedVoltage, convolution1)
{
    auto slices = Context::Slice;
    auto beam = Context::Beam;

    slices->track();
    auto epsilon = 1e-8;

    auto indVoltTime = new InducedVoltageTime(slices, {resonator},
            InducedVoltageTime::time_or_freq::time_domain);
    indVoltTime->induced_voltage_generation(beam);
    auto res = indVoltTime->fInducedVoltage;

    std::string params = std::string(TEST_FILES "/Impedances/") +
                         "InducedVoltage/InducedVoltageTime/convolution1/";

    f_vector_t v;
    util::read_vector_from_file(v, params + "induced_voltage.txt");

    ASSERT_EQ(v.size(), res.size());

    double max = *max_element(res.begin(), res.end(), [](double i, double j) {
        return std::abs(i) < std::abs(j);
    });
    max = std::abs(max);

    for (uint i = 0; i < v.size(); ++i) {
        double ref = v[i];
        double real = res[i];
        ASSERT_NEAR(ref, real, epsilon * max)
                << "Testing of indVoltTime->fInducedVoltage failed on i " << i
                << std::endl;
    }
    delete indVoltTime;
}



TEST_F(testInducedVoltage, track1)
{
    auto slices = Context::Slice;
    auto beam = Context::Beam;

    std::string params = std::string(TEST_FILES "/Impedances/") +
                         "InducedVoltage/InducedVoltageTime/track1/";
    auto epsilon = 1e-8;

    slices->track();

    auto indVoltTime = new InducedVoltageTime(slices, {resonator});
    indVoltTime->track(beam);

    f_vector_t v;
    util::read_vector_from_file(v, params + "dE.txt");
    ASSERT_EQ(v.size(), beam->dE.size());

    for (uint i = 0; i < v.size(); ++i) {
        double ref = v[i];
        double real = beam->dE[i];
        ASSERT_NEAR(ref, real, epsilon * std::max(std::abs(ref), std::abs(real)))
                << "Testing of beam->dE failed on i " << i << std::endl;
    }

    delete indVoltTime;
}


TEST_F(testTotalInducedVoltage, sum1)
{
    auto slices = Context::Slice;
    auto beam = Context::Beam;
    auto epsilon = 1e-8;
    auto params = std::string(TEST_FILES "/Impedances/") +
                  "InducedVoltage/TotalInducedVoltage/sum1/";

    slices->track();

    auto indVoltTime = new InducedVoltageTime(slices, {resonator});
    auto totVol = new TotalInducedVoltage(beam, slices, {indVoltTime});

    f_vector_t res = totVol->induced_voltage_sum(beam, 200);

    f_vector_t v;
    util::read_vector_from_file(v, params + "extIndVolt.txt");
    ASSERT_EQ(v.size(), res.size());

    double max = *max_element(res.begin(), res.end(), [](double i, double j) {
        return std::abs(i) < std::abs(j);
    });
    max = std::abs(max);
    // warning checking only the first 100 elems
    for (uint i = 0; i < v.size(); ++i) {
        double ref = v[i];
        double real = res[i];

        ASSERT_NEAR(ref, real, epsilon * max)
                << "Testing of extIndVolt failed on i "
                << i << std::endl;
    }

    res.clear();

    res = totVol->fInducedVoltage;
    util::read_vector_from_file(v, params + "induced_voltage.txt");
    ASSERT_EQ(v.size(), res.size());

    max = *max_element(res.begin(), res.end(),
    [](double i, double j) { return std::abs(i) < std::abs(j); });
    max = std::abs(max);
    // warning checking only the first 100 elems
    for (uint i = 0; i < v.size(); ++i) {
        double ref = v[i];
        double real = res[i];

        ASSERT_NEAR(ref, real, epsilon * max)
                << "Testing of fInducedVoltage failed on i " << i << std::endl;
    }

    delete indVoltTime;
    delete totVol;
}


TEST_F(testTotalInducedVoltage, sum2)
{
    auto slices = Context::Slice;
    auto beam = Context::Beam;
    auto epsilon = 1e-8;
    auto params = std::string(TEST_FILES "/Impedances/") +
                  "InducedVoltage/TotalInducedVoltage/sum2/";

    slices->track();

    auto indVoltTime1 = new InducedVoltageTime(slices, {resonator});
    auto indVoltTime2 = new InducedVoltageTime(slices, {resonator});
    auto totVol = new TotalInducedVoltage(beam, slices, {indVoltTime1, indVoltTime2});

    f_vector_t res = totVol->induced_voltage_sum(beam, 200);

    f_vector_t v;
    util::read_vector_from_file(v, params + "extIndVolt.txt");
    ASSERT_EQ(v.size(), res.size());

    double max = *max_element(res.begin(), res.end(), [](double i, double j) {
        return std::abs(i) < std::abs(j);
    });
    max = std::abs(max);
    // warning checking only the first 100 elems
    for (uint i = 0; i < v.size(); ++i) {
        double ref = v[i];
        double real = res[i];

        ASSERT_NEAR(ref, real, epsilon * max)
                << "Testing of extIndVolt failed on i "
                << i << std::endl;
    }

    res.clear();

    res = totVol->fInducedVoltage;
    util::read_vector_from_file(v, params + "induced_voltage.txt");
    ASSERT_EQ(v.size(), res.size());

    max = *max_element(res.begin(), res.end(),
    [](double i, double j) { return std::abs(i) < std::abs(j); });
    max = std::abs(max);
    for (uint i = 0; i < v.size(); ++i) {
        double ref = v[i];
        double real = res[i];

        ASSERT_NEAR(ref, real, epsilon * max)
                << "Testing of fInducedVoltage failed on i " << i << std::endl;
    }

    delete indVoltTime1;
    delete indVoltTime2;
    delete totVol;
}





TEST_F(testTotalInducedVoltage, track1)
{
    auto slices = Context::Slice;
    auto beam = Context::Beam;

    auto epsilon = 1e-8;
    auto params = std::string(TEST_FILES "/Impedances/") +
                  "InducedVoltage/TotalInducedVoltage/track1/";


    auto indVoltTime = new InducedVoltageTime(slices, {resonator});
    auto totVol = new TotalInducedVoltage(beam, slices, {indVoltTime});

    for (int i = 0; i < 1000; i++) totVol->track(beam);

    f_vector_t v;
    util::read_vector_from_file(v, params + "dE.txt");
    ASSERT_EQ(v.size(), beam->dE.size());
    for (uint i = 0; i < v.size(); ++i) {
        double ref = v[i];
        double real = beam->dE[i];
        ASSERT_NEAR(ref, real, epsilon * std::max(std::abs(ref), std::abs(real)))
                << "Testing of beam->dE failed on i " << i << std::endl;
    }

    delete indVoltTime;
    delete totVol;
}


TEST_F(testTotalInducedVoltage, track2)
{
    auto slices = Context::Slice;
    auto beam = Context::Beam;

    auto epsilon = 1e-8;
    auto params = std::string(TEST_FILES "/Impedances/") +
                  "InducedVoltage/TotalInducedVoltage/track2/";


    auto indVoltTime = new InducedVoltageTime(slices, {resonator});
    auto indVoltTime2 = new InducedVoltageTime(slices, {resonator});
    auto totVol = new TotalInducedVoltage(beam, slices, {indVoltTime, indVoltTime2});

    for (int i = 0; i < 100; i++) totVol->track(beam);

    f_vector_t v;
    util::read_vector_from_file(v, params + "dE.txt");
    ASSERT_EQ(v.size(), beam->dE.size());
    for (uint i = 0; i < v.size(); ++i) {
        double ref = v[i];
        double real = beam->dE[i];
        ASSERT_NEAR(ref, real, epsilon * std::max(std::abs(ref), std::abs(real)))
                << "Testing of beam->dE failed on i " << i << std::endl;
    }

    delete indVoltTime;
    delete indVoltTime2;
    delete totVol;
}


int main(int ac, char *av[])
{
    ::testing::InitGoogleTest(&ac, av);
    return RUN_ALL_TESTS();
}

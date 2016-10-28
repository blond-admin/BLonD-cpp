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

Resonators *resonator;


class testTC5 : public ::testing::Test {

protected:

    const long long N_b = 1e10; // Intensity
    const double tau_0 = 2e-9;  // Initial bunch length, 4 sigma [s]
    const double C = 6911.56;   // Machine circumference [m]
    const double p_i = 25.92e9; // Synchronous momentum [eV/c]
    const double h = 4620;  // Harmonic number
    const double V = 0.9e6;     // RF voltage [V]
    const double dphi = 0;      // Phase modulation/offset
    const double gamma_t = 1 / std::sqrt(0.00192); // Transition gamma
    const double alpha =
        1.0 / gamma_t / gamma_t; // First order mom. comp. factor
    const int alpha_order = 1;
    const int n_sections = 1;
    uint N_t = 100; // Number of turns to track
    uint N_p = 10000; // Macro-particles
    uint N_slices = 1 << 8; // = (2^8)
    const std::string datafiles = DEMO_FILES "/TC5_Wake_impedance/";


    virtual void SetUp()
    {

        omp_set_num_threads(4);

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


        longitudinal_bigaussian(tau_0 / 4, 0, -1, false);

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


// TODO: Why so bad precision?
TEST_F(testTC5, timeTrack)
{
    auto Beam = Context::Beam;
    auto Slice = Context::Slice;
    auto epsilon = 1e-3;
    auto params = std::string(TEST_FILES "/TC5_final/time/");

    auto long_tracker = new RingAndRfSection();
    std::vector<Intensity *> wakeSourceList({resonator});
    auto indVoltTime = new InducedVoltageTime(Slice, wakeSourceList);
    std::vector<InducedVoltage *> indVoltList({indVoltTime});
    auto totVol = new TotalInducedVoltage(Beam, Slice, indVoltList);

    for (uint i = 0; i < N_t; ++i) {
        totVol->track(Beam);
        long_tracker->track();
        Slice->track();
    }

    f_vector_t v;
    util::read_vector_from_file(v, params + "dE.txt");
    ASSERT_EQ(v.size(), Beam->dE.size());
    for (uint i = 0; i < v.size(); ++i) {
        double ref = v[i];
        double real = Beam->dE[i];
        ASSERT_NEAR(ref, real, epsilon * std::max(std::abs(ref), std::abs(real)))
                << "Testing of Beam->dE failed on i " << i << std::endl;
    }

    util::read_vector_from_file(v, params + "dt.txt");
    ASSERT_EQ(v.size(), Beam->dt.size());
    for (uint i = 0; i < v.size(); ++i) {
        double ref = v[i];
        double real = Beam->dt[i];
        ASSERT_NEAR(ref, real, epsilon * std::max(std::abs(ref), std::abs(real)))
                << "Testing of Beam->dt failed on i " << i << std::endl;
    }

    util::read_vector_from_file(v, params + "n_macroparticles.txt");
    ASSERT_EQ(v.size(), Slice->n_macroparticles.size());
    for (uint i = 0; i < v.size(); ++i) {
        int ref = v[i];
        int real = Slice->n_macroparticles[i];
        ASSERT_EQ(ref, real)
                << "Testing of Slice->n_macroparticles failed on i " << i
                << std::endl;
    }

    delete indVoltTime;
    delete totVol;
    delete long_tracker;
}


// TODO: Why so bad precision?
TEST_F(testTC5, freqTrack)
{
    auto Beam = Context::Beam;
    auto Slice = Context::Slice;
    auto epsilon = 1e-3;
    auto params = std::string(TEST_FILES "/TC5_final/freq/");

    auto long_tracker = new RingAndRfSection();
    std::vector<Intensity *> ImpSourceList({resonator});
    auto indVoltFreq = new InducedVoltageFreq(Slice, ImpSourceList, 1e5);
    std::vector<InducedVoltage *> indVoltList({indVoltFreq});
    auto totVol = new TotalInducedVoltage(Beam, Slice, indVoltList);

    for (uint i = 0; i < N_t; ++i) {
        totVol->track(Beam);
        long_tracker->track();
        Slice->track();
    }
    // util::dump(Beam->dE, "dE\n", 20);

    f_vector_t v;
    util::read_vector_from_file(v, params + "dE.txt");
    ASSERT_EQ(v.size(), Beam->dE.size());
    for (uint i = 0; i < v.size(); ++i) {
        double ref = v[i];
        double real = Beam->dE[i];
        ASSERT_NEAR(ref, real, epsilon * std::max(std::abs(ref), std::abs(real)))
                << "Testing of Beam->dE failed on i " << i << std::endl;
    }

    util::read_vector_from_file(v, params + "dt.txt");
    ASSERT_EQ(v.size(), Beam->dt.size());
    for (uint i = 0; i < v.size(); ++i) {
        double ref = v[i];
        double real = Beam->dt[i];
        ASSERT_NEAR(ref, real, epsilon * std::max(std::abs(ref), std::abs(real)))
                << "Testing of Beam->dt failed on i " << i << std::endl;
    }

    util::read_vector_from_file(v, params + "n_macroparticles.txt");
    ASSERT_EQ(v.size(), Slice->n_macroparticles.size());
    for (uint i = 0; i < v.size(); ++i) {
        int ref = v[i];
        int real = Slice->n_macroparticles[i];
        ASSERT_EQ(ref, real)
                << "Testing of Slice->n_macroparticles failed on i " << i
                << std::endl;
    }

    delete indVoltFreq;
    delete totVol;
    delete long_tracker;
}

int main(int ac, char *av[])
{
    ::testing::InitGoogleTest(&ac, av);
    return RUN_ALL_TESTS();
}

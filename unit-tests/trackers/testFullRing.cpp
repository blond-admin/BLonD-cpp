#include <iostream>

#include <blond/beams/Distributions.h>
#include <blond/constants.h>
#include <blond/input_parameters/GeneralParameters.h>
#include <blond/math_functions.h>
#include <blond/trackers/Tracker.h>
#include <blond/utilities.h>
#include <gtest/gtest.h>


class testFullRing : public ::testing::Test {

protected:
    // Bunch parameters
    const uint N_b = 0; // Intensity

    // Machine and RF parameters
    const double radius = 25;
    const double C = 2 * constant::pi * radius; // Machine circumference [m]
    const double p_i = 310891054.809;           // Synchronous momentum [eV/c]
    const uint h = 1;                          // Harmonic number
    const double V = 8000;                      // RF voltage [V]
    const double dphi = -constant::pi;          // Phase modulation/offset
    const double gamma_t = 4.076750841;         // Transition gamma
    const double alpha =
        1.0 / gamma_t / gamma_t; // First order mom. comp. factor
    const uint alpha_order = 1;
    const uint n_sections = 1;
    // Tracking details

    uint N_t = 1000;  // Number of turns to track
    uint N_p = 1000; // Macro-particles

    uint N_slices = 100; // = (2^8)

    RfParameters *RfP1, *RfP2;
    RingAndRfSection *long_tracker1, *long_tracker2;


    virtual void SetUp()
    {
        omp_set_num_threads(1);

        f_vector_2d_t momentumVec(n_sections, f_vector_t(N_t + 1, p_i));

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

        Context::RfP = new RfParameters(GP, n_sections, hVec,
                                        voltageVec, dphiVec);

        RfP1 = new RfParameters(GP, n_sections, hVec, voltageVec, dphiVec);
        RfP2 = new RfParameters(GP, n_sections, hVec, voltageVec, dphiVec);
        long_tracker1 = new RingAndRfSection(RfP1);
        long_tracker2 = new RingAndRfSection(RfP2);
        // long_tracker = new RingAndRfSection(RfP);

        // Context::Slice = new Slices(RfP, Beam, N_slices, 0, -constant::pi, constant::pi,
        //                             cuts_unit_type::rad);
    }

    virtual void TearDown()
    {
        // Code here will be called immediately after each test
        // (right before the destructor).
        delete Context::GP;
        delete Context::Beam;
        delete Context::RfP;
        delete RfP1;
        delete RfP2;
        delete long_tracker1;
        delete long_tracker2;
    }
};

TEST_F(testFullRing, constructor1)
{
    auto Beam = Context::Beam;
    auto GP = Context::GP;
    auto RfP = Context::RfP;
    auto epsilon = 1e-8;
    auto params = std::string(TEST_FILES "/FullRing/constructor1/");

    longitudinal_bigaussian(GP, RfP, Beam, 200e-9, 1e6, -1, false);

    auto long_tracker = new RingAndRfSection(RfP);
    auto fullRing = new FullRingAndRf({long_tracker, long_tracker});

    f_vector_t v;
    util::read_vector_from_file(v, params + "ring_radius.txt");

    ASSERT_EQ(v.size(), 1);
    for (uint i = 0; i < v.size(); ++i) {
        auto ref = v[i];
        auto real = fullRing->fRingRadius;
        ASSERT_NEAR(ref, real, epsilon * std::max(std::abs(ref), std::abs(real)))
                << "Testing of ring_radius failed on i " << i << std::endl;
    }

    delete long_tracker;
    delete fullRing;
}

TEST_F(testFullRing, track1)
{
    auto Beam = Context::Beam;
    auto GP = Context::GP;
    auto RfP = Context::RfP;

    auto epsilon = 1e-8;

    auto params = std::string(TEST_FILES "/FullRing/track1/");

    longitudinal_bigaussian(GP, RfP, Beam, 200e-9, 1e6, -1, false);

    auto long_tracker = new RingAndRfSection(RfP);
    auto fullRing = new FullRingAndRf({long_tracker, long_tracker});

    for (uint i = 0; i < 100; ++i) fullRing->track();

    f_vector_t v;
    util::read_vector_from_file(v, params + "dE.txt");
    ASSERT_EQ(v.size(), Beam->dE.size());
    for (uint i = 0; i < v.size(); ++i) {
        auto ref = v[i];
        auto real = Beam->dE[i];
        ASSERT_NEAR(ref, real, epsilon * std::max(std::abs(ref), std::abs(real)))
                << "Testing of dE failed on i " << i << std::endl;
    }

    util::read_vector_from_file(v, params + "dt.txt");
    ASSERT_EQ(v.size(), Beam->dt.size());
    for (uint i = 0; i < v.size(); ++i) {
        auto ref = v[i];
        auto real = Beam->dt[i];
        ASSERT_NEAR(ref, real, epsilon * std::max(std::abs(ref), std::abs(real)))
                << "Testing of dt failed on i " << i << std::endl;
    }

    delete long_tracker;
    delete fullRing;
}

TEST_F(testFullRing, track2)
{
    auto Beam = Context::Beam;
    auto GP = Context::GP;
    auto RfP = Context::RfP;

    auto epsilon = 1e-8;
    auto params = std::string(TEST_FILES "/FullRing/track2/");

    longitudinal_bigaussian(GP, RfP, Beam, 200e-9, 1e6, -1, false);

    auto fullRing = new FullRingAndRf({long_tracker1, long_tracker2});

    for (uint i = 0; i < N_t; ++i) fullRing->track();

    f_vector_t v;
    util::read_vector_from_file(v, params + "dE.txt");
    ASSERT_EQ(v.size(), Beam->dE.size());
    for (uint i = 0; i < v.size(); ++i) {
        auto ref = v[i];
        auto real = Beam->dE[i];
        ASSERT_NEAR(ref, real, epsilon * std::max(std::abs(ref), std::abs(real)))
                << "Testing of dE failed on i " << i << std::endl;
    }

    util::read_vector_from_file(v, params + "dt.txt");
    ASSERT_EQ(v.size(), Beam->dt.size());
    for (uint i = 0; i < v.size(); ++i) {
        auto ref = v[i];
        auto real = Beam->dt[i];
        ASSERT_NEAR(ref, real, epsilon * std::max(std::abs(ref), std::abs(real)))
                << "Testing of dt failed on i " << i << std::endl;
    }

    delete fullRing;
}



TEST_F(testFullRing, potential_well_generation1)
{
    auto Beam = Context::Beam;
    auto GP = Context::GP;
    auto RfP = Context::RfP;

    auto params = std::string(TEST_FILES "/FullRing/potential_well_generation1/");
    auto epsilon = 1e-8;

    auto long_tracker = new RingAndRfSection(RfP);
    auto fullRing = new FullRingAndRf({long_tracker});

    fullRing->potential_well_generation(0, 1000);

    f_vector_t v;
    util::read_vector_from_file(v, params + "potential_well.txt");
    ASSERT_EQ(v.size(), fullRing->fPotentialWell.size());
    for (uint i = 0; i < v.size(); ++i) {
        auto ref = v[i];
        auto real = fullRing->fPotentialWell[i];
        ASSERT_NEAR(ref, real, epsilon * std::max(std::abs(ref), std::abs(real)))
                << "Testing of fPotentialWell failed on i " << i << std::endl;
    }

    util::read_vector_from_file(v, params + "potential_well_coordinates.txt");
    ASSERT_EQ(v.size(), fullRing->fPotentialWellCoordinates.size());
    for (uint i = 0; i < v.size(); ++i) {
        auto ref = v[i];
        auto real = fullRing->fPotentialWellCoordinates[i];
        ASSERT_NEAR(ref, real, epsilon * std::max(std::abs(ref), std::abs(real)))
                << "Testing of fPotentialWellCoordinates failed on i " << i << std::endl;
    }

    delete long_tracker;
    delete fullRing;
}



TEST_F(testFullRing, potential_well_generation2)
{
    auto Beam = Context::Beam;
    auto GP = Context::GP;
    auto RfP = Context::RfP;

    auto params = std::string(TEST_FILES "/FullRing/potential_well_generation2/");
    auto epsilon = 1e-8;

    auto long_tracker = new RingAndRfSection(RfP);
    auto fullRing = new FullRingAndRf({long_tracker, long_tracker});

    fullRing->potential_well_generation(10, 1000, 1);

    f_vector_t v;
    util::read_vector_from_file(v, params + "potential_well.txt");
    ASSERT_EQ(v.size(), fullRing->fPotentialWell.size());
    for (uint i = 0; i < v.size(); ++i) {
        auto ref = v[i];
        auto real = fullRing->fPotentialWell[i];
        ASSERT_NEAR(ref, real, epsilon * std::max(std::abs(ref), std::abs(real)))
                << "Testing of fPotentialWell failed on i " << i << std::endl;
    }

    delete long_tracker;
    delete fullRing;
}



int main(int ac, char *av[])
{
    ::testing::InitGoogleTest(&ac, av);
    return RUN_ALL_TESTS();
}

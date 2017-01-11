#include <blond/blond.h>
#include <testing_utilities.h>
#include <gtest/gtest.h>
using namespace std;

class testTWC : public ::testing::Test {

protected:
    // Simulation parameters
    // --------------------------------------------------------
    // Bunch parameters
    const long long N_b = 1e10; // Intensity
    const double tau_0 = 2e-9;  // Initial bunch length, 4 sigma [s]
    // Machine and RF parameters
    const double C = 6911.56;   // Machine circumference [m]
    const long long h = 4620;                     // Harmonic number
    const double V = 0.9e6;                        // RF voltage [V]
    const double dphi = 0;                         // Phase modulation/offset
    const double gamma_t = 1 / std::sqrt(0.00192); // Transition gamma
    // Derived parameters
    const double sync_momentum = 25.92e9; // # [eV / c]
    const double alpha = 1.0 / gamma_t / gamma_t;  // First order mom. comp. factor
    const int alpha_order = 1;
    const int n_sections = 1;
    // Tracking details

    int N_t = 10;       // Number of turns to track
    int N_p = 1000; // Macro-particles

    int N_slices = 10; // = (2^8)
    virtual void SetUp()
    {

        omp_set_num_threads(1);

        f_vector_2d_t momentumVec(n_sections, f_vector_t(N_t + 1, sync_momentum));

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


        // RingAndRfSection *long_tracker = new RingAndRfSection();

        longitudinal_bigaussian(GP, RfP, Beam, tau_0 / 4, 0, -1, false);

        Context::Slice = new Slices(RfP, Beam, N_slices, 0, 0,
                                    2 * constant::pi, Slices::cuts_unit_t::rad);


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


TEST_F(testTWC, wake_calc1)
{
    auto Slice = Context::Slice;
    std::string params = TEST_FILES "/Impedances/Intensity/TWC/wake_calc1/";
    f_vector_t v;
    double epsilon = 1e-8;

    f_vector_t RS, frequency, a_factor;
    f_vector_t random;
    util::read_vector_from_file(random, TEST_FILES "/normal_distribution.dat");

    RS = slice(random, 0, 10) * 1e6;
    apply_f_in_place(RS, std::abs<double>);
    frequency = slice(random, 10, 20) * 1e9;
    apply_f_in_place(frequency, std::abs<double>);
    a_factor = slice(random, 20, 30);
    apply_f_in_place(a_factor, std::abs<double>);

    f_vector_t timeArray;
    timeArray = Slice->bin_centers - Slice->bin_centers[0];

    auto TWC = TravelingWaveCavity(RS, frequency, a_factor);
    TWC.wake_calc(timeArray);

    util::read_vector_from_file(v, params + "wake.txt");
    ASSERT_NEAR_LOOP(v, TWC.fWake, "wake", epsilon);
}

TEST_F(testTWC, wake_calc2)
{
    auto Slice = Context::Slice;
    std::string params = TEST_FILES "/Impedances/Intensity/TWC/wake_calc2/";
    f_vector_t v;
    double epsilon = 1e-8;

    f_vector_t RS, frequency, a_factor;
    f_vector_t random;
    util::read_vector_from_file(random, TEST_FILES "/normal_distribution.dat");

    RS = slice(random, 0, 100) * 1e6;
    apply_f_in_place(RS, std::abs<double>);
    frequency = slice(random, 100, 200) * 1e9;
    apply_f_in_place(frequency, std::abs<double>);
    a_factor = slice(random, 200, 300);
    apply_f_in_place(a_factor, std::abs<double>);

    f_vector_t timeArray;
    timeArray = Slice->bin_centers - Slice->bin_centers[0];

    auto TWC = TravelingWaveCavity(RS, frequency, a_factor);
    TWC.wake_calc(timeArray);

    util::read_vector_from_file(v, params + "wake.txt");
    ASSERT_NEAR_LOOP(v, TWC.fWake, "wake", epsilon);
}

TEST_F(testTWC, imped_calc1)
{
    auto Slice = Context::Slice;
    std::string params = TEST_FILES "/Impedances/Intensity/TWC/imped_calc1/";
    f_vector_t v;
    double epsilon = 1e-8;

    f_vector_t RS, frequency, a_factor;
    f_vector_t random;
    util::read_vector_from_file(random, TEST_FILES "/normal_distribution.dat");

    RS = slice(random, 0, 10) * 1e6;
    apply_f_in_place(RS, std::abs<double>);
    frequency = slice(random, 10, 20) * 1e9;
    apply_f_in_place(frequency, std::abs<double>);
    a_factor = slice(random, 20, 30);
    apply_f_in_place(a_factor, std::abs<double>);

    f_vector_t timeArray;
    timeArray = Slice->bin_centers - Slice->bin_centers[0];

    auto TWC = TravelingWaveCavity(RS, frequency, a_factor);
    TWC.imped_calc(timeArray);
    
    util::read_vector_from_file(v, params + "impedance_norm.txt");
    f_vector_t norm;
    for(auto &z: TWC.fImpedance) norm.push_back(abs(z));
    ASSERT_NEAR_LOOP(v, norm, "impedance_norm", epsilon);

}

TEST_F(testTWC, imped_calc2)
{
    auto Slice = Context::Slice;
    std::string params = TEST_FILES "/Impedances/Intensity/TWC/imped_calc2/";
    f_vector_t v;
    double epsilon = 1e-8;

    f_vector_t RS, frequency, a_factor;
    f_vector_t random;
    util::read_vector_from_file(random, TEST_FILES "/normal_distribution.dat");

    RS = slice(random, 0, 100) * 1e6;
    apply_f_in_place(RS, std::abs<double>);
    frequency = slice(random, 100, 200) * 1e9;
    apply_f_in_place(frequency, std::abs<double>);
    a_factor = slice(random, 200, 300);
    apply_f_in_place(a_factor, std::abs<double>);

    f_vector_t timeArray;
    timeArray = Slice->bin_centers - Slice->bin_centers[0];

    auto TWC = TravelingWaveCavity(RS, frequency, a_factor);
    TWC.imped_calc(timeArray);

    util::read_vector_from_file(v, params + "impedance_norm.txt");
    f_vector_t norm;
    for(auto &z: TWC.fImpedance) norm.push_back(abs(z));
    ASSERT_NEAR_LOOP(v, norm, "impedance_norm", epsilon);
}


int main(int ac, char *av[])
{
    ::testing::InitGoogleTest(&ac, av);
    return RUN_ALL_TESTS();
}

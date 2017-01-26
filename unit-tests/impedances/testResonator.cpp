#include <blond/blond.h>

#include <gtest/gtest.h>

const std::string datafiles = DEMO_FILES "/TC5_Wake_impedance/";

// Simulation parameters
// --------------------------------------------------------
// Bunch parameters
const long long N_b = 1e10; // Intensity
const double tau_0 = 2e-9;  // Initial bunch length, 4 sigma [s]
// const particle_type particle = proton;
// Machine and RF parameters
const double C = 6911.56;   // Machine circumference [m]
const double p_i = 25.92e9; // Synchronous momentum [eV/c]
// const double p_f = 460.005e9;                  // Synchronous momentum, final
const long long h = 4620;                     // Harmonic number
const double V = 0.9e6;                        // RF voltage [V]
const double dphi = 0;                         // Phase modulation/offset
const double gamma_t = 1 / std::sqrt(0.00192); // Transition gamma
const double alpha = 1.0 / gamma_t / gamma_t;  // First order mom. comp. factor
const int alpha_order = 1;
const int n_sections = 1;
// Tracking details

int N_t = 2;       // Number of turns to track
int N_p = 5000000; // Macro-particles

int N_slices = 1 << 8; // = (2^8)

// RingAndRfSection *long_tracker;
Resonators *resonator;

class testResonator : public ::testing::Test {

protected:
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

        auto RfP = Context::RfP = new RfParameters(GP, n_sections, hVec,
                voltageVec, dphiVec);


        // RingAndRfSection *long_tracker = new RingAndRfSection();

        longitudinal_bigaussian(GP, RfP, Beam, tau_0 / 4, 0, -1, false);

        Context::Slice = new Slices(RfP, Beam, N_slices, 0, 0,
                                    2 * constant::pi, Slices::cuts_unit_t::rad);
        // util::dump(Slice->bin_centers, 10, "bin_centers\n");

        std::vector<double> v;
        util::read_vector_from_file(v, datafiles + "TC5_new_HQ_table.dat");
        assert(v.size() % 3 == 0);

        std::vector<double> R_shunt, f_res, Q_factor;

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
        // delete long_tracker;
    }
};

TEST_F(testResonator, initializations)
{
    std::string params = TEST_FILES "/Impedances/Intensity/";

    std::vector<double> v;
    util::read_vector_from_file(v, params + "R_shunt.txt");

    ASSERT_EQ(v.size(), resonator->fRS.size());

    double epsilon = 1e-8;

    for (unsigned int i = 0; i < v.size(); ++i) {
        double ref = v[i];
        double real = resonator->fRS[i];
        ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
                << "Testing of fRS failed on i " << i << std::endl;
    }
    v.clear();

    util::read_vector_from_file(v, params + "Q_factor.txt");
    ASSERT_EQ(v.size(), resonator->fQ.size());

    epsilon = 1e-8;
    for (unsigned int i = 0; i < v.size(); ++i) {
        double ref = v[i];
        double real = resonator->fQ[i];
        ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
                << "Testing of fQ failed on i " << i << std::endl;
    }
    v.clear();

    util::read_vector_from_file(v, params + "f_res.txt");
    ASSERT_EQ(v.size(), resonator->fFrequencyR.size());

    epsilon = 1e-8;
    for (unsigned int i = 0; i < v.size(); ++i) {
        double ref = v[i];
        double real = resonator->fFrequencyR[i];
        ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
                << "Testing of fFrequencyR failed on i " << i << std::endl;
    }
    v.clear();

    util::read_vector_from_file(v, params + "omega_R.txt");
    ASSERT_EQ(v.size(), resonator->fOmegaR.size());

    epsilon = 1e-8;
    for (unsigned int i = 0; i < v.size(); ++i) {
        double ref = v[i];
        double real = resonator->fOmegaR[i];
        ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
                << "Testing of fOmegaR failed on i " << i << std::endl;
    }
    v.clear();
}

TEST_F(testResonator, wake_calc)
{
    auto Slice = Context::Slice;
    std::string params = TEST_FILES "/Impedances/Intensity/";

    std::vector<double> v;
    util::read_vector_from_file(v, params + "Wake.txt");

    std::vector<double> timeArray;
    timeArray.reserve(N_slices);
    for (int i = 0; i < N_slices; ++i) {
        timeArray.push_back(Slice->bin_centers[i] - Slice->bin_centers[0]);
    }
    resonator->wake_calc(timeArray);
    // util::dump(&resonator->fWake[0], 10, "Wake\n");

    ASSERT_EQ(v.size(), resonator->fWake.size());

    double epsilon = 1e-8;

    for (unsigned int i = 0; i < v.size(); ++i) {
        double ref = v[i];
        double real = resonator->fWake[i];
        ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
                << "Testing of fWake failed on i " << i << std::endl;
    }
    v.clear();
}

TEST_F(testResonator, imped_calc)
{
    auto Slice = Context::Slice;

    std::string params = TEST_FILES "/Impedances/Intensity/";

    std::vector<double> v;
    util::read_vector_from_file(v, params + "Impedance.txt");

    std::vector<double> timeArray;
    timeArray.reserve(N_slices);
    for (int i = 0; i < N_slices; ++i) {
        timeArray.push_back((Slice->bin_centers[i] - Slice->bin_centers[0]) *
                            1e10);
    }

    // util::dump(&timeArray[0], 10, "timeArray\n");

    resonator->imped_calc(timeArray);
    // util::dump(&resonator->fImpedance[0], 10, "Impedance\n");

    ASSERT_EQ(v.size(), resonator->fImpedance.size());

    double epsilon = 1e-6;

    for (unsigned int i = 0; i < v.size(); ++i) {
        double ref = v[i];
        double real = std::abs((complex_t)resonator->fImpedance[i]);
        ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
                << "Testing of fImpedance failed on i " << i << std::endl;
    }
    v.clear();
}

int main(int ac, char *av[])
{
    ::testing::InitGoogleTest(&ac, av);
    return RUN_ALL_TESTS();
}

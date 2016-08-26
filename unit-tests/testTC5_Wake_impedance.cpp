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
#include <complex>
#include <gtest/gtest.h>
#include <omp.h>
#include <stdio.h>

const std::string datafiles = DEMO_FILES "/TC5_Wake_impedance/";

// Simulation parameters
// --------------------------------------------------------
// Bunch parameters
const long long int N_b = (long int)1e10; // Intensity
const ftype tau_0 = 2e-9;                 // Initial bunch length, 4 sigma [s]
// const particle_type particle = proton;
// Machine and RF parameters
const ftype C = 6911.56;   // Machine circumference [m]
const ftype p_i = 25.92e9; // Synchronous momentum [eV/c]
// const ftype p_f = 460.005e9;                  // Synchronous momentum, final
const long long h = 4620;                     // Harmonic number
const ftype V = 0.9e6;                        // RF voltage [V]
const ftype dphi = 0;                         // Phase modulation/offset
const ftype gamma_t = 1 / std::sqrt(0.00192); // Transition gamma
const ftype alpha = 1.0 / gamma_t / gamma_t;  // First order mom. comp. factor
const int alpha_order = 1;
const int n_sections = 1;
// Tracking details

unsigned N_t = 100; // Number of turns to track
unsigned N_p = 500; // Macro-particles

unsigned N_slices = 1 << 8; // = (2^8)

RingAndRfSection* long_tracker;
Resonators* resonator;

class testTC5 : public ::testing::Test {

  protected:
    virtual void SetUp() {

        omp_set_num_threads(Context::n_threads);

        f_vector_2d_t momentumVec(n_sections, f_vector_t(N_t + 1, p_i));

        // f_vector_2d_t alphaVec(alpha_order + 1, f_vector_t(n_sections,
        // alpha));
        f_vector_2d_t alphaVec(n_sections, f_vector_t(alpha_order + 1, alpha));

        f_vector_t CVec(n_sections, C);

        f_vector_2d_t hVec(n_sections, f_vector_t(N_t + 1, h));

        f_vector_2d_t voltageVec(n_sections, f_vector_t(N_t + 1, V));

        f_vector_2d_t dphiVec(n_sections, f_vector_t(N_t + 1, dphi));

        Context::GP = new GeneralParameters(N_t, CVec, alphaVec, alpha_order,
                                            momentumVec, proton);

        Context::Beam = new Beams(N_p, N_b);

        Context::RfP = new RfParameters(n_sections, hVec, voltageVec, dphiVec);

        long_tracker = new RingAndRfSection();

        longitudinal_bigaussian(tau_0 / 4, 0, -1, false);

        Context::Slice = new Slices(N_slices, 0, 0, 2 * constant::pi, rad);
        // util::dump(Slice->bin_centers, 10, "bin_centers\n");

        std::vector<ftype> v;
        util::read_vector_from_file(v, datafiles + "TC5_new_HQ_table.dat");
        assert(v.size() % 3 == 0);

        std::vector<ftype> R_shunt, f_res, Q_factor;

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

    virtual void TearDown() {
        // Code here will be called immediately after each test
        // (right before the destructor).
        delete Context::GP;
        delete Context::Beam;
        delete Context::RfP;
        delete Context::Slice;
        delete long_tracker;
        delete resonator;
    }
};

TEST_F(testTC5, timeTrack) {
    auto Beam = Context::Beam;
    auto Slice = Context::Slice;

    std::vector<Intensity*> wakeSourceList({resonator});
    InducedVoltageTime* indVoltTime = new InducedVoltageTime(wakeSourceList);
    std::vector<InducedVoltage*> indVoltList({indVoltTime});

    TotalInducedVoltage* totVol = new TotalInducedVoltage(indVoltList);

    for (unsigned i = 0; i < N_t; ++i) {
        totVol->track();
        long_tracker->track();
        Slice->track();
    }

    auto params = std::string(TEST_FILES "/") + "TC5_final/time/";

    std::vector<ftype> v;
    util::read_vector_from_file(v, params + "dE.txt");

    // WARNING checking only the fist 500 elems
    std::vector<ftype> res = Beam->dE;
    res.resize(500);
    ASSERT_EQ(v.size(), res.size());

    ftype epsilon = 1e-8;
    for (unsigned int i = 0; i < v.size(); ++i) {
        ftype ref = v[i];
        ftype real = res[i];

        ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
            << "Testing of Beam->dE failed on i " << i << std::endl;
    }

    v.clear();
    res.clear();
    util::read_vector_from_file(v, params + "dt.txt");

    // WARNING checking only the fist 500 elems
    res = std::vector<ftype>(Beam->dt.data(), Beam->dt.data() + 500);
    ASSERT_EQ(v.size(), res.size());

    epsilon = 1e-8;
    for (unsigned int i = 0; i < v.size(); ++i) {
        ftype ref = v[i];
        ftype real = res[i];

        ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
            << "Testing of Beam->dt failed on i " << i << std::endl;
    }

    v.clear();
    res.clear();
    util::read_vector_from_file(v, params + "n_macroparticles.txt");

    res = f_vector_t(Slice->n_macroparticles.begin(),
                     Slice->n_macroparticles.end());
    ASSERT_EQ(v.size(), res.size());

    epsilon = 1e-8;
    // warning checking only the first 100 elems
    for (unsigned int i = 0; i < v.size(); ++i) {
        ftype ref = v[i];
        ftype real = Slice->n_macroparticles[i];

        ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
            << "Testing of Slice->n_macroparticles failed on i " << i
            << std::endl;
    }

    delete indVoltTime;
    delete totVol;
}

TEST_F(testTC5, freqTrack) {
    auto Beam = Context::Beam;
    auto Slice = Context::Slice;

    std::vector<Intensity*> ImpSourceList({resonator});
    auto indVoltFreq = new InducedVoltageFreq(ImpSourceList, 1e5);
    std::vector<InducedVoltage*> indVoltList({indVoltFreq});

    TotalInducedVoltage* totVol = new TotalInducedVoltage(indVoltList);

    for (unsigned i = 0; i < N_t; ++i) {
        totVol->track();
        long_tracker->track();
        Slice->track();
    }

    auto params = std::string(TEST_FILES "/") + "TC5_final/freq/";

    std::vector<ftype> v;
    util::read_vector_from_file(v, params + "dE.txt");

    // WARNING checking only the fist 500 elems
    std::vector<ftype> res(Beam->dE.data(), Beam->dE.data() + 500);
    ASSERT_EQ(v.size(), res.size());

    ftype epsilon = 1e-8;
    // warning checking only the first 100 elems
    for (unsigned int i = 0; i < v.size(); ++i) {
        ftype ref = v[i];
        ftype real = res[i];

        ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
            << "Testing of Beam->dE failed on i " << i << std::endl;
    }

    v.clear();
    res.clear();
    util::read_vector_from_file(v, params + "dt.txt");

    // WARNING checking only the fist 500 elems
    res = std::vector<ftype>(Beam->dt.data(), Beam->dt.data() + 500);
    ASSERT_EQ(v.size(), res.size());

    epsilon = 1e-8;
    for (unsigned int i = 0; i < v.size(); ++i) {
        ftype ref = v[i];
        ftype real = res[i];

        ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
            << "Testing of Beam->dt failed on i " << i << std::endl;
    }

    v.clear();
    res.clear();
    util::read_vector_from_file(v, params + "n_macroparticles.txt");

    res = f_vector_t(Slice->n_macroparticles.begin(),
                     Slice->n_macroparticles.end());
    ASSERT_EQ(v.size(), res.size());

    epsilon = 1e-8;
    for (unsigned int i = 0; i < v.size(); ++i) {
        ftype ref = v[i];
        ftype real = Slice->n_macroparticles[i];

        ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
            << "Testing of Slice->n_macroparticles failed on i " << i
            << std::endl;
    }

    delete indVoltFreq;
    delete totVol;
}

int main(int ac, char* av[]) {
    ::testing::InitGoogleTest(&ac, av);
    return RUN_ALL_TESTS();
}

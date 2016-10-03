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
#include <blond/fft.h>
#include <gtest/gtest.h>
#include <stdio.h>


Resonators *resonator;

class testInducedVoltageFreq : public ::testing::Test {

protected:
    const long long int N_b = 1e10; // Intensity
    const ftype tau_0 = 2e-9;  // Initial bunch length, 4 sigma [s]
    const ftype C = 6911.56;   // Machine circumference [m]
    const ftype p_i = 25.92e9; // Synchronous momentum [eV/c]
    const long long h = 4620;  // Harmonic number
    const ftype V = 0.9e6;     // RF voltage [V]
    const ftype dphi = 0;      // Phase modulation/offset
    const ftype gamma_t = 1 / std::sqrt(0.00192); // Transition gamma
    const ftype alpha =
        1.0 / gamma_t / gamma_t; // First order mom. comp. factor
    const int alpha_order = 1;
    const int n_sections = 1;

    int N_t = 2;       // Number of turns to track
    int N_p = 5000000; // Macro-particles

    int N_slices = 1 << 8; // = (2^8)
    const std::string datafiles = DEMO_FILES "/TC5_Wake_impedance/";

    virtual void SetUp()
    {
        omp_set_num_threads(2);

        f_vector_2d_t momentumVec(n_sections, f_vector_t(N_t + 1, p_i));

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

        Context::Slice = new Slices(N_slices, 0, 0, 2 * constant::pi, rad);

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

TEST_F(testInducedVoltageFreq, constructor1)
{
    auto epsilon = 1e-8;

    std::vector<Intensity *> ImpSourceList({resonator});

    auto indVoltFreq = new InducedVoltageFreq(ImpSourceList, 1e5);

    auto params = std::string(TEST_FILES "/Impedances/") +
                  "InducedVoltage/InducedVoltageFreq/constructor1/";

    f_vector_t v;

    util::read_vector_from_file(v, params + "n_fft_sampling.txt");
    // ASSERT_EQ(v.size(), indVoltTime->fTotalWake.size());

    for (uint i = 0; i < v.size(); ++i) {
        auto ref = v[i];
        auto real = indVoltFreq->fNFFTSampling;
        ASSERT_NEAR(ref, real, epsilon * std::max(std::abs(ref), std::abs(real)))
                << "Testing of indVoltFreq->fNFFTSampling failed on i " << i
                << std::endl;
    }

    util::read_vector_from_file(v, params + "frequency_resolution.txt");
    for (uint i = 0; i < v.size(); ++i) {
        auto ref = v[i];
        auto real = indVoltFreq->fFreqResolution;
        ASSERT_NEAR(ref, real, epsilon * std::max(std::abs(ref), std::abs(real)))
                << "Testing of indVoltFreq->fFreqResolution failed on i " << i
                << std::endl;
    }

    util::read_vector_from_file(v, params + "frequency_array.txt");
    for (uint i = 0; i < v.size(); ++i) {
        auto ref = v[i];
        auto real = indVoltFreq->fFreqArray[i];
        ASSERT_NEAR(ref, real, epsilon * std::max(std::abs(ref), std::abs(real)))
                << "Testing of indVoltFreq->fFreqArray failed on i " << i
                << std::endl;
    }

    delete indVoltFreq;
}

TEST_F(testInducedVoltageFreq, constructor2)
{
    auto epsilon = 1e-8;

    std::vector<Intensity *> ImpSourceList({resonator});

    auto indVoltFreq =
        new InducedVoltageFreq(ImpSourceList, 1e5, round_option, 100);

    auto params = std::string(TEST_FILES "/Impedances/") +
                  "InducedVoltage/InducedVoltageFreq/constructor2/";

    f_vector_t v;

    util::read_vector_from_file(v, params + "n_turns_memory.txt");
    // ASSERT_EQ(v.size(), indVoltTime->fTotalWake.size());

    for (uint i = 0; i < v.size(); ++i) {
        auto ref = v[i];
        auto real = indVoltFreq->fNTurnsMem;
        ASSERT_NEAR(ref, real, epsilon * std::max(std::abs(ref), std::abs(real)))
                << "Testing of indVoltFreq->fNTurnsMem failed on i " << i
                << std::endl;
    }

    util::read_vector_from_file(v, params + "len_array_memory.txt");

    for (uint i = 0; i < v.size(); ++i) {
        auto ref = v[i];
        auto real = indVoltFreq->fLenArrayMem;
        ASSERT_NEAR(ref, real, epsilon * std::max(std::abs(ref), std::abs(real)))
                << "Testing of indVoltFreq->fLenArrayMem failed on i " << i
                << std::endl;
    }

    util::read_vector_from_file(v, params + "len_array_memory_extended.txt");

    for (uint i = 0; i < v.size(); ++i) {
        auto ref = v[i];
        auto real = indVoltFreq->fLenArrayMemExt;
        ASSERT_NEAR(ref, real, epsilon * std::max(std::abs(ref), std::abs(real)))
                << "Testing of indVoltFreq->fLenArrayMemExt failed on i " << i
                << std::endl;
    }

    util::read_vector_from_file(v, params + "n_points_fft.txt");

    for (uint i = 0; i < v.size(); ++i) {
        auto ref = v[i];
        auto real = indVoltFreq->fNPointsFFT;
        ASSERT_NEAR(ref, real, epsilon * std::max(std::abs(ref), std::abs(real)))
                << "Testing of indVoltFreq->fNPointsFFT failed on i " << i
                << std::endl;
    }

    util::read_vector_from_file(v, params + "frequency_array_memory.txt");
    // Only fist 1k elements of frequency_array_memory are tested
    for (uint i = 0; i < v.size(); ++i) {
        auto ref = v[i];
        auto real = indVoltFreq->fFreqArrayMem[i];
        ASSERT_NEAR(ref, real, epsilon * std::max(std::abs(ref), std::abs(real)))
                << "Testing of indVoltFreq->fFreqArrayMem failed on i " << i
                << std::endl;
    }

    util::read_vector_from_file(v, params + "time_array_memory.txt");
    // Only fist 100 elements of frequency_array_memory are tested
    for (uint i = 0; i < v.size(); ++i) {
        auto ref = v[i];
        auto real = indVoltFreq->fTimeArrayMem[i];
        ASSERT_NEAR(ref, real, epsilon * std::max(std::abs(ref), std::abs(real)))
                << "Testing of indVoltFreq->fTimeArrayMem failed on i " << i
                << std::endl;
    }

    util::read_vector_from_file(v, params + "total_impedance_memory.txt");
    // Only fist 1000 elements of total_impedance_memory are tested
    for (uint i = 0; i < v.size(); ++i) {
        auto ref = v[i];
        auto real = std::abs(indVoltFreq->fTotalImpedanceMem[i]);
        ASSERT_NEAR(ref, real, epsilon * std::max(std::abs(ref), std::abs(real)))
                << "Testing of indVoltFreq->fTotalImpedanceMem failed on i " << i
                << std::endl;
    }

    delete indVoltFreq;
}

TEST_F(testInducedVoltageFreq, sum_impedances1)
{
    auto epsilon = 1e-8;

    std::vector<Intensity *> ImpSourceList({resonator});

    auto indVoltFreq = new InducedVoltageFreq(ImpSourceList, 1e5);

    auto freq_array = fft::rfftfreq(Context::Slice->n_slices);

    indVoltFreq->sum_impedances(freq_array);

    auto params = std::string(TEST_FILES "/Impedances/") +
                  "InducedVoltage/InducedVoltageFreq/sum_impedances/";

    f_vector_t v;

    util::read_vector_from_file(v, params + "total_impedance.txt");

    ASSERT_EQ(v.size(), indVoltFreq->fTotalImpedance.size());

    for (uint i = 0; i < v.size(); ++i) {
        auto ref = v[i];
        auto real = std::abs(indVoltFreq->fTotalImpedance[i]);
        ASSERT_NEAR(ref, real, epsilon * std::max(std::abs(ref), std::abs(real)))
                << "Testing of indVoltFreq->fTotalImpedance failed on i " << i
                << std::endl;
    }

    delete indVoltFreq;
}

TEST_F(testInducedVoltageFreq, sum_impedances2)
{
    auto epsilon = 1e-8;
    auto Slice = Context::Slice;

    std::vector<Intensity *> ImpSourceList({resonator});

    auto indVoltFreq = new InducedVoltageFreq(ImpSourceList, 1e5);

    auto freq_array = fft::rfftfreq(Slice->n_slices, Slice->bin_centers[1] -
                                    Slice->bin_centers[0]);

    indVoltFreq->sum_impedances(freq_array);

    auto params = std::string(TEST_FILES "/Impedances/") +
                  "InducedVoltage/InducedVoltageFreq/sum_impedances2/";

    f_vector_t v;

    util::read_vector_from_file(v, params + "total_impedance.txt");
    ASSERT_EQ(v.size(), indVoltFreq->fTotalImpedance.size());
    for (uint i = 0; i < v.size(); ++i) {
        auto ref = v[i];
        auto real = std::abs(indVoltFreq->fTotalImpedance[i]);
        ASSERT_NEAR(ref, real, epsilon * std::max(std::abs(ref), std::abs(real)))
                << "Testing of indVoltFreq->fTotalImpedance failed on i " << i
                << std::endl;
    }

    delete indVoltFreq;
}

TEST_F(testInducedVoltageFreq, reprocess1)
{
    auto epsilon = 1e-8;
    auto Slice = Context::Slice;

    std::vector<Intensity *> ImpSourceList({resonator});

    auto indVoltFreq = new InducedVoltageFreq(ImpSourceList, 1e5);

    for (uint i = 0; i < Slice->n_slices; ++i)
        Slice->bin_centers[i] = 1.1 * Slice->bin_centers[i];

    indVoltFreq->reprocess();

    auto params = std::string(TEST_FILES "/Impedances/") +
                  "InducedVoltage/InducedVoltageFreq/reprocess1/";

    f_vector_t v;

    util::read_vector_from_file(v, params + "n_fft_sampling.txt");
    for (uint i = 0; i < v.size(); ++i) {
        auto ref = v[i];
        auto real = indVoltFreq->fNFFTSampling;
        ASSERT_NEAR(ref, real, epsilon * std::max(std::abs(ref), std::abs(real)))
                << "Testing of indVoltFreq->fNFFTSampling failed on i " << i
                << std::endl;
    }

    util::read_vector_from_file(v, params + "frequency_resolution.txt");
    for (uint i = 0; i < v.size(); ++i) {
        auto ref = v[i];
        auto real = indVoltFreq->fFreqResolution;
        ASSERT_NEAR(ref, real, epsilon * std::max(std::abs(ref), std::abs(real)))
                << "Testing of indVoltFreq->fFreqResolution failed on i " << i
                << std::endl;
    }

    util::read_vector_from_file(v, params + "frequency_array.txt");
    // only first 1k elements of frequency_array are tested
    for (uint i = 0; i < v.size(); ++i) {
        auto ref = v[i];
        auto real = indVoltFreq->fFreqArray[i];
        ASSERT_NEAR(ref, real, epsilon * std::max(std::abs(ref), std::abs(real)))
                << "Testing of indVoltFreq->fFreqArray failed on i " << i
                << std::endl;
    }

    delete indVoltFreq;
}


TEST_F(testInducedVoltageFreq, induced_voltage_generation1)
{
    auto epsilon = 1e-8;
    std::vector<Intensity *> ImpSourceList({resonator});

    auto indVoltFreq = new InducedVoltageFreq(ImpSourceList, 1e5);
    Context::Slice->track();

    indVoltFreq->induced_voltage_generation();
    auto params =
        std::string(TEST_FILES "/Impedances/") +
        "InducedVoltage/InducedVoltageFreq/induced_voltage_generation1/";

    f_vector_t v;

    util::read_vector_from_file(v, params + "induced_voltage.txt");

    ASSERT_EQ(v.size(), indVoltFreq->fInducedVoltage.size());

    for (uint i = 0; i < v.size(); ++i) {
        auto ref = v[i];
        ftype real = indVoltFreq->fInducedVoltage[i];
        ASSERT_NEAR(ref, real, epsilon * std::max(std::abs(ref), std::abs(real)))
                << "Testing of indVoltFreq->fInducedVoltage failed on i " << i
                << std::endl;
    }

    delete indVoltFreq;
}


TEST_F(testInducedVoltageFreq, induced_voltage_generation2)
{
    auto epsilon = 1e-8;
    std::vector<Intensity *> ImpSourceList({resonator});

    auto indVoltFreq = new InducedVoltageFreq(ImpSourceList, 1e4);
    Context::Slice->track();

    auto res = indVoltFreq->induced_voltage_generation(50);
    auto params =
        std::string(TEST_FILES "/Impedances/") +
        "InducedVoltage/InducedVoltageFreq/induced_voltage_generation2/";

    f_vector_t v;

    util::read_vector_from_file(v, params + "induced_voltage.txt");

    ASSERT_EQ(v.size(), res.size());

    for (uint i = 0; i < v.size(); ++i) {
        auto ref = v[i];
        ftype real = res[i];
        ASSERT_NEAR(ref, real, epsilon * std::max(std::abs(ref), std::abs(real)))
                << "Testing of indVoltFreq->fInducedVoltage failed on i " << i
                << std::endl;
    }

    delete indVoltFreq;
}


TEST_F(testInducedVoltageFreq, track1)
{
    auto epsilon = 1e-8;
    std::vector<Intensity *> ImpSourceList({resonator});

    auto indVoltFreq = new InducedVoltageFreq(ImpSourceList, 1e5);
    Context::Slice->track();

    indVoltFreq->track();
    auto params = std::string(TEST_FILES "/Impedances/") +
                  "InducedVoltage/InducedVoltageFreq/track1/";

    f_vector_t v;

    util::read_vector_from_file(v, params + "dE.txt");

    // ASSERT_EQ(v.size(), Beam->dE.size());
    // only testing 1k particles

    for (uint i = 0; i < v.size(); ++i) {
        auto ref = v[i];
        ftype real = Context::Beam->dE[i];
        ASSERT_NEAR(ref, real, epsilon * std::max(std::abs(ref), std::abs(real)))
                << "Testing of Beam->dE failed on i " << i << std::endl;
    }

    delete indVoltFreq;
}



TEST_F(testInducedVoltageFreq, track2)
{

    std::vector<Intensity *> ImpSourceList({resonator});

    auto indVoltFreq = new InducedVoltageFreq(ImpSourceList, 1e5);
    Context::Slice->track();

    for (auto i = 0; i < 10; i++)
        indVoltFreq->track();

    auto params = std::string(TEST_FILES "/Impedances/") +
                  "InducedVoltage/InducedVoltageFreq/track2/";

    f_vector_t v;

    util::read_vector_from_file(v, params + "dE.txt");

    // ASSERT_EQ(v.size(), Beam->dE.size());
    // only testing 1k particles

    auto epsilon = 1e-8;
    for (uint i = 0; i < v.size(); ++i) {
        auto ref = v[i];
        ftype real = Context::Beam->dE[i];
        ASSERT_NEAR(ref, real, epsilon * std::max(std::abs(ref), std::abs(real)))
                << "Testing of Beam->dE failed on i " << i << std::endl;
    }

    delete indVoltFreq;
}

int main(int ac, char *av[])
{
    ::testing::InitGoogleTest(&ac, av);
    return RUN_ALL_TESTS();
}

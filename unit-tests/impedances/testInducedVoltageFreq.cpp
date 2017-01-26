#include <blond/blond.h>
#include <gtest/gtest.h>

Resonators *resonator;

class testInducedVoltageFreq : public ::testing::Test {

protected:
    const long long int N_b = 1e10; // Intensity
    const double tau_0 = 2e-9;  // Initial bunch length, 4 sigma [s]
    const double C = 6911.56;   // Machine circumference [m]
    const double p_i = 25.92e9; // Synchronous momentum [eV/c]
    const long long h = 4620;  // Harmonic number
    const double V = 0.9e6;     // RF voltage [V]
    const double dphi = 0;      // Phase modulation/offset
    const double gamma_t = 1 / std::sqrt(0.00192); // Transition gamma
    const double alpha =
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

TEST_F(testInducedVoltageFreq, constructor1)
{
    auto epsilon = 1e-8;
    auto slices = Context::Slice;
    auto indVoltFreq = new InducedVoltageFreq(slices, {resonator}, 1e5);

    auto params = std::string(TEST_FILES "/Impedances/") +
                  "InducedVoltage/InducedVoltageFreq/constructor1/";

    f_vector_t v;

    util::read_vector_from_file(v, params + "n_fft_sampling.txt");
    // ASSERT_EQ(v.size(), indVoltTime->fTotalWake.size());

    for (uint i = 0; i < v.size(); ++i) {
        int ref = v[i];
        int real = indVoltFreq->fNFFTSampling;
        ASSERT_EQ(ref, real)
                << "Testing of indVoltFreq->fNFFTSampling failed\n";
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
    auto slices = Context::Slice;
    auto indVoltFreq = new InducedVoltageFreq(slices, {resonator},
            1e5, InducedVoltageFreq::round_option, 100);

    auto params = std::string(TEST_FILES "/Impedances/") +
                  "InducedVoltage/InducedVoltageFreq/constructor2/";

    f_vector_t v;

    util::read_vector_from_file(v, params + "n_turns_memory.txt");
    // ASSERT_EQ(v.size(), indVoltTime->fTotalWake.size());

    for (uint i = 0; i < v.size(); ++i) {
        int ref = v[i];
        int real = indVoltFreq->fNTurnsMem;
        ASSERT_EQ(ref, real)
                << "Testing of indVoltFreq->fNTurnsMem failed\n";
    }

    util::read_vector_from_file(v, params + "len_array_memory.txt");

    for (uint i = 0; i < v.size(); ++i) {
        int ref = v[i];
        int real = indVoltFreq->fLenArrayMem;
        ASSERT_EQ(ref, real)
                << "Testing of indVoltFreq->fLenArrayMem failed on i " << i
                << std::endl;
    }

    util::read_vector_from_file(v, params + "len_array_memory_extended.txt");

    for (uint i = 0; i < v.size(); ++i) {
        int ref = v[i];
        int real = indVoltFreq->fLenArrayMemExt;
        ASSERT_EQ(ref, real)
                << "Testing of indVoltFreq->fLenArrayMemExt failed\n";
    }

    util::read_vector_from_file(v, params + "n_points_fft.txt");

    for (uint i = 0; i < v.size(); ++i) {
        int ref = v[i];
        int real = indVoltFreq->fNPointsFFT;
        ASSERT_EQ(ref, real)
                << "Testing of indVoltFreq->fNPointsFFT failed\n";
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
    auto slices = Context::Slice;

    auto indVoltFreq = new InducedVoltageFreq(slices, {resonator}, 1e5);

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

    auto indVoltFreq = new InducedVoltageFreq(Slice, {resonator}, 1e5);

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

    auto indVoltFreq = new InducedVoltageFreq(Slice, {resonator}, 1e5);

    for (int i = 0; i < Slice->n_slices; ++i)
        Slice->bin_centers[i] = 1.1 * Slice->bin_centers[i];

    indVoltFreq->reprocess(Slice);

    auto params = std::string(TEST_FILES "/Impedances/") +
                  "InducedVoltage/InducedVoltageFreq/reprocess1/";

    f_vector_t v;

    util::read_vector_from_file(v, params + "n_fft_sampling.txt");
    for (uint i = 0; i < v.size(); ++i) {
        int ref = v[i];
        int real = indVoltFreq->fNFFTSampling;
        ASSERT_EQ(ref, real)
                << "Testing of indVoltFreq->fNFFTSampling failed\n";
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
    auto beam = Context::Beam;
    auto slices = Context::Slice;
    auto indVoltFreq = new InducedVoltageFreq(slices, {resonator}, 1e5);
    slices->track();

    indVoltFreq->induced_voltage_generation(beam);
    auto params =
        std::string(TEST_FILES "/Impedances/") +
        "InducedVoltage/InducedVoltageFreq/induced_voltage_generation1/";

    f_vector_t v;

    util::read_vector_from_file(v, params + "induced_voltage.txt");

    ASSERT_EQ(v.size(), indVoltFreq->fInducedVoltage.size());

    for (uint i = 0; i < v.size(); ++i) {
        auto ref = v[i];
        double real = indVoltFreq->fInducedVoltage[i];
        ASSERT_NEAR(ref, real, epsilon * std::max(std::abs(ref), std::abs(real)))
                << "Testing of indVoltFreq->fInducedVoltage failed on i " << i
                << std::endl;
    }

    delete indVoltFreq;
}


TEST_F(testInducedVoltageFreq, induced_voltage_generation2)
{
    auto epsilon = 1e-8;
    auto slices = Context::Slice;
    auto beam = Context::Beam;

    auto indVoltFreq = new InducedVoltageFreq(slices, {resonator}, 1e4);
    slices->track();

    auto res = indVoltFreq->induced_voltage_generation(beam, 50);
    auto params =
        std::string(TEST_FILES "/Impedances/") +
        "InducedVoltage/InducedVoltageFreq/induced_voltage_generation2/";

    f_vector_t v;

    util::read_vector_from_file(v, params + "induced_voltage.txt");

    ASSERT_EQ(v.size(), res.size());

    for (uint i = 0; i < v.size(); ++i) {
        auto ref = v[i];
        double real = res[i];
        ASSERT_NEAR(ref, real, epsilon * std::max(std::abs(ref), std::abs(real)))
                << "Testing of indVoltFreq->fInducedVoltage failed on i " << i
                << std::endl;
    }

    delete indVoltFreq;
}


TEST_F(testInducedVoltageFreq, track1)
{
    auto slices = Context::Slice;
    auto beam = Context::Beam;

    auto epsilon = 1e-8;

    auto indVoltFreq = new InducedVoltageFreq(slices, {resonator}, 1e5);
    slices->track();

    indVoltFreq->track(beam);
    auto params = std::string(TEST_FILES "/Impedances/") +
                  "InducedVoltage/InducedVoltageFreq/track1/";

    f_vector_t v;

    util::read_vector_from_file(v, params + "dE.txt");

    // ASSERT_EQ(v.size(), Beam->dE.size());
    // only testing 1k particles

    for (uint i = 0; i < v.size(); ++i) {
        auto ref = v[i];
        double real = Context::Beam->dE[i];
        ASSERT_NEAR(ref, real, epsilon * std::max(std::abs(ref), std::abs(real)))
                << "Testing of Beam->dE failed on i " << i << std::endl;
    }

    delete indVoltFreq;
}



TEST_F(testInducedVoltageFreq, track2)
{
    auto slices = Context::Slice;
    auto beam = Context::Beam;

    auto indVoltFreq = new InducedVoltageFreq(slices, {resonator}, 1e5);
    slices->track();

    for (auto i = 0; i < 10; i++)
        indVoltFreq->track(beam);

    auto params = std::string(TEST_FILES "/Impedances/") +
                  "InducedVoltage/InducedVoltageFreq/track2/";

    f_vector_t v;

    util::read_vector_from_file(v, params + "dE.txt");

    // ASSERT_EQ(v.size(), Beam->dE.size());
    // only testing 1k particles

    auto epsilon = 1e-8;
    for (uint i = 0; i < v.size(); ++i) {
        auto ref = v[i];
        double real = Context::Beam->dE[i];
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

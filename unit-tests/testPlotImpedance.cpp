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
#include <blond/plots/plot_impedance.h>
#include <gtest/gtest.h>
#include <stdio.h>

using namespace std;


class testPlotImpedance : public ::testing::Test {

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

        Context::GP = new GeneralParameters(N_t, CVec, alphaVec, alpha_order,
                                            momentumVec, proton);

        Context::Beam = new Beams(N_p, N_b);

        Context::RfP = new RfParameters(n_sections, hVec, voltageVec, dphiVec);

        longitudinal_bigaussian(tau_0 / 4, 0, -1, false);

        Context::Slice = new Slices(N_slices, 0, 0, 2 * constant::pi,
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

TEST_F(testPlotImpedance, plot_impedance_vs_frequency1)
{
    auto slice = Context::Slice;

    std::vector<Intensity *> ImpSourceList({resonator});
    auto indVoltFreq = new InducedVoltageFreq(ImpSourceList, 1e5);

    ASSERT_EQ(plot_impedance_vs_frequency(0, indVoltFreq, slice), 1);

    delete indVoltFreq;
}


TEST_F(testPlotImpedance, plot_impedance_vs_frequency2)
{
    auto slice = Context::Slice;

    f_vector_t timeArray;
    timeArray.reserve(N_slices);
    for (int i = 0; i < N_slices; ++i)
        timeArray.push_back(slice->bin_centers[i]
                            - slice->bin_centers[0]);

    resonator->wake_calc(timeArray);
    auto inputTable = new InputTable(timeArray, resonator->fWake,
                                     resonator->fWake);

    std::vector<Intensity *> ImpSourceList({resonator});
    auto indVoltFreq = new InducedVoltageFreq(ImpSourceList, 1e5);

    ASSERT_EQ(plot_impedance_vs_frequency(0, indVoltFreq, slice, "sum",
                                          "no_spectrum", "freq_table"), 1);


    delete indVoltFreq;
    delete inputTable;
}

TEST_F(testPlotImpedance, plot_induced_voltage_vs_bin_centers)
{
    auto slice = Context::Slice;

    std::vector<Intensity *> wakeSourceList({resonator});
    auto indVoltTime = new InducedVoltageTime(wakeSourceList);
    std::vector<InducedVoltage *> indVoltList({indVoltTime});
    auto totVol = new TotalInducedVoltage(indVoltList);
    totVol->track();

    ASSERT_EQ(plot_induced_voltage_vs_bin_centers(0, totVol, slice), 1);


    delete indVoltTime;
    delete totVol;

}



int main(int ac, char *av[])
{
    python::initialize();
    ::testing::InitGoogleTest(&ac, av);
    auto ret = RUN_ALL_TESTS();
    python::finalize();
    return ret;
}

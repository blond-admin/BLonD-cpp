#include <blond/beams/Distributions.h>
#include <blond/constants.h>
#include <blond/input_parameters/GeneralParameters.h>
#include <blond/math_functions.h>
#include <blond/trackers/Tracker.h>
#include <blond/utilities.h>
#include <gtest/gtest.h>
#include <testing_utilities.h>
#include <blond/plots/plot_slices.h>
#include <blond/vector_math.h>

using namespace std;
using namespace mymath;


class testDistributions : public ::testing::Test {

protected:
    // Bunch parameters
    const long long N_b = 0; // Intensity

    // Machine and RF parameters
    const double tau_0 = 0.4e-9;
    const double radius = 25;
    const double C = 2 * constant::pi * radius; // Machine circumference [m]
    const double p_i = 310891054.809;           // Synchronous momentum [eV/c]
    const double h = 1;                          // Harmonic number
    const double V = 8000;                      // RF voltage [V]
    const double dphi = -constant::pi;          // Phase modulation/offset
    const double gamma_t = 4.076750841;         // Transition gamma
    const double alpha =
        1.0 / gamma_t / gamma_t; // First order mom. comp. factor
    const int alpha_order = 1;
    const int n_sections = 1;
    // Tracking details

    int N_t = 1000;  // Number of turns to track
    int N_p = 1000; // Macro-particles

    int N_slices = 100; // = (2^8)

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


        auto GP = Context::GP;
        auto Beam = Context::Beam = new Beams(GP, N_p, N_b);

        auto RfP = Context::RfP = new RfParameters(GP, n_sections, hVec,
                voltageVec, dphiVec);

    }

    virtual void TearDown()
    {
        // Code here will be called immediately after each test
        // (right before the destructor).
        delete Context::GP;
        delete Context::Beam;
        delete Context::RfP;
    }
};

TEST(testHelpers, line_density_function1)
{
    auto epsilon = 1e-8;
    auto params = string(TEST_FILES "/Distributions/line_density_function1/");

    int size = 100;
    f_vector_t coord_array(size);
    for (int i = 0; i < size; i++) {
        coord_array[i] = cos(i);
    }

    string dist_type = "waterbag";
    double bunch_length = 1.5;
    auto ret = line_density_function(coord_array, dist_type, bunch_length);

    f_vector_t v;
    util::read_vector_from_file(v, params + "density.txt");

    ASSERT_EQ(v.size(), ret.size());
    for (uint i = 0; i < v.size(); ++i) {
        auto ref = v[i];
        auto real = ret[i];
        ASSERT_NEAR(ref, real, epsilon * max(abs(ref), abs(real)))
                << "Testing of density failed on i " << i << endl;
    }


}


TEST(testHelpers, line_density_function2)
{
    auto epsilon = 1e-8;
    auto params = string(TEST_FILES "/Distributions/line_density_function2/");
    int size = 100;

    f_vector_t coord_array(size);
    for (int i = 0; i < size; i++) {
        coord_array[i] = cos(i);
    }

    string dist_type = "binomial";
    double bunch_length = 1.;
    auto ret = line_density_function(coord_array, dist_type, bunch_length,
                                     0.5, 0.6);
    f_vector_t v;
    util::read_vector_from_file(v, params + "density.txt");

    ASSERT_EQ(v.size(), ret.size());
    for (uint i = 0; i < v.size(); ++i) {
        auto ref = v[i];
        auto real = ret[i];
        ASSERT_NEAR(ref, real, epsilon * max(abs(ref), abs(real)))
                << "Testing of density failed on i " << i << endl;
    }

}

TEST(testHelpers, line_density_function3)
{
    auto epsilon = 1e-8;
    auto params = string(TEST_FILES "/Distributions/line_density_function3/");
    int size = 100;

    f_vector_t coord_array(size);
    for (int i = 0; i < size; i++) {
        coord_array[i] = cos(i);
    }

    string dist_type = "parabolic_line";
    double bunch_length = 2.;
    auto ret = line_density_function(coord_array, dist_type, bunch_length);
    f_vector_t v;
    util::read_vector_from_file(v, params + "density.txt");

    ASSERT_EQ(v.size(), ret.size());
    for (uint i = 0; i < v.size(); ++i) {
        auto ref = v[i];
        auto real = ret[i];
        ASSERT_NEAR(ref, real, epsilon * max(abs(ref), abs(real)))
                << "Testing of density failed on i " << i << endl;
    }

}

TEST(testHelpers, line_density_function4)
{
    auto epsilon = 1e-8;
    auto params = string(TEST_FILES "/Distributions/line_density_function4/");
    int size = 100;

    f_vector_t coord_array(size);
    for (int i = 0; i < size; i++) {
        coord_array[i] = cos(i);
    }

    string dist_type = "gaussian";
    double bunch_length = 1.2;
    auto ret = line_density_function(coord_array, dist_type, bunch_length, 1.);
    f_vector_t v;
    util::read_vector_from_file(v, params + "density.txt");

    ASSERT_EQ(v.size(), ret.size());
    for (uint i = 0; i < v.size(); ++i) {
        auto ref = v[i];
        auto real = ret[i];
        ASSERT_NEAR(ref, real, epsilon * max(abs(ref), abs(real)))
                << "Testing of density failed on i " << i << endl;
    }

}


TEST(testHelpers, line_density_function5)
{
    auto epsilon = 1e-8;
    auto params = string(TEST_FILES "/Distributions/line_density_function5/");
    int size = 100;

    f_vector_t coord_array(size);
    for (int i = 0; i < size; i++) {
        coord_array[i] = cos(i);
    }

    string dist_type = "cosine_squared";
    double bunch_length = 1.2;
    auto ret = line_density_function(coord_array, dist_type, bunch_length, 1.);
    f_vector_t v;
    util::read_vector_from_file(v, params + "density.txt");

    ASSERT_EQ(v.size(), ret.size());
    for (uint i = 0; i < v.size(); ++i) {
        auto ref = v[i];
        auto real = ret[i];
        ASSERT_NEAR(ref, real, epsilon * max(abs(ref), abs(real)))
                << "Testing of density failed on i " << i << endl;
    }

}

TEST(testHelpers, line_density_function_deathtest1)
{
    f_vector_t x;
    ASSERT_DEATH(line_density_function(x, "blabla", 0.),
                 "[line_density_function]\\s*");
}


TEST(testHelpers, distribution_density_function_deathtest1)
{
    f_vector_t x;
    ASSERT_DEATH(distribution_density_function(x, "blabla", 0.),
                 "[distribution_density_function]\\s*");
}


TEST(testHelpers, distribution_density_function1)
{
    auto epsilon = 1e-8;
    auto params = string(TEST_FILES "/Distributions/distribution_density_function1/");
    int size = 100;

    f_vector_t coord_array(size);
    for (int i = 0; i < size; i++) {
        coord_array[i] = cos(i);
    }

    string dist_type = "binomial";
    double bunch_length = 2.;
    auto ret = distribution_density_function(coord_array, dist_type, bunch_length, 0.7);
    f_vector_t v;
    util::read_vector_from_file(v, params + "density.txt");

    ASSERT_EQ(v.size(), ret.size());
    for (uint i = 0; i < v.size(); ++i) {
        auto ref = v[i];
        auto real = ret[i];
        ASSERT_NEAR(ref, real, epsilon * max(abs(ref), abs(real)))
                << "Testing of density failed on i " << i << endl;
    }

}


TEST(testHelpers, distribution_density_function2)
{
    auto epsilon = 1e-8;
    auto params = string(TEST_FILES "/Distributions/distribution_density_function2/");
    int size = 100;

    f_vector_t coord_array(size);
    for (int i = 0; i < size; i++) {
        coord_array[i] = cos(i);
    }

    string dist_type = "waterbag";
    double bunch_length = 2.;
    auto ret = distribution_density_function(coord_array, dist_type, bunch_length);
    f_vector_t v;
    util::read_vector_from_file(v, params + "density.txt");

    ASSERT_EQ(v.size(), ret.size());
    for (uint i = 0; i < v.size(); ++i) {
        auto ref = v[i];
        auto real = ret[i];
        ASSERT_NEAR(ref, real, epsilon * max(abs(ref), abs(real)))
                << "Testing of density failed on i " << i << endl;
    }

}

TEST(testHelpers, distribution_density_function3)
{
    auto epsilon = 1e-8;
    auto params = string(TEST_FILES "/Distributions/distribution_density_function3/");
    int size = 100;

    f_vector_t coord_array(size);
    for (int i = 0; i < size; i++) {
        coord_array[i] = cos(i);
    }

    string dist_type = "parabolic_line";
    double bunch_length = 2.;
    auto ret = distribution_density_function(coord_array, dist_type, bunch_length, 1.1);
    f_vector_t v;
    util::read_vector_from_file(v, params + "density.txt");

    ASSERT_EQ(v.size(), ret.size());
    for (uint i = 0; i < v.size(); ++i) {
        auto ref = v[i];
        auto real = ret[i];
        ASSERT_NEAR(ref, real, epsilon * max(abs(ref), abs(real)))
                << "Testing of density failed on i " << i << endl;
    }

}


TEST(testHelpers, distribution_density_function4)
{
    auto epsilon = 1e-8;
    auto params = string(TEST_FILES "/Distributions/distribution_density_function4/");
    int size = 100;

    f_vector_t coord_array(size);
    for (int i = 0; i < size; i++) {
        coord_array[i] = cos(i);
    }

    string dist_type = "gaussian";
    double bunch_length = 3.;
    auto ret = distribution_density_function(coord_array, dist_type, bunch_length);
    f_vector_t v;
    util::read_vector_from_file(v, params + "density.txt");

    ASSERT_EQ(v.size(), ret.size());
    for (uint i = 0; i < v.size(); ++i) {
        auto ref = v[i];
        auto real = ret[i];
        ASSERT_NEAR(ref, real, epsilon * max(abs(ref), abs(real)))
                << "Testing of density failed on i " << i << endl;
    }

}


TEST_F(testDistributions, matched_from_line_density1)
{
    auto RfP = Context::RfP;
    auto Beam = Context::Beam;
    auto Slice = new Slices(RfP, Beam, N_slices);
    auto epsilon = 1e-8;
    auto params = string(TEST_FILES "/Distributions/line_density1/");
    string datafiles = DEMO_FILES "/TC5_Wake_impedance/";


    f_vector_t v;
    util::read_vector_from_file(v, datafiles + "TC5_new_HQ_table.dat");

    f_vector_t R_shunt, f_res, Q_factor;

    R_shunt.reserve(v.size() / 3);
    f_res.reserve(v.size() / 3);
    Q_factor.reserve(v.size() / 3);

    for (uint i = 0; i < v.size(); i += 3) {
        f_res.push_back(v[i] * 1e9);
        Q_factor.push_back(v[i + 1]);
        R_shunt.push_back(v[i + 2] * 1e6);
    }

    auto resonator = new Resonators(R_shunt, f_res, Q_factor);
    auto indVoltTime = new InducedVoltageTime(Slice, {resonator});
    auto totVolt = new TotalInducedVoltage(Beam, Slice, {indVoltTime});

    auto long_tracker = new RingAndRfSection(RfP);
    auto fullRing = new FullRingAndRf({long_tracker});

    map<string, string> line_density_opt;
    line_density_opt["type"] = "gaussian";
    line_density_opt["bunch_length"] = "100e-9";
    line_density_opt["density_variable"] = "density_from_J";
    line_density_opt["exponent"] = "1.5";

    auto ret = matched_from_line_density(Beam, fullRing, line_density_opt,
                                         FullRingAndRf::lowest_freq, totVolt);

    util::read_vector_from_file(v, params + "hamiltonian_coord.txt");
    ASSERT_NEAR_LOOP(v, ret.hamiltonian_coord, "hamiltonian_coord", epsilon);

    util::read_vector_from_file(v, params + "density_function.txt");
    ASSERT_NEAR_LOOP(v, ret.density_function, "density_function", epsilon);

    util::read_vector_from_file(v, params + "time_line_den.txt");
    ASSERT_NEAR_LOOP(v, ret.time_line_den, "time_line_den", epsilon);

    util::read_vector_from_file(v, params + "line_density.txt");
    ASSERT_NEAR_LOOP(v, ret.line_density, "line_density", epsilon);

    delete totVolt;
    delete indVoltTime;
    delete resonator;
    delete Slice;
    delete long_tracker;
    delete fullRing;
}



TEST_F(testDistributions, matched_from_line_density2)
{
    auto RfP = Context::RfP;
    auto Beam = Context::Beam;
    auto Slice = new Slices(RfP, Beam, N_slices);
    auto epsilon = 1e-8;
    auto params = string(TEST_FILES "/Distributions/line_density2/");
    string datafiles = DEMO_FILES "/TC5_Wake_impedance/";


    f_vector_t v;
    util::read_vector_from_file(v, datafiles + "TC5_new_HQ_table.dat");

    f_vector_t R_shunt, f_res, Q_factor;

    R_shunt.reserve(v.size() / 3);
    f_res.reserve(v.size() / 3);
    Q_factor.reserve(v.size() / 3);

    for (uint i = 0; i < v.size(); i += 3) {
        f_res.push_back(v[i] * 1e9);
        Q_factor.push_back(v[i + 1]);
        R_shunt.push_back(v[i + 2] * 1e6);
    }

    auto resonator = new Resonators(R_shunt, f_res, Q_factor);
    auto indVoltTime = new InducedVoltageTime(Slice, {resonator});
    auto totVolt = new TotalInducedVoltage(Beam, Slice, {indVoltTime});

    auto long_tracker = new RingAndRfSection(RfP);
    auto fullRing = new FullRingAndRf({long_tracker});

    map<string, string> line_density_opt;
    line_density_opt["type"] = "parabolic_amplitude";
    line_density_opt["bunch_length"] = "200e-9";
    line_density_opt["density_variable"] = "density_from_J";
    // line_density_opt["exponent"] = "2.5";

    auto ret = matched_from_line_density(Beam, fullRing, line_density_opt,
                                         FullRingAndRf::lowest_freq, totVolt);

    util::read_vector_from_file(v, params + "hamiltonian_coord.txt");
    ASSERT_NEAR_LOOP(v, ret.hamiltonian_coord, "hamiltonian_coord", epsilon);

    util::read_vector_from_file(v, params + "density_function.txt");
    ASSERT_NEAR_LOOP(v, ret.density_function, "density_function", epsilon);

    util::read_vector_from_file(v, params + "time_line_den.txt");
    ASSERT_NEAR_LOOP(v, ret.time_line_den, "time_line_den", epsilon);

    util::read_vector_from_file(v, params + "line_density.txt");
    ASSERT_NEAR_LOOP(v, ret.line_density, "line_density", epsilon);

    delete totVolt;
    delete indVoltTime;
    delete resonator;
    delete Slice;
    delete long_tracker;
    delete fullRing;
}

TEST_F(testDistributions, matched_from_line_density3)
{
    auto RfP = Context::RfP;
    auto Beam = Context::Beam;
    auto Slice = new Slices(RfP, Beam, N_slices);
    auto epsilon = 1e-8;
    auto params = string(TEST_FILES "/Distributions/line_density3/");
    string datafiles = DEMO_FILES "/TC5_Wake_impedance/";


    f_vector_t v;
    util::read_vector_from_file(v, datafiles + "TC5_new_HQ_table.dat");

    f_vector_t R_shunt, f_res, Q_factor;

    R_shunt.reserve(v.size() / 3);
    f_res.reserve(v.size() / 3);
    Q_factor.reserve(v.size() / 3);

    for (uint i = 0; i < v.size(); i += 3) {
        f_res.push_back(v[i] * 1e9);
        Q_factor.push_back(v[i + 1]);
        R_shunt.push_back(v[i + 2] * 1e6);
    }

    auto resonator = new Resonators(R_shunt, f_res, Q_factor);
    auto indVoltTime = new InducedVoltageTime(Slice, {resonator});
    auto totVolt = new TotalInducedVoltage(Beam, Slice, {indVoltTime});

    auto long_tracker = new RingAndRfSection(RfP);
    auto fullRing = new FullRingAndRf({long_tracker});

    map<string, string> line_density_opt;
    line_density_opt["type"] = "cosine_squared";
    line_density_opt["bunch_length"] = "150e-9";
    line_density_opt["density_variable"] = "density_from_J";
    // line_density_opt["exponent"] = "2.5";

    auto ret = matched_from_line_density(Beam, fullRing, line_density_opt,
                                         FullRingAndRf::lowest_freq, totVolt,
                                         "savefig", "fig", "both");

    util::read_vector_from_file(v, params + "hamiltonian_coord.txt");
    ASSERT_NEAR_LOOP(v, ret.hamiltonian_coord, "hamiltonian_coord", epsilon);

    util::read_vector_from_file(v, params + "density_function.txt");
    ASSERT_NEAR_LOOP(v, ret.density_function, "density_function", epsilon);

    util::read_vector_from_file(v, params + "time_line_den.txt");
    ASSERT_NEAR_LOOP(v, ret.time_line_den, "time_line_den", epsilon);

    util::read_vector_from_file(v, params + "line_density.txt");
    ASSERT_NEAR_LOOP(v, ret.line_density, "line_density", epsilon);

    delete totVolt;
    delete indVoltTime;
    delete resonator;
    delete Slice;
    delete long_tracker;
    delete fullRing;
}

TEST_F(testDistributions, matched_from_line_density4)
{
    auto RfP = Context::RfP;
    auto Beam = Context::Beam;
    auto Slice = new Slices(RfP, Beam, N_slices);
    auto epsilon = 1e-8;
    auto params = string(TEST_FILES "/Distributions/line_density4/");
    string datafiles = DEMO_FILES "/TC5_Wake_impedance/";


    f_vector_t v;
    util::read_vector_from_file(v, datafiles + "TC5_new_HQ_table.dat");

    f_vector_t R_shunt, f_res, Q_factor;

    R_shunt.reserve(v.size() / 3);
    f_res.reserve(v.size() / 3);
    Q_factor.reserve(v.size() / 3);

    for (uint i = 0; i < v.size(); i += 3) {
        f_res.push_back(v[i] * 1e9);
        Q_factor.push_back(v[i + 1]);
        R_shunt.push_back(v[i + 2] * 1e6);
    }

    auto resonator = new Resonators(R_shunt, f_res, Q_factor);
    auto indVoltTime = new InducedVoltageTime(Slice, {resonator});
    auto totVolt = new TotalInducedVoltage(Beam, Slice, {indVoltTime});

    auto long_tracker = new RingAndRfSection(RfP);
    auto fullRing = new FullRingAndRf({long_tracker});

    map<string, string> line_density_opt;
    line_density_opt["type"] = "waterbag";
    line_density_opt["bunch_length"] = "150e-9";
    line_density_opt["density_variable"] = "density_from_J";
    // line_density_opt["exponent"] = "2.5";

    auto ret = matched_from_line_density(Beam, fullRing, line_density_opt,
                                         FullRingAndRf::lowest_freq, totVolt,
                                         "savefig", "fig", "second");

    util::read_vector_from_file(v, params + "hamiltonian_coord.txt");
    ASSERT_NEAR_LOOP(v, ret.hamiltonian_coord, "hamiltonian_coord", epsilon);

    util::read_vector_from_file(v, params + "density_function.txt");
    ASSERT_NEAR_LOOP(v, ret.density_function, "density_function", epsilon);

    util::read_vector_from_file(v, params + "time_line_den.txt");
    ASSERT_NEAR_LOOP(v, ret.time_line_den, "time_line_den", epsilon);

    util::read_vector_from_file(v, params + "line_density.txt");
    ASSERT_NEAR_LOOP(v, ret.line_density, "line_density", epsilon);

    delete totVolt;
    delete indVoltTime;
    delete resonator;
    delete Slice;
    delete long_tracker;
    delete fullRing;
}


TEST_F(testDistributions, matched_from_line_density5)
{
    auto RfP = Context::RfP;
    auto Beam = Context::Beam;
    auto Slice = new Slices(RfP, Beam, N_slices);
    auto epsilon = 1e-8;
    auto params = string(TEST_FILES "/Distributions/line_density5/");

    f_vector_t v;

    auto long_tracker = new RingAndRfSection(RfP);
    auto fullRing = new FullRingAndRf({long_tracker});

    map<string, string> line_density_opt;
    line_density_opt["type"] = "binomial";
    line_density_opt["bunch_length"] = "150e-9";
    line_density_opt["density_variable"] = "density_from_H";
    line_density_opt["exponent"] = "0.7";

    map<string, f_vector_t> extra_volt_dict;
    auto time_array = linspace(0.0, 1e-3, 10000);
    auto voltage_array = apply_f(time_array, [this](double x) {return sin(x) * V / 10.;});

    extra_volt_dict["time_array"] = time_array;
    extra_volt_dict["voltage_array"] = voltage_array;

    auto ret = matched_from_line_density(Beam, fullRing, line_density_opt,
                                         FullRingAndRf::lowest_freq, nullptr,
                                         "", "", "first", extra_volt_dict);

    util::read_vector_from_file(v, params + "hamiltonian_coord.txt");
    ASSERT_NEAR_LOOP(v, ret.hamiltonian_coord, "hamiltonian_coord", epsilon);

    util::read_vector_from_file(v, params + "density_function.txt");
    ASSERT_NEAR_LOOP(v, ret.density_function, "density_function", epsilon);

    util::read_vector_from_file(v, params + "time_line_den.txt");
    ASSERT_NEAR_LOOP(v, ret.time_line_den, "time_line_den", epsilon);

    util::read_vector_from_file(v, params + "line_density.txt");
    ASSERT_NEAR_LOOP(v, ret.line_density, "line_density", epsilon);

    delete Slice;
    delete long_tracker;
    delete fullRing;
}

TEST_F(testDistributions, matched_from_line_density_deathtest1)
{
    auto RfP = Context::RfP;
    auto Beam = Context::Beam;
    auto Slice = new Slices(RfP, Beam, N_slices);
    auto epsilon = 1e-8;

    auto long_tracker = new RingAndRfSection(RfP);
    auto fullRing = new FullRingAndRf({long_tracker});

    map<string, string> line_density_opt;


    ASSERT_DEATH(matched_from_line_density(Beam, fullRing, line_density_opt),
                 "\\w* was not recognized\n");


    delete Slice;
    delete long_tracker;
    delete fullRing;
}


TEST_F(testDistributions, matched_from_distribution_density1)
{
    auto RfP = Context::RfP;
    auto Beam = Context::Beam;
    auto Slice = new Slices(RfP, Beam, N_slices);
    auto epsilon = 1e-8;
    auto params = string(TEST_FILES "/Distributions/distribution_density1/");
    string datafiles = DEMO_FILES "/TC5_Wake_impedance/";


    f_vector_t v;
    util::read_vector_from_file(v, datafiles + "TC5_new_HQ_table.dat");

    f_vector_t R_shunt, f_res, Q_factor;

    R_shunt.reserve(v.size() / 3);
    f_res.reserve(v.size() / 3);
    Q_factor.reserve(v.size() / 3);

    for (uint i = 0; i < v.size(); i += 3) {
        f_res.push_back(v[i] * 1e9);
        Q_factor.push_back(v[i + 1]);
        R_shunt.push_back(v[i + 2] * 1e6);
    }

    auto resonator = new Resonators(R_shunt, f_res, Q_factor);
    auto indVoltTime = new InducedVoltageTime(Slice, {resonator});
    auto totVolt = new TotalInducedVoltage(Beam, Slice, {indVoltTime});

    auto long_tracker = new RingAndRfSection(RfP);
    auto fullRing = new FullRingAndRf({long_tracker});

    map<string, multi_t> distribution_opt;
    distribution_opt["type"] = {"gaussian"};
    distribution_opt["bunch_length"] = {250e-9};
    distribution_opt["density_variable"] = {"density_from_J"};
    distribution_opt["exponent"] = {2.5};
    // distribution_opt["bunch_length_fit"] = {"gauss"};


    // map<string, f_vector_t> extra_volt_dict;
    // extra_volt_dict["time_array"] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    // extra_volt_dict["voltage_array"] = {112, 2, 34, 41, 5, 6, 71, 8, 9, 10};

    auto ret = matched_from_distribution_density(
                   Beam, fullRing, distribution_opt,
                   FullRingAndRf::lowest_freq, totVolt);

    util::read_vector_from_file(v, params + "time_coord_low_res.txt");
    ASSERT_NEAR_LOOP(v, ret.time_coord_low_res, "time_coord_low_res", epsilon);

    util::read_vector_from_file(v, params + "line_density.txt");
    ASSERT_NEAR_LOOP(v, ret.line_density, "line_density", epsilon);

    delete totVolt;
    delete indVoltTime;
    delete resonator;
    delete Slice;
    delete long_tracker;
    delete fullRing;
}


TEST_F(testDistributions, matched_from_distribution_density2)
{
    auto RfP = Context::RfP;
    auto Beam = Context::Beam;
    auto Slice = new Slices(RfP, Beam, N_slices);
    auto epsilon = 1e-8;
    auto params = string(TEST_FILES "/Distributions/distribution_density2/");
    string datafiles = DEMO_FILES "/TC5_Wake_impedance/";


    f_vector_t v;
    util::read_vector_from_file(v, datafiles + "TC5_new_HQ_table.dat");

    f_vector_t R_shunt, f_res, Q_factor;

    R_shunt.reserve(v.size() / 3);
    f_res.reserve(v.size() / 3);
    Q_factor.reserve(v.size() / 3);

    for (uint i = 0; i < v.size(); i += 3) {
        f_res.push_back(v[i] * 1e9);
        Q_factor.push_back(v[i + 1]);
        R_shunt.push_back(v[i + 2] * 1e6);
    }

    auto resonator = new Resonators(R_shunt, f_res, Q_factor);
    auto indVoltTime = new InducedVoltageTime(Slice, {resonator});
    auto totVolt = new TotalInducedVoltage(Beam, Slice, {indVoltTime});

    auto long_tracker = new RingAndRfSection(RfP);
    auto fullRing = new FullRingAndRf({long_tracker});

    map<string, multi_t> distribution_opt;
    distribution_opt["type"] = {"parabolic_amplitude"};
    distribution_opt["bunch_length"] = {200e-9};
    distribution_opt["density_variable"] = {"density_from_H"};
    // distribution_opt["exponent"] = {2.5};
    // distribution_opt["bunch_length_fit"] = {"gauss"};


    // map<string, f_vector_t> extra_volt_dict;
    // extra_volt_dict["time_array"] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    // extra_volt_dict["voltage_array"] = {112, 2, 34, 41, 5, 6, 71, 8, 9, 10};

    auto ret = matched_from_distribution_density(
                   Beam, fullRing, distribution_opt,
                   FullRingAndRf::lowest_freq, totVolt);

    util::read_vector_from_file(v, params + "time_coord_low_res.txt");
    ASSERT_NEAR_LOOP(v, ret.time_coord_low_res, "time_coord_low_res", epsilon);

    util::read_vector_from_file(v, params + "line_density.txt");
    ASSERT_NEAR_LOOP(v, ret.line_density, "line_density", epsilon);

    delete totVolt;
    delete indVoltTime;
    delete resonator;
    delete Slice;
    delete long_tracker;
    delete fullRing;
}




TEST_F(testDistributions, matched_from_distribution_density3)
{
    auto RfP = Context::RfP;
    auto Beam = Context::Beam;
    auto Slice = new Slices(RfP, Beam, N_slices);
    auto epsilon = 1e-8;
    auto params = string(TEST_FILES "/Distributions/distribution_density3/");
    string datafiles = DEMO_FILES "/TC5_Wake_impedance/";


    f_vector_t v;
    util::read_vector_from_file(v, datafiles + "TC5_new_HQ_table.dat");

    f_vector_t R_shunt, f_res, Q_factor;

    R_shunt.reserve(v.size() / 3);
    f_res.reserve(v.size() / 3);
    Q_factor.reserve(v.size() / 3);

    for (uint i = 0; i < v.size(); i += 3) {
        f_res.push_back(v[i] * 1e9);
        Q_factor.push_back(v[i + 1]);
        R_shunt.push_back(v[i + 2] * 1e6);
    }

    auto resonator = new Resonators(R_shunt, f_res, Q_factor);
    auto indVoltTime = new InducedVoltageTime(Slice, {resonator});
    auto totVolt = new TotalInducedVoltage(Beam, Slice, {indVoltTime});

    auto long_tracker = new RingAndRfSection(RfP);
    auto fullRing = new FullRingAndRf({long_tracker});

    map<string, multi_t> distribution_opt;
    distribution_opt["type"] = {"waterbag"};
    distribution_opt["bunch_length"] = {200e-9};
    distribution_opt["density_variable"] = {"density_from_H"};
    // distribution_opt["exponent"] = {2.5};
    // distribution_opt["bunch_length_fit"] = {"gauss"};


    map<string, f_vector_t> extra_volt_dict;
    auto time_array = linspace(0.0, 1e-3, 10000);
    auto voltage_array = apply_f(time_array, [this](double x) {return sin(x) * V / 10.;});

    extra_volt_dict["time_array"] = time_array;
    extra_volt_dict["voltage_array"] = voltage_array;

    auto ret = matched_from_distribution_density(
                   Beam, fullRing, distribution_opt,
                   FullRingAndRf::lowest_freq, totVolt,
                   extra_volt_dict);

    util::read_vector_from_file(v, params + "time_coord_low_res.txt");
    ASSERT_NEAR_LOOP(v, ret.time_coord_low_res, "time_coord_low_res", epsilon);

    util::read_vector_from_file(v, params + "line_density.txt");
    ASSERT_NEAR_LOOP(v, ret.line_density, "line_density", epsilon);

    delete totVolt;
    delete indVoltTime;
    delete resonator;
    delete Slice;
    delete long_tracker;
    delete fullRing;
}

TEST_F(testDistributions, matched_from_distribution_density4)
{
    auto RfP = Context::RfP;
    auto Beam = Context::Beam;
    auto Slice = new Slices(RfP, Beam, N_slices);
    auto epsilon = 1e-8;
    auto params = string(TEST_FILES "/Distributions/distribution_density4/");
    string datafiles = DEMO_FILES "/TC5_Wake_impedance/";


    f_vector_t v;
    util::read_vector_from_file(v, datafiles + "TC5_new_HQ_table.dat");

    f_vector_t R_shunt, f_res, Q_factor;

    R_shunt.reserve(v.size() / 3);
    f_res.reserve(v.size() / 3);
    Q_factor.reserve(v.size() / 3);

    for (uint i = 0; i < v.size(); i += 3) {
        f_res.push_back(v[i] * 1e9);
        Q_factor.push_back(v[i + 1]);
        R_shunt.push_back(v[i + 2] * 1e6);
    }

    auto resonator = new Resonators(R_shunt, f_res, Q_factor);
    auto indVoltTime = new InducedVoltageTime(Slice, {resonator});
    auto totVolt = new TotalInducedVoltage(Beam, Slice, {indVoltTime});

    auto long_tracker = new RingAndRfSection(RfP);
    auto fullRing = new FullRingAndRf({long_tracker});

    map<string, multi_t> distribution_opt;
    distribution_opt["type"] = {"parabolic_line"};
    distribution_opt["bunch_length"] = {300e-9};
    distribution_opt["density_variable"] = {"density_from_J"};
    // distribution_opt["exponent"] = {2.5};
    distribution_opt["bunch_length_fit"] = {"gauss"};


    // map<string, f_vector_t> extra_volt_dict;
    // extra_volt_dict["time_array"] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    // extra_volt_dict["voltage_array"] = {112, 2, 34, 41, 5, 6, 71, 8, 9, 10};

    auto ret = matched_from_distribution_density(
                   Beam, fullRing, distribution_opt,
                   FullRingAndRf::lowest_freq, totVolt);

    util::read_vector_from_file(v, params + "time_coord_low_res.txt");
    ASSERT_NEAR_LOOP(v, ret.time_coord_low_res, "time_coord_low_res", epsilon);

    util::read_vector_from_file(v, params + "line_density.txt");
    ASSERT_NEAR_LOOP(v, ret.line_density, "line_density", epsilon);

    delete totVolt;
    delete indVoltTime;
    delete resonator;
    delete Slice;
    delete long_tracker;
    delete fullRing;
}

TEST_F(testDistributions, matched_from_distribution_density5)
{
    auto RfP = Context::RfP;
    auto Beam = Context::Beam;
    auto Slice = new Slices(RfP, Beam, N_slices);
    auto epsilon = 1e-8;
    auto params = string(TEST_FILES "/Distributions/distribution_density5/");
    f_vector_t v;

    auto long_tracker = new RingAndRfSection(RfP);
    auto fullRing = new FullRingAndRf({long_tracker});

    map<string, multi_t> distribution_opt;
    distribution_opt["type"] = {"gaussian"};
    distribution_opt["bunch_length"] = {300e-9};
    distribution_opt["density_variable"] = {"density_from_J"};
    distribution_opt["exponent"] = {0.7};
    distribution_opt["bunch_length_fit"] = {"end_to_end"};

    auto ret = matched_from_distribution_density(
                   Beam, fullRing, distribution_opt,
                   FullRingAndRf::lowest_freq);

    util::read_vector_from_file(v, params + "time_coord_low_res.txt");
    ASSERT_NEAR_LOOP(v, ret.time_coord_low_res, "time_coord_low_res", epsilon);

    util::read_vector_from_file(v, params + "line_density.txt");
    ASSERT_NEAR_LOOP(v, ret.line_density, "line_density", epsilon);
    delete Slice;
    delete long_tracker;
    delete fullRing;
}

TEST_F(testDistributions, matched_from_distribution_density6)
{
    auto RfP = Context::RfP;
    auto Beam = Context::Beam;
    auto Slice = new Slices(RfP, Beam, N_slices);
    auto epsilon = 1e-8;
    auto params = string(TEST_FILES "/Distributions/distribution_density6/");
    f_vector_t v;

    auto long_tracker = new RingAndRfSection(RfP);
    auto fullRing = new FullRingAndRf({long_tracker});

    map<string, multi_t> distribution_opt;
    distribution_opt["type"] = {"gaussian"};
    distribution_opt["bunch_length"] = {300e-9};
    distribution_opt["density_variable"] = {"density_from_J"};
    distribution_opt["exponent"] = {0.7};
    distribution_opt["bunch_length_fit"] = {"fwhm"};

    auto ret = matched_from_distribution_density(
                   Beam, fullRing, distribution_opt,
                   FullRingAndRf::lowest_freq);

    util::read_vector_from_file(v, params + "time_coord_low_res.txt");
    ASSERT_NEAR_LOOP(v, ret.time_coord_low_res, "time_coord_low_res", epsilon);

    util::read_vector_from_file(v, params + "line_density.txt");
    ASSERT_NEAR_LOOP(v, ret.line_density, "line_density", epsilon);

    delete Slice;
    delete long_tracker;
    delete fullRing;
}

TEST_F(testDistributions, matched_from_distribution_density_deathtest1)
{
    auto RfP = Context::RfP;
    auto Beam = Context::Beam;
    auto Slice = new Slices(RfP, Beam, N_slices);
    auto epsilon = 1e-8;

    auto long_tracker = new RingAndRfSection(RfP);
    auto fullRing = new FullRingAndRf({long_tracker});

    map<string, multi_t> distribution_opt;
    distribution_opt["type"] = {"gaussian"};
    distribution_opt["bunch_length"] = {300e-9};
    distribution_opt["density_variable"] = {""};
    distribution_opt["exponent"] = {0.7};

    ASSERT_DEATH(matched_from_distribution_density(Beam, fullRing, distribution_opt),
                 "\\w* was not recognized\n");


    delete Slice;
    delete long_tracker;
    delete fullRing;
}


class testBigaussian : public ::testing::Test {

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

        // Context::Beam = new Beams(N_p, N_b);

        auto GP = Context::GP;
        auto Beam = Context::Beam = new Beams(GP, N_p, N_b);

        auto RfP = Context::RfP = new RfParameters(GP, n_sections, hVec,
                voltageVec, dphiVec);

    }

    virtual void TearDown()
    {
        // Code here will be called immediately after each test
        // (right before the destructor).
        delete Context::GP;
        delete Context::Beam;
        delete Context::RfP;
    }

};

class testBigaussian2 : public ::testing::Test {

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

    const int N_t = 500;    // Number of turns to track
    const int N_p = 500000; // Macro-particles
    const int N_slices = 10;

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

    }

    virtual void TearDown()
    {
        // Code here will be called immediately after each test
        // (right before the destructor).
        delete Context::GP;
        delete Context::Beam;
        delete Context::RfP;
    }

};

TEST_F(testBigaussian, sigmas1)
{
    auto GP = Context::GP;
    auto Beam = Context::Beam;
    auto RfP = Context::RfP;
    string params = TEST_FILES "/Distributions/Bigaussian/sigmas1/";
    auto epsilon = 1e-8;

    longitudinal_bigaussian(GP, RfP, Beam, tau_0 / 4 , 0, -1, false);

    f_vector_t v;

    util::read_vector_from_file(v, params + "sigma_dE.txt");
    auto ref = v[0];
    auto real = Context::Beam->sigma_dE;
    ASSERT_NEAR(ref, real, epsilon * max(abs(ref), abs(real)));

    util::read_vector_from_file(v, params + "sigma_dt.txt");
    ref = v[0];
    real = Context::Beam->sigma_dt;
    ASSERT_NEAR(ref, real, epsilon * max(abs(ref), abs(real)));
}

TEST_F(testBigaussian, dE1)
{
    auto epsilon = 1e-8;

    string params = TEST_FILES "/Distributions/Bigaussian/dE1/";

    auto GP = Context::GP;
    auto Beam = Context::Beam;
    auto RfP = Context::RfP;

    longitudinal_bigaussian(GP, RfP, Beam, 0.1 * tau_0, 0, -1, false);

    f_vector_t v;
    util::read_vector_from_file(v, params + "dE.txt");
    for (uint i = 0; i < v.size(); ++i) {
        auto ref = v[i];
        auto real = Context::Beam->dE[i];
        ASSERT_NEAR(ref, real, epsilon * max(abs(ref), abs(real)))
                << "Testing of dE failed on i " << i << '\n';

    }
}

TEST_F(testBigaussian, dE2)
{
    auto epsilon = 1e-8;

    string params = TEST_FILES "/Distributions/Bigaussian/dE2/";
    auto GP = Context::GP;
    auto Beam = Context::Beam;
    auto RfP = Context::RfP;


    longitudinal_bigaussian(GP, RfP, Beam, tau_0, 1e3, -1, false);

    f_vector_t v;
    util::read_vector_from_file(v, params + "dE.txt");
    for (uint i = 0; i < v.size(); ++i) {
        auto ref = v[i];
        auto real = Context::Beam->dE[i];
        ASSERT_NEAR(ref, real, epsilon * max(abs(ref), abs(real)))
                << "Testing of dE failed on i " << i << '\n';

    }
}

TEST_F(testBigaussian, dt1)
{
    auto GP = Context::GP;
    auto Beam = Context::Beam;
    auto RfP = Context::RfP;

    auto epsilon = 1e-8;

    string params = TEST_FILES "/Distributions/Bigaussian/dt1/";

    longitudinal_bigaussian(GP, RfP, Beam, tau_0, 1e3, -1, false);

    f_vector_t v;
    util::read_vector_from_file(v, params + "dt.txt");
    for (uint i = 0; i < v.size(); ++i) {
        double ref = v[i];
        double real = Beam->dt[i];
        ASSERT_NEAR(ref, real, epsilon * max(abs(ref), abs(real)))
                << "Testing of dt failed on i " << i << '\n';

    }
}

TEST_F(testBigaussian2, dE3)
{
    auto epsilon = 1e-5;
    auto GP = Context::GP;
    auto Beam = Context::Beam;
    auto RfP = Context::RfP;


    string params = TEST_FILES "/Distributions/Bigaussian/dE3/";

    longitudinal_bigaussian(GP, RfP, Beam, 1e-9, 5e6, -1, false);

    f_vector_t v;

    util::read_vector_from_file(v, params + "dE_mean.txt");
    auto ref = v[0];
    auto real = mymath::mean(Beam->dE.data(), Beam->dE.size());
    ASSERT_NEAR(ref, real, epsilon * max(abs(ref), abs(real)))
            << "Testing of Beam->dE_mean failed\n";

    util::read_vector_from_file(v, params + "dE_std.txt");
    ref = v[0];
    real = mymath::standard_deviation(Beam->dE.data(), Beam->dE.size());
    ASSERT_NEAR(ref, real, epsilon * max(abs(ref), abs(real)))
            << "Testing of Beam->dE_std failed\n";
}

TEST_F(testBigaussian2, dt2)
{
    auto epsilon = 1e-5;
    auto GP = Context::GP;
    auto Beam = Context::Beam;
    auto RfP = Context::RfP;


    string params = TEST_FILES "/Distributions/Bigaussian/dt2/";

    longitudinal_bigaussian(GP, RfP, Beam, 1e-7, 5e4, -1, false);

    f_vector_t v;

    util::read_vector_from_file(v, params + "dt_mean.txt");
    auto ref = v[0];
    auto real = mymath::mean(Beam->dt.data(), Beam->dt.size());
    ASSERT_NEAR(ref, real, epsilon * max(abs(ref), abs(real)))
            << "Testing of Beam->dt_mean failed\n";

    util::read_vector_from_file(v, params + "dt_std.txt");
    ref = v[0];
    real = mymath::standard_deviation(Beam->dt.data(), Beam->dt.size());
    ASSERT_NEAR(ref, real, epsilon * max(abs(ref), abs(real)))
            << "Testing of Beam->dt_std failed\n";
}



int main(int ac, char *av[])
{
    python::initialize();
    ::testing::InitGoogleTest(&ac, av);
    auto ret = RUN_ALL_TESTS();
    python::finalize();
    return ret;
}

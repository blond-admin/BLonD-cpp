#include <blond/globals.h>
#include <blond/utilities.h>
#include <blond/math_functions.h>
#include <omp.h>
#include <stdio.h>
#include <blond/input_parameters/GeneralParameters.h>
#include <blond/input_parameters/RfParameters.h>
#include <blond/beams/Beams.h>
#include <blond/beams/Slices.h>
#include <blond/beams/Distributions.h>
#include <blond/trackers/Tracker.h>
#include <blond/impedances/Intensity.h>
#include <gtest/gtest.h>
#include <complex>


const std::string datafiles =
        "../demos/input_files/TC5_Wake_impedance/";

// Simulation parameters --------------------------------------------------------
// Bunch parameters
const int N_b = (int) 1e10;                          // Intensity
const ftype tau_0 = 2e-9;                       // Initial bunch length, 4 sigma [s]
// const particle_type particle = proton;
// Machine and RF parameters
const ftype C = 6911.56;                        // Machine circumference [m]
const ftype p_i = 25.92e9;                      // Synchronous momentum [eV/c]
//const ftype p_f = 460.005e9;                  // Synchronous momentum, final
const long h = 4620;                            // Harmonic number
const ftype V = 0.9e6;                          // RF voltage [V]
const ftype dphi = 0;                           // Phase modulation/offset
const ftype gamma_t = 1 / std::sqrt(0.00192);   // Transition gamma
const ftype alpha = 1.0 / gamma_t / gamma_t;    // First order mom. comp. factor
const int alpha_order = 1;
const int n_sections = 1;
// Tracking details

int N_t = 2;    // Number of turns to track
int N_p = 5000000;         // Macro-particles

int n_threads = 1;
int N_slices = 1 << 8; // = (2^8)

GeneralParameters *GP;
Beams *Beam;
Slices *Slice;
RfParameters *RfP;
//RingAndRfSection *long_tracker;
Resonators *resonator;

class testInputTableIntensity : public ::testing::Test {

protected:

    virtual void SetUp() {

        omp_set_num_threads(n_threads);

        f_vector_t momentum(N_t + 1);
        std::fill_n(momentum.begin(), N_t + 1, p_i);

        ftype *alpha_array = new ftype[(alpha_order + 1) * n_sections];
        std::fill_n(alpha_array, (alpha_order + 1) * n_sections, alpha);

        ftype *C_array = new ftype[n_sections];
        std::fill_n(C_array, n_sections, C);

        ftype *h_array = new ftype[n_sections * (N_t + 1)];
        std::fill_n(h_array, (N_t + 1) * n_sections, h);

        ftype *V_array = new ftype[n_sections * (N_t + 1)];
        std::fill_n(V_array, (N_t + 1) * n_sections, V);

        ftype *dphi_array = new ftype[n_sections * (N_t + 1)];
        std::fill_n(dphi_array, (N_t + 1) * n_sections, dphi);

        GP = new GeneralParameters(N_t, C_array, alpha_array, alpha_order, momentum.data(),
                                   proton);

        Beam = new Beams(N_p, N_b);

        RfP = new RfParameters(n_sections, h_array, V_array, dphi_array);

        //RingAndRfSection *long_tracker = new RingAndRfSection();

        longitudinal_bigaussian(tau_0 / 4, 0, 1, false);

        Slice = new Slices(N_slices, 0, 0, 2 * constant::pi, rad);
        //util::dump(Slice->bin_centers, 10, "bin_centers\n");

        std::vector<ftype> v;
        util::read_vector_from_file(v, datafiles +
                                       "TC5_new_HQ_table.dat");
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
        //printf("TEST_F Resonator Sizes %lu %lu\n", resonator->fRS.size(),
        //       resonator->fQ.size());



    }


    virtual void TearDown() {
        // Code here will be called immediately after each test
        // (right before the destructor).
        delete GP;
        delete Beam;
        delete RfP;
        delete Slice;
        delete resonator;
        //delete long_tracker;
    }


};


TEST_F(testInputTableIntensity, wake_calc
)
{

std::vector<ftype> timeArray;
timeArray.
reserve(N_slices);
for (
int i = 0;
i<N_slices;
++i) {
timeArray.
push_back(Slice
->bin_centers[i] - Slice->bin_centers[0]);
}
resonator->
wake_calc(timeArray);
std::vector<ftype> v1;
InputTable *inputTable = new InputTable(timeArray,
                                        resonator->fWake, v1);


for (
uint i = 0;
i<timeArray.

size();

++i) {
timeArray[i] = 1.1 * Slice->bin_centers[i] - Slice->bin_centers[0];
}
inputTable->
wake_calc(timeArray);

std::string params = "../unit-tests/references/Impedances/Intensity/InputTable/";

std::vector<ftype> v;
util::read_vector_from_file(v, params
+ "Wake.txt");

ASSERT_EQ(v
.

size(), inputTable

->fWake.

size()

);

ftype epsilon = 1e-8;

for (
unsigned int i = 0;
i<v.

size();

++i) {
ftype ref = v[i];
ftype real = inputTable->fWake[i];
ASSERT_NEAR(ref, real, epsilon
*
std::max(fabs(ref), fabs(real)
))
<< "Testing of fWake failed on i "
<< i <<
std::endl;
}
v.

clear();

delete
inputTable;
}


TEST_F(testInputTableIntensity, imped_calc
)
{

std::vector<ftype> timeArray;
timeArray.
reserve(N_slices);
for (
int i = 0;
i<N_slices;
++i) {
timeArray.
push_back(Slice
->bin_centers[i] - Slice->bin_centers[0]);
}
resonator->
wake_calc(timeArray);

std::transform(timeArray
.

begin(), timeArray

.

end(), timeArray

.

begin(),

[](
ftype a
) {
return a * 1e10;});
resonator->
imped_calc(timeArray);

std::vector<ftype> Re;
Re.
resize(resonator
->fImpedance.

size()

);
std::vector<ftype> Im;
Im.
resize(resonator
->fImpedance.

size()

);

std::transform(resonator
->fImpedance.

begin(), resonator

->fImpedance.

end(),
        Re

.

begin(),

[](
complex_t a
) {
return a.

real();

});
std::transform(resonator
->fImpedance.

begin(), resonator

->fImpedance.

end(),
        Im

.

begin(),

[](
complex_t a
) {
return a.

imag();

});

InputTable *inputTable = new InputTable(timeArray, Re, Im);

inputTable->
imped_calc(timeArray);

std::string params = "../unit-tests/references/Impedances/Intensity/InputTable/";

std::vector<ftype> v;
util::read_vector_from_file(v, params
+ "Impedance.txt");

ASSERT_EQ(v
.

size(), inputTable

->fImpedance.

size()

);

ftype epsilon = 1e-6;

for (
unsigned int i = 0;
i<v.

size();

++i) {
ftype ref = v[i];
ftype real = std::abs(inputTable->fImpedance[i]);
ASSERT_NEAR(ref, real, epsilon
*
std::max(fabs(ref), fabs(real)
))
<< "Testing of fImpedance failed on i "
<< i <<
std::endl;
}
v.

clear();

delete
inputTable;

}


int main(int ac, char *av[]) {
    ::testing::InitGoogleTest(&ac, av);
    return RUN_ALL_TESTS();
}
#include <iostream>
#include <string>
#include <list>

#include <unistd.h>

#include <gtest/gtest.h>
#include <blond/math_functions.h>
#include <blond/utilities.h>
#include <blond/beams/Distributions.h>
#include <blond/input_parameters/GeneralParameters.h>
#include <blond/trackers/Tracker.h>
#include <blond/llrf/PhaseLoop.h>

const ftype epsilon = 1e-7;
const std::string params = "../unit-tests/references/PL/PL_params/";

GeneralParameters *GP;
Beams *Beam;
RfParameters *RfP;
Slices *Slice;
LHC *PL;
RingAndRfSection *long_tracker;
int n_threads = 1;


class testPL : public ::testing::Test {

protected:

    //const ftype tau_0 = 0.4e-9;          // Initial bunch length, 4 sigma [s]


    virtual void SetUp() {
        //printf("ok here\n");

        std::vector<ftype> v;
        util::read_vector_from_file(v, datafiles + "LHC_momentum_programme");

        // optional
        v.erase(v.begin(), v.begin() + from_line);

        int remaining = N_t + 1 - v.size();
        for (int i = 0; i < remaining; ++i) {
            v.push_back(6.5e12);
        }
        assert((int) v.size() == N_t + 1);
        ftype *ps = &v[0];  //new ftype[v.size()];


        ftype *V_array = new ftype[N_t + 1];
        mymath::linspace(V_array, 6e6, 10e6, 13563374, 13e6);
        std::fill_n(&V_array[563374], 436627, 10e6);

        // Define general parameters

        ftype *alpha_array = new ftype[(alpha_order + 1) * n_sections];
        std::fill_n(alpha_array, (alpha_order + 1) * n_sections, alpha);

        ftype *C_array = new ftype[n_sections];
        std::fill_n(C_array, n_sections, C);

        GP = new GeneralParameters(N_t, C_array, alpha_array, alpha_order, ps,
                                   proton);

        // Define rf_params
        ftype *dphi_array = new ftype[n_sections * (N_t + 1)];
        std::fill_n(dphi_array, (N_t + 1) * n_sections, dphi);

        ftype *h_array = new ftype[n_sections * (N_t + 1)];
        std::fill_n(h_array, (N_t + 1) * n_sections, h);

        RfP = new RfParameters(n_sections, h_array, V_array, dphi_array);

        // Define beam and distribution: Load matched, filamented distribution
        Beam = new Beams(N_p, N_b);
        std::vector<ftype> v2;
        util::read_vector_from_file(v2, datafiles + "coords_13000001.dat");

        int k = 0;
        for (unsigned int i = 0; i < v2.size(); i += 3) {
            Beam->dt[k] = v2[i] * 1e-9; // [s]
            Beam->dE[k] = v2[i + 1] * 1e6; // [eV]
            Beam->id[k] = v2[i + 2];
            k++;
        }

        Slice = new Slices(N_slices, 0, -0.5e-9, 3e-9);

        // Define phase loop and frequency loop gain
        ftype PL_gain = 1 / (5 * GP->t_rev[0]);
        ftype SL_gain = PL_gain / 10;

        ftype *PL_gain_array = new ftype[N_t + 1];
        std::fill_n(PL_gain_array, N_t + 1, PL_gain);

        PL = new LHC(PL_gain_array, SL_gain);

        long_tracker = new RingAndRfSection(simple, PL);

    }


    virtual void TearDown() {
        // Code here will be called immediately after each test
        // (right before the destructor).
        delete GP;
        delete Beam;
        delete RfP;
        delete Slice;
        delete PL;
        delete long_tracker;
    }


private:


    // Machine and RF parameters
    const float C = 26658.883;        // Machine circumference [m]
    const int h = 35640;            // Harmonic number
    const float dphi = 0.;            // Phase modulation/offset
    const float gamma_t = 55.759505;  // Transition gamma
    const float alpha = 1. / gamma_t / gamma_t;     // First order mom. comp. factor

    // Tracking details
    const int N_p = 100000;         // Macro-particles
    const long N_b = 1e9;           // Intensity
    const int alpha_order = 1;
    const int n_sections = 1;
    const int N_t = 1000000;        // Number of turns to track; full ramp: 8700001
    const int bl_target = 1.25e-9;  // 4 sigma r.m.s. target bunch length in [s]

    const int N_slices = 151;

    const std::string datafiles =
            "../demos/input_files/LHC_restart/";

    const int from_line = 0;


};


TEST_F(testPL, lhc_a
)
{

std::vector<ftype> v;
util::read_vector_from_file(v, params
+ "lhc_a");
// Only check 1 out of 10 elements
// otherwise refernece file too big
ASSERT_EQ(v
.

size(), GP

->n_turns / 10 + 1);
for (
unsigned int i = 0;
i<v.

size();

++i) {
//printf("ok here \n");

ftype ref = v[i];
ftype real = PL->lhc_a[i * 10];
ASSERT_NEAR(ref, real, epsilon
*
std::max(fabs(ref), fabs(real)
));
}

//printf("ok here\n");
}

TEST_F(testPL, lhc_t
)
{

std::vector<ftype> v;
util::read_vector_from_file(v, params
+ "lhc_t");
// Only check 1 out of 10 elements
// otherwise refernece file too big
ASSERT_EQ(v
.

size(), GP

->n_turns / 10 + 1);
for (
unsigned int i = 0;
i<v.

size();

++i) {
//printf("%d\n", i);
ftype ref = v[i];
ftype real = PL->lhc_t[i * 10];
ASSERT_NEAR(ref, real, epsilon
*
std::max(fabs(ref), fabs(real)
));
}
}

TEST_F(testPL, phi_beam
)
{
Slice->

track();

PL->

beam_phase();

std::vector<ftype> v;
util::read_vector_from_file(v, params
+ "phi_beam");
ftype ref = v[0];
ftype real = PL->phi_beam;
ASSERT_NEAR(ref, real, epsilon
*
std::max(fabs(ref), fabs(real)
));
}

TEST_F(testPL, dphi
)
{
Slice->

track();

PL->

beam_phase();

PL->

phase_difference();

std::vector<ftype> v;
util::read_vector_from_file(v, params
+ "dphi");
ftype ref = v[0];
ftype real = PL->dphi;
ASSERT_NEAR(ref, real, epsilon
*
std::max(fabs(ref), fabs(real)
));
}

TEST_F(testPL, domega_RF
)
{
Slice->

track();

PL->

beam_phase();

PL->

phase_difference();

long_tracker->

track();

std::vector<ftype> v;
util::read_vector_from_file(v, params
+ "domega_RF");
ftype ref = v[0];
ftype real = PL->domega_RF;
ASSERT_NEAR(ref, real, epsilon
*
std::max(fabs(ref), fabs(real)
));
}

TEST_F(testPL, lhc_y
)
{
Slice->

track();

PL->

beam_phase();

PL->

phase_difference();

long_tracker->

track();

std::vector<ftype> v;
util::read_vector_from_file(v, params
+ "lhc_y");
ftype ref = v[0];
ftype real = PL->lhc_y;
ASSERT_NEAR(ref, real, epsilon
*
std::max(fabs(ref), fabs(real)
));
}

TEST_F(testPL, omega_RF
)
{
Slice->

track();

PL->

beam_phase();

PL->

phase_difference();

long_tracker->

track();

//RfP->counter++;
int counter = RfP->counter;

std::vector<ftype> v;
util::read_vector_from_file(v, params
+ "omega_RF");
ftype ref = v[0];
ftype real = RfP->omega_RF[counter];
ASSERT_NEAR(ref, real, epsilon
*
std::max(fabs(ref), fabs(real)
));
}

TEST_F(testPL, dphi_RF
)
{
Slice->

track();

PL->

beam_phase();

PL->

phase_difference();

long_tracker->

track();
//RfP->counter++;

std::vector<ftype> v;
util::read_vector_from_file(v, params
+ "dphi_RF");
ftype ref = v[0];
ftype real = RfP->dphi_RF[0];
ASSERT_NEAR(ref, real, epsilon
*
std::max(fabs(ref), fabs(real)
));
}

TEST_F(testPL, phi_RF
)
{
Slice->

track();

PL->

beam_phase();

PL->

phase_difference();

long_tracker->

track();
//RfP->counter++;

int counter = RfP->counter;

std::vector<ftype> v;
util::read_vector_from_file(v, params
+ "phi_RF");
ftype ref = v[0];
ftype real = RfP->phi_RF[counter];
ASSERT_NEAR(ref, real, epsilon
*
std::max(fabs(ref), fabs(real)
));
}


int main(int ac, char *av[]) {
    ::testing::InitGoogleTest(&ac, av);
    return RUN_ALL_TESTS();
}
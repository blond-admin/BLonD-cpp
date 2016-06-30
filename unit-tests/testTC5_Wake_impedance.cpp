#include "globals.h"
#include "utilities.h"
#include "math_functions.h"
#include <omp.h>
#include <stdio.h>
#include "../input_parameters/GeneralParameters.h"
#include "../input_parameters/RfParameters.h"
#include "../beams/Beams.h"
#include "../beams/Slices.h"
#include "../beams/Distributions.h"
#include "../trackers/Tracker.h"
#include "../impedances/InducedVoltage.h"
#include <gtest/gtest.h>
#include <complex>


const std::string datafiles =
   "../tests/input_files/TC5_Wake_impedance/";

// Simulation parameters --------------------------------------------------------
// Bunch parameters
const long int N_b = (long int) 1e10;                          // Intensity
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

unsigned N_t = 100;    // Number of turns to track
unsigned N_p = 500;         // Macro-particles

int n_threads = 1;
unsigned N_slices = 1 << 8; // = (2^8)

GeneralParameters *GP;
Beams *Beam;
Slices *Slice;
RfParameters *RfP;
RingAndRfSection *long_tracker;
Resonators *resonator;

class testTC5 : public ::testing::Test {

protected:

   virtual void SetUp()
   {

      omp_set_num_threads(n_threads);

      f_vector_2d_t momentumVec(n_sections, f_vector_t(N_t + 1, p_i));

      //f_vector_2d_t alphaVec(alpha_order + 1, f_vector_t(n_sections, alpha));
      f_vector_2d_t alphaVec(n_sections, f_vector_t(alpha_order +1, alpha));

      f_vector_t CVec(n_sections, C);

      f_vector_2d_t hVec(n_sections , f_vector_t(N_t + 1, h));

      f_vector_2d_t voltageVec(n_sections , f_vector_t(N_t + 1, V));

      f_vector_2d_t dphiVec(n_sections , f_vector_t(N_t + 1, dphi));

      GP = new GeneralParameters(N_t, CVec, alphaVec, alpha_order,
                                 momentumVec, proton);

      Beam = new Beams(N_p, N_b);

      RfP = new RfParameters(n_sections, hVec, voltageVec, dphiVec);

      long_tracker = new RingAndRfSection();

      longitudinal_bigaussian(tau_0 / 4, 0, -1, false);

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

   }


   virtual void TearDown()
   {
      // Code here will be called immediately after each test
      // (right before the destructor).
      delete GP;
      delete Beam;
      delete RfP;
      delete Slice;
      delete long_tracker;
      delete resonator;
   }


};



TEST_F(testTC5, timeTrack)
{

   std::vector<Intensity *> wakeSourceList({resonator});
   InducedVoltageTime *indVoltTime = new InducedVoltageTime(wakeSourceList);
   std::vector<InducedVoltage *> indVoltList({indVoltTime});

   TotalInducedVoltage *totVol = new TotalInducedVoltage(indVoltList);


   for (unsigned i = 0; i < N_t; ++i) {
      totVol->track();
      long_tracker->track();
      Slice->track();
   }

   auto params = std::string("../unit-tests/references/")
                 + "TC5_final/time/";

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
            << "Testing of Beam->dE failed on i "
            << i << std::endl;
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
            << "Testing of Beam->dt failed on i "
            << i << std::endl;
   }

   v.clear();
   res.clear();
   util::read_vector_from_file(v, params + "n_macroparticles.txt");

   res = f_vector_t(Slice->n_macroparticles.begin(), Slice->n_macroparticles.end());
   ASSERT_EQ(v.size(), res.size());

   epsilon = 1e-8;
   // warning checking only the first 100 elems
   for (unsigned int i = 0; i < v.size(); ++i) {
      ftype ref = v[i];
      ftype real = Slice->n_macroparticles[i];

      ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
            << "Testing of Slice->n_macroparticles failed on i "
            << i << std::endl;
   }

   delete indVoltTime;
   delete totVol;

}



TEST_F(testTC5, freqTrack)
{
   std::vector<Intensity *> ImpSourceList({resonator});
   auto indVoltFreq = new InducedVoltageFreq(ImpSourceList, 1e5);
   std::vector<InducedVoltage *> indVoltList({indVoltFreq});

   TotalInducedVoltage *totVol = new TotalInducedVoltage(indVoltList);

   for (unsigned i = 0; i < N_t; ++i) {
      totVol->track();
      long_tracker->track();
      Slice->track();
   }

   auto params = std::string("../unit-tests/references/")
                 + "TC5_final/freq/";

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
            << "Testing of Beam->dE failed on i "
            << i << std::endl;
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
            << "Testing of Beam->dt failed on i "
            << i << std::endl;
   }

   v.clear();
   res.clear();
   util::read_vector_from_file(v, params + "n_macroparticles.txt");

   res = f_vector_t(Slice->n_macroparticles.begin(), Slice->n_macroparticles.end());
   ASSERT_EQ(v.size(), res.size());

   epsilon = 1e-8;
   for (unsigned int i = 0; i < v.size(); ++i) {
      ftype ref = v[i];
      ftype real = Slice->n_macroparticles[i];

      ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
            << "Testing of Slice->n_macroparticles failed on i "
            << i << std::endl;
   }

   delete indVoltFreq;
   delete totVol;

}





int main(int ac, char *av[])
{
   ::testing::InitGoogleTest(&ac, av);
   return RUN_ALL_TESTS();
}

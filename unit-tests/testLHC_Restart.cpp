#include <iostream>
#include <string>
#include <list>

#include <unistd.h>

#include <gtest/gtest.h>
#include "math_functions.h"
#include "utilities.h"
#include "../beams/Distributions.h"
#include "../input_parameters/GeneralParameters.h"
#include "../trackers/Tracker.h"
#include "../llrf/PhaseLoop.h"
#include <omp.h>


//const ftype epsilon = 1e-3;
const std::string params = "../unit-tests/references/PL/LHC_restart_params/";

GeneralParameters *GP;
Beams *Beam;
RfParameters *RfP;
Slices *Slice;
LHC *PL;
RingAndRfSection *long_tracker;
int n_threads = 1;


class testLHC_Restart : public ::testing::Test {

protected:
   const int N_p = 100000;         // Macro-particles
   //const ftype tau_0 = 0.4e-9;          // Initial bunch length, 4 sigma [s]
   const int N_t = 1000000;        // Number of turns to track; full ramp: 8700001

   const int N_slices = 151;

   virtual void SetUp()
   {
      //printf("ok here\n");
      omp_set_num_threads(n_threads);

      std::vector < ftype > v;
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
      std::vector < ftype > v2;
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


   virtual void TearDown()
   {
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

   const long N_b = 1e9;           // Intensity
   const int alpha_order = 1;
   const int n_sections = 1;
   const int bl_target = 1.25e-9;  // 4 sigma r.m.s. target bunch length in [s]


   const std::string datafiles =
      "../tests/input_files/LHC_restart/";

   const int from_line = 0;


};



TEST_F(testLHC_Restart, dphi_RF_and_dphi)
{



   std::vector<ftype> real1, real2;


   for (int i = 0; i < 1000; ++i) {


      if (RfP->counter < 570000)
         PL->reference = 0.5236;
      else
         PL->reference = 1.0472;

      Slice->track();

      long_tracker->track();

      //RfP->counter++;
      //printf("   Beam energy %.6e eV\n", GP->energy[0]);
      //printf("   RF phase %.6e rad\n", RfP->dphi_RF[0]);
      //printf("   PL phase correction %.6e rad\n", PL->dphi);
      //ftype ref = v1[i];
      //ftype real = RfP->dphi_RF[RfP->counter];
      real1.push_back(RfP->dphi_RF[0]);
      real2.push_back(PL->dphi);

   }

   std::vector<ftype> v;


   util::read_vector_from_file(v, params + "dphi_RF[0]");

   ASSERT_EQ(v.size(), real1.size());
   ftype epsilon = 1e-6;
   for (unsigned int i = 0; i < v.size(); ++i) {
      //printf("%d\n", i);
      ftype ref = v[i];
      ftype real = real1[i];
      ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
            << "Testing of dphi_RF failed on i "
            << i << std::endl;
   }

   epsilon = 1e-3;
   v.clear();
   util::read_vector_from_file(v, params + "dphi");
   ASSERT_EQ(v.size(), real2.size());
   for (unsigned int i = 0; i < v.size(); ++i) {
      //printf("%d\n", i);
      ftype ref = v[i];
      ftype real = real2[i];
      ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
            << "Testing of dphi failed on i "
            << i << std::endl;
   }


   epsilon = 5e-1;
   v.clear();
   util::read_vector_from_file(v, params + "dE");
   ASSERT_EQ(v.size(), Beam->n_macroparticles);
   for (unsigned int i = 0; i < v.size(); ++i) {
      //printf("%d\n", i);
      ftype ref = v[i];
      ftype real = Beam->dE[i];
      ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
            << "Testing of dE failed on i "
            << i << std::endl;
   }

   epsilon = 1e-6;
   v.clear();
   util::read_vector_from_file(v, params + "dt");
   ASSERT_EQ(v.size(), Beam->n_macroparticles);
   for (unsigned int i = 0; i < v.size(); ++i) {
      //printf("%d\n", i);
      ftype ref = v[i];
      ftype real = Beam->dt[i];
      ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
            << "Testing of dt failed on i "
            << i << std::endl;
   }

   epsilon = 1e-8;
   v.clear();
   util::read_vector_from_file(v, params + "n_macroparticles");
   ASSERT_EQ(v.size(), Slice->n_slices);
   for (unsigned int i = 0; i < v.size(); ++i) {
      //printf("%d\n", i);
      ftype ref = v[i];
      ftype real = Slice->n_macroparticles[i];
      ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
            << "Testing of n_macroparticles failed on i "
            << i << std::endl;
   }

}



int main(int ac, char *av[])
{
   ::testing::InitGoogleTest(&ac, av);
   return RUN_ALL_TESTS();
}
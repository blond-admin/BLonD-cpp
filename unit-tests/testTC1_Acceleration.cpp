#include <iostream>

#include <omp.h>
#include <gtest/gtest.h>
#include <blond/math_functions.h>
#include <blond/utilities.h>
#include <blond/beams/Distributions.h>
#include <blond/trackers/Tracker.h>
using namespace blond;

const ftype epsilon = 1e-8;
const std::string params = "../unit-tests/references/TC1_final/TC1_final_params/";

class testTC1 : public ::testing::Test {

protected:
   const long N_b = 1e9;           // Intensity
   const ftype tau_0 = 0.4e-9;          // Initial bunch length, 4 sigma [s]

   const int N_t = 2000;    // Number of turns to track
   const int N_p = 100;         // Macro-particles
   const int N_slices = 10;

   virtual void SetUp()
   {
      ftype *momentum = new ftype[N_t + 1];
      mymath::linspace(momentum, p_i, p_f, N_t + 1);

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

      context.GP = new GeneralParameters(N_t, C_array, alpha_array, alpha_order, momentum,
                                 proton);

      context.Beam = new Beams(N_p, N_b);

      context.RfP = new RfParameters(n_sections, h_array, V_array, dphi_array);

      longitudinal_bigaussian(tau_0 / 4, 0, 1, false);

      context.Slice = new Slices(N_slices);

   }


   virtual void TearDown()
   {
      // Code here will be called immediately after each test
      // (right before the destructor).
      delete context.GP;
      delete context.Beam;
      delete context.RfP;
      delete context.Slice;
   }


private:

   // Machine and RF parameters
   const ftype C = 26658.883;          // Machine circumference [m]
   const long p_i = 450e9;          // Synchronous momentum [eV/c]
   const ftype p_f = 460.005e9;          // Synchronous momentum, final
   const long h = 35640;          // Harmonic number
   const ftype V = 6e6;          // RF voltage [V]
   const ftype dphi = 0;          // Phase modulation/offset
   const ftype gamma_t = 55.759505;          // Transition gamma
   const ftype alpha = 1.0 / gamma_t / gamma_t; // First order mom. comp. factor
   const int alpha_order = 1;
   const int n_sections = 1;
   // Tracking details



};


TEST_F(testTC1, phaseSpace)
{
   auto Beam = context.Beam;
   omp_set_num_threads(context.n_threads);
   RingAndRfSection *long_tracker = new RingAndRfSection();
   for (int i = 0; i < N_t; ++i) {
      long_tracker->track();
      context.Slice->track();
      context.RfP->counter++;
      //beam->losses_longitudinal_cut(beam->dt, 0, 2.5e-9, beam->id);
   }

   std::vector<ftype> v;
   util::read_vector_from_file(v, params + "dE");
   for (unsigned int i = 0; i < v.size(); ++i) {
      ftype ref_dE = v[i];
      ftype real_dE = Beam->dE[i];
      ASSERT_NEAR(ref_dE, real_dE, epsilon * std::max(fabs(ref_dE), fabs(real_dE)));
   }
   v.clear();
   util::read_vector_from_file(v, params + "dt");
   for (unsigned int i = 0; i < v.size(); ++i) {
      ftype ref_dt = v[i];
      ftype real_dt = Beam->dt[i];
      ASSERT_NEAR(ref_dt, real_dt, epsilon * std::max(fabs(ref_dt), fabs(real_dt)));
   }

   delete long_tracker;



}



int main(int ac, char *av[])
{
   ::testing::InitGoogleTest(&ac, av);
   return RUN_ALL_TESTS();
}
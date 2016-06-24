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
#include "constants.h"

const ftype epsilon = 1e-8;
const std::string track_params = "../unit-tests/references/Slices/Slices_track_params/";
const std::string set_cuts_params = "../unit-tests/references/Slices/Slices_set_cuts_params/";
const std::string sort_particles_params = "../unit-tests/references/Slices/Slices_sort_particles_params/";

GeneralParameters *GP;
Beams *Beam;
RfParameters *RfP;
Slices *Slice;
int n_threads = 1;


class testSlices : public ::testing::Test {

protected:
   const long N_b = 1e9;           // Intensity
   const ftype tau_0 = 0.4e-9;          // Initial bunch length, 4 sigma [s]

   virtual void SetUp()
   {
      f_vector_2d_t momentumVec(n_sections, f_vector_t(N_t + 1));
      for (auto &v : momentumVec)
         mymath::linspace(v.data(), p_i, p_f, N_t + 1);

      f_vector_2d_t alphaVec(n_sections , f_vector_t(alpha_order+1, alpha));

      f_vector_t CVec(n_sections, C);

      f_vector_2d_t hVec(n_sections , f_vector_t(N_t + 1, h));

      f_vector_2d_t voltageVec(n_sections , f_vector_t(N_t + 1, V));

      f_vector_2d_t dphiVec(n_sections , f_vector_t(N_t + 1, dphi));


      GP = new GeneralParameters(N_t, CVec, alphaVec, alpha_order, momentumVec,
                                 proton);

      Beam = new Beams(N_p, N_b);

      RfP = new RfParameters(n_sections, hVec, voltageVec, dphiVec);


      longitudinal_bigaussian(tau_0 / 4, 0, 1, false);

      Slice = new Slices(N_slices);

   }


   virtual void TearDown()
   {
      // Code here will be called immediately after each test
      // (right before the destructor).
      delete GP;
      delete Beam;
      delete RfP;
      delete Slice;
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

   const int N_t = 2000;    // Number of turns to track
   const int N_p = 100;         // Macro-particles
   const int N_slices = 10;


};



TEST_F(testSlices, set_cuts_left)
{
   //putenv("FIXED_PARTICLES=1");

   std::vector<ftype> v;

   util::read_vector_from_file(v, set_cuts_params + "cut_left");
   ftype ref = v[0];
   ftype real = Slice->cut_left;
   ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)));
}

TEST_F(testSlices, set_cuts_right)
{
   //putenv("FIXED_PARTICLES=1");

   std::vector<ftype> v;

   util::read_vector_from_file(v, set_cuts_params + "cut_right");
   ftype ref = v[0];
   ftype real = Slice->cut_right;
   ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)));
}

TEST_F(testSlices, set_cuts_bin_centers)
{
   //putenv("FIXED_PARTICLES=1");

   std::vector<ftype> v;

   util::read_vector_from_file(v, set_cuts_params + "bin_centers");
   for (unsigned int i = 0; i < v.size(); ++i) {
      ftype ref = v[i];
      ftype real = Slice->bin_centers[i];
      ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)));
   }
}

// TODO test set_cuts a little further


TEST_F(testSlices, sort_particles_dE)
{

   std::vector<ftype> v;
   util::read_vector_from_file(v, sort_particles_params + "dE");
   for (unsigned int i = 0; i < v.size(); ++i) {
      ftype ref = v[i];
      ftype real = Beam->dE[i];
      ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)));
   }
}

TEST_F(testSlices, sort_particles_dt)
{

   std::vector<ftype> v;
   util::read_vector_from_file(v, sort_particles_params + "dt");
   for (unsigned int i = 0; i < v.size(); ++i) {
      ftype ref = v[i];
      ftype real = Beam->dt[i];
      ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)));
   }
}

TEST_F(testSlices, sort_particles_id)
{

   std::vector<ftype> v;
   util::read_vector_from_file(v, sort_particles_params + "id");
   for (unsigned int i = 0; i < v.size(); ++i) {
      ftype ref = v[i];
      ftype real = Beam->id[i];
      ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)));
   }
}

TEST_F(testSlices, n_macroparticles)
{

   //RingAndRfSection *long_tracker = new RingAndRfSection();
   //long_tracker->track(0, Beam->n_macroparticles);
   Slice->track();
   //util::dump(Slice->n_macroparticles, 100, "something\n");
   std::vector<ftype> v;
   util::read_vector_from_file(v, track_params + "n_macroparticles");
   for (unsigned int i = 0; i < v.size(); ++i) {
      ftype ref = v[i];
      ftype real = Slice->n_macroparticles[i];
      ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)));
   }
}



int main(int ac, char *av[])
{
   ::testing::InitGoogleTest(&ac, av);
   return RUN_ALL_TESTS();
}

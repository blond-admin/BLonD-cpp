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
const std::string track_params = "../unit-tests/references/Tracker/Tracker_track_params/";

GeneralParameters *GP;
Beams *Beam;
RfParameters *RfP;
Slices *Slice;
int n_threads = 1;


class testTracker : public ::testing::Test {

protected:
   const long N_b = 1e9;           // Intensity
   const ftype tau_0 = 0.4e-9;          // Initial bunch length, 4 sigma [s]

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

      GP = new GeneralParameters(N_t, C_array, alpha_array, alpha_order, momentum,
                                 proton);

      Beam = new Beams(N_p, N_b);

      RfP = new RfParameters(n_sections, h_array, V_array, dphi_array);

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


class testTrackerPeriodicity : public ::testing::Test {

protected:
   const long N_b = 1e9;           // Intensity
   const ftype tau_0 = 0.4e-9;          // Initial bunch length, 4 sigma [s]

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

      GP = new GeneralParameters(N_t, C_array, alpha_array, alpha_order, momentum,
                                 proton);

      Beam = new Beams(N_p, N_b);

      RfP = new RfParameters(n_sections, h_array, V_array, dphi_array);

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


TEST_F(testTracker, track_dE)
{

   RingAndRfSection *long_tracker = new RingAndRfSection();
   long_tracker->track();

   std::vector<ftype> v;
   util::read_vector_from_file(v, track_params + "dE");
   for (unsigned int i = 0; i < v.size(); ++i) {
      ftype ref = v[i];
      ftype real = Beam->dE[i];
      ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)));
   }
   delete long_tracker;
}

TEST_F(testTracker, track_dt)
{

   RingAndRfSection *long_tracker = new RingAndRfSection();

   long_tracker->track();

   std::vector<ftype> v;
   util::read_vector_from_file(v, track_params + "dt");
   for (unsigned int i = 0; i < v.size(); ++i) {
      ftype ref = v[i];
      ftype real = Beam->dt[i];
      ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)));
   }
   delete long_tracker;

}


TEST_F(testTrackerPeriodicity, kick)
{
   auto params = std::string("../unit-tests/references/Tracker/periodicity/kick/");
   RingAndRfSection *long_tracker = new RingAndRfSection(simple, NULL, NULL, true, 0.0);

   int_vector_t indices = {1, 2, 4, 8, 16, 32, 64};

   for (int i = 0; i < 100; ++i) {
      long_tracker->kick(indices, RfP->counter);
      //long_tracker->track();
   }
   //util::dump(Beam->dE.data(), 100, "Beam->dE\n");
   std::vector<ftype> v;
   util::read_vector_from_file(v, params + "beam_dE.txt");

   auto res = std::vector<ftype>(Beam->dE.data(), Beam->dE.data() + Beam->n_macroparticles);

   ASSERT_EQ(v.size(), res.size());
   auto epsilon = 1e-8;

   for (unsigned int i = 0; i < v.size(); ++i) {
      ftype ref = v[i];
      ftype real = res[i];
      ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
            << "Testing of Beam->dE failed on i "
            << i << std::endl;
   }

   delete long_tracker;

}


TEST_F(testTrackerPeriodicity, drift)
{
   auto params = std::string("../unit-tests/references/Tracker/periodicity/drift/");
   RingAndRfSection *long_tracker = new RingAndRfSection(simple, NULL, NULL, true, 0.0);

   int_vector_t indices = {0, 1, 2, 3, 4, 98, 99};

   for (int i = 0; i < 100; ++i) {
      long_tracker->drift(indices, RfP->counter);
      //long_tracker->track();
   }
   //util::dump(Beam->dE.data(), 100, "Beam->dE\n");
   std::vector<ftype> v;
   util::read_vector_from_file(v, params + "beam_dt.txt");

   auto res = std::vector<ftype>(Beam->dt.data(), Beam->dt.data() + Beam->n_macroparticles);

   ASSERT_EQ(v.size(), res.size());
   auto epsilon = 1e-8;

   for (unsigned int i = 0; i < v.size(); ++i) {
      ftype ref = v[i];
      ftype real = res[i];
      ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
            << "Testing of Beam->dt failed on i "
            << i << std::endl;
   }

   delete long_tracker;

}

TEST_F(testTrackerPeriodicity, set_periodicity)
{
   auto params = std::string("../unit-tests/references/Tracker/periodicity/set_periodicity/");
   RingAndRfSection *long_tracker = new RingAndRfSection(simple, NULL, NULL, true, 0.0);

   ftype mean = mymath::mean<ftype>(Beam->dt.data(), Beam->dt.size());

   GP->t_rev[RfP->counter + 1 ] = mean;

   long_tracker->set_periodicity();


   std::vector<ftype> v;
   util::read_vector_from_file(v, params + "indices_right_outside.txt");

   auto res = long_tracker->indices_right_outside;

   ASSERT_EQ(v.size(), res.size());
   auto epsilon = 1e-8;

   for (unsigned int i = 0; i < v.size(); ++i) {
      ftype ref = v[i];
      ftype real = res[i];
      ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
            << "Testing of Beam->indices_right_outside failed on i "
            << i << std::endl;
   }

   v.clear();
   res.clear();
   util::read_vector_from_file(v, params + "indices_inside_frame.txt");

   res = long_tracker->indices_inside_frame;

   ASSERT_EQ(v.size(), res.size());
   epsilon = 1e-8;

   for (unsigned int i = 0; i < v.size(); ++i) {
      ftype ref = v[i];
      ftype real = res[i];
      ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
            << "Testing of Beam->indices_inside_frame failed on i "
            << i << std::endl;
   }


   delete long_tracker;

}


TEST_F(testTrackerPeriodicity, track1)
{
   auto params = std::string("../unit-tests/references/Tracker/periodicity/track1/");
   RingAndRfSection *long_tracker = new RingAndRfSection(simple, NULL, NULL, true, 0.0);

   ftype mean = mymath::mean<ftype>(Beam->dt.data(), Beam->dt.size());

   for (uint i = 0; i < Beam->dt.size(); ++i) {
      if (i % 5 == 0)
         Beam->dt[i] += mean;
   }

   //util::dump(Beam->dt.data(), Beam->dt.size(), "Beam->dt before tracking\n");


   for (int i = 0; i < 100; ++i) {
      long_tracker->track();
   }
   std::vector<ftype> v;
   util::read_vector_from_file(v, params + "beam_dt.txt");

   //util::dump(Beam->dt.data(), Beam->dt.size(), "Beam->dt\n");
   std::vector<ftype> res(Beam->dt.data(), Beam->dt.data() + Beam->n_macroparticles);
   ASSERT_EQ(v.size(), res.size());
   ftype epsilon = 1e-5;

   for (unsigned int i = 0; i < v.size(); ++i) {
      ftype ref = v[i];
      ftype real = res[i];
      ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
            << "Testing of Beam->dt failed on i "
            << i << std::endl;
   }


   v.clear();
   res.clear();

   res = std::vector<ftype>(Beam->dE.data(), Beam->dE.data() + Beam->n_macroparticles);

   util::read_vector_from_file(v, params + "beam_dE.txt");

   ASSERT_EQ(v.size(), res.size());
   epsilon = 1e-5;

   for (unsigned int i = 0; i < v.size(); ++i) {
      ftype ref = v[i];
      ftype real = res[i];
      ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
            << "Testing of Beam->dE failed on i "
            << i << std::endl;
   }

   delete long_tracker;

}


TEST_F(testTrackerPeriodicity, track2)
{
   auto params = std::string("../unit-tests/references/Tracker/periodicity/track2/");
   RingAndRfSection *long_tracker = new RingAndRfSection(simple, NULL, NULL, true, 0.0);

   ftype mean = mymath::mean<ftype>(Beam->dt.data(), Beam->dt.size());

   for (uint i = 0; i < Beam->dt.size(); ++i) {
      if (i % 2 == 0)
         Beam->dt[i] = std::max( Beam->dt[i] - mean, 0.0);
      else if (i % 3 == 0)
         Beam->dt[i] += mean;
   }


   for (int i = 0; i < 100; ++i) {
      long_tracker->track();
   }
   std::vector<ftype> v;
   util::read_vector_from_file(v, params + "beam_dt.txt");

   //util::dump(Beam->dt.data(), Beam->dt.size(), "Beam->dt\n");
   std::vector<ftype> res(Beam->dt.data(), Beam->dt.data() + Beam->n_macroparticles);
   ASSERT_EQ(v.size(), res.size());
   ftype epsilon = 1e-6;

   for (unsigned int i = 0; i < v.size(); ++i) {
      ftype ref = v[i];
      ftype real = res[i];
      ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
            << "Testing of Beam->dt failed on i "
            << i << std::endl;
   }


   v.clear();
   res.clear();

   res = std::vector<ftype>(Beam->dE.data(), Beam->dE.data() + Beam->n_macroparticles);

   util::read_vector_from_file(v, params + "beam_dE.txt");

   ASSERT_EQ(v.size(), res.size());
   epsilon = 1e-5;

   for (unsigned int i = 0; i < v.size(); ++i) {
      ftype ref = v[i];
      ftype real = res[i];
      ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
            << "Testing of Beam->dE failed on i "
            << i << std::endl;
   }

   delete long_tracker;

}




int main(int ac, char *av[])
{
   ::testing::InitGoogleTest(&ac, av);
   return RUN_ALL_TESTS();
}
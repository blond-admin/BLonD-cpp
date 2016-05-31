#include <iostream>

#include <gtest/gtest.h>
#include <blond/math_functions.h>
#include <blond/utilities.h>
#include <blond/input_parameters/GeneralParameters.h>
#include <blond/beams/Beams.h>
#include <blond/input_parameters/RfParameters.h>
#include <blond/beams/Slices.h>
//#include "constants.h"

const ftype epsilon = 1e-8;
GeneralParameters *GP;
Beams *Beam;
RfParameters *RfP;
Slices *Slice;
int n_threads = 1;

class testGP : public ::testing::Test {

protected:
   virtual void SetUp()
   {
      ftype *momentum = new ftype[N_t + 1];
      mymath::linspace(momentum, p_i, p_f, N_t + 1);

      ftype *alpha_array = new ftype[(alpha_order + 1) * n_sections];

      std::fill_n(alpha_array, (alpha_order + 1) * n_sections, alpha);

      ftype *C_array = new ftype[n_sections];
      std::fill_n(C_array, n_sections, C);

      GP = new GeneralParameters(N_t, C_array, alpha_array, alpha_order, momentum,
                                 proton);
   }


   virtual void TearDown()
   {
      // Code here will be called immediately after each test
      // (right before the destructor).
      delete GP;
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


TEST_F(testGP, test_charge)
{
   std::string GP_params = "../unit-tests/references/GP/GP_params/";
   std::vector<ftype> v;
   util::read_vector_from_file(v, GP_params + "charge");
   //std::cout << v[0];
   ftype ref = v[0];
   ftype real = GP->charge;
   ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)));
}

TEST_F(testGP, test_mass)
{
   std::string GP_params = "../unit-tests/references/GP/GP_params/";
   std::vector<ftype> v;
   util::read_vector_from_file(v, GP_params + "mass");
   //std::cout << v[0];
   ftype ref = v[0];
   ftype real = GP->mass;
   ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)));
}

TEST_F(testGP, test_ring_radius)
{
   std::string GP_params = "../unit-tests/references/GP/GP_params/";
   std::vector<ftype> v;
   util::read_vector_from_file(v, GP_params + "ring_radius");
   //std::cout << v[0];
   ftype ref = v[0];
   ftype real = GP->ring_radius;
   ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)));
}

TEST_F(testGP, test_t_rev)
{
   std::string GP_params = "../unit-tests/references/GP/GP_params/";
   std::vector<ftype> v;
   util::read_vector_from_file(v, GP_params + "t_rev");
   //std::cout << v.size() << std::endl;
   //ASSERT_EQ(v.size(), GP->n_turns+1 );
   for (unsigned int i = 0; i < v.size(); ++i) {
      ftype ref = v[i];
      ftype real = GP->t_rev[i];
      ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)));
   }

}

TEST_F(testGP, test_cycle_time)
{
   std::string GP_params = "../unit-tests/references/GP/GP_params/";
   std::vector<ftype> v;
   util::read_vector_from_file(v, GP_params + "cycle_time");
   //ASSERT_EQ(v.size(), GP->n_turns+1);
   for (unsigned int i = 0; i < v.size(); ++i) {
      ftype ref = v[i];
      ftype real = GP->cycle_time[i];
      ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)));
   }
}

TEST_F(testGP, test_omega_rev)
{
   std::string GP_params = "../unit-tests/references/GP/GP_params/";
   std::vector<ftype> v;
   util::read_vector_from_file(v, GP_params + "omega_rev");
   //std::cout << v.size() << std::endl;
   //ASSERT_EQ(v.size(), GP->n_turns+1);
   for (unsigned int i = 0; i < v.size(); ++i) {
      ftype ref = v[i];
      ftype real = GP->omega_rev[i];
      ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)));
   }
}

TEST_F(testGP, test_eta_0)
{
   std::string GP_params = "../unit-tests/references/GP/GP_params/";
   std::vector<ftype> v;
   util::read_vector_from_file(v, GP_params + "eta_0[0]");
   //std::cout << v.size() << std::endl;
   //ASSERT_EQ(v.size(), GP->n_turns +1);
   for (unsigned int i = 0; i < v.size(); ++i) {
      ftype ref = v[i];
      ftype real = GP->eta_0[i];
      ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)));
   }
}

int main(int ac, char *av[])
{
   ::testing::InitGoogleTest(&ac, av);
   return RUN_ALL_TESTS();
}
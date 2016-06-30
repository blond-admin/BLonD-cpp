#include <iostream>
#include <string>
#include <list>

#include <unistd.h>

#include <gtest/gtest.h>
#include "math_functions.h"
#include "utilities.h"
#include "../beams/Distributions.h"
#include "../input_parameters/GeneralParameters.h"
#include "constants.h"

const ftype epsilon = 1e-8;
const std::string fixed_params = "../unit-tests/references/Bigaussian/Bigaussian_fixed_params/";

GeneralParameters *GP;
Beams *Beam;
RfParameters *RfP;
Slices *Slice;
int n_threads = 1;

class testBigaussian : public ::testing::Test {

protected:
   const long N_b = 1e9;           // Intensity
   const ftype tau_0 = 0.4e-9;          // Initial bunch length, 4 sigma [s]

   virtual void SetUp()
   {
      f_vector_2d_t momentumVec(n_sections, f_vector_t(N_t + 1));
      for (auto &v : momentumVec)
         mymath::linspace(v.data(), p_i, p_f, N_t + 1);

      f_vector_2d_t alphaVec(n_sections, f_vector_t(alpha_order+1, alpha));

      f_vector_t CVec(n_sections, C);

      f_vector_2d_t hVec(n_sections , f_vector_t(N_t + 1, h));

      f_vector_2d_t voltageVec(n_sections , f_vector_t(N_t + 1, V));

      f_vector_2d_t dphiVec(n_sections , f_vector_t(N_t + 1, dphi));


      GP = new GeneralParameters(N_t, CVec, alphaVec, alpha_order, momentumVec,
                                 proton);

      Beam = new Beams(N_p, N_b);

      RfP = new RfParameters(n_sections, hVec, voltageVec, dphiVec);


   }


   virtual void TearDown()
   {
      // Code here will be called immediately after each test
      // (right before the destructor).
      delete GP;
      delete Beam;
      delete RfP;
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



TEST_F(testBigaussian, test_sigma_dE)
{

   longitudinal_bigaussian(tau_0 / 4, 0, -1, false);

   std::vector<ftype> v;

   util::read_vector_from_file(v, fixed_params + "sigma_dE");
   ftype ref = v[0];
   ftype real = Beam->sigma_dE;
   ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)));
}

TEST_F(testBigaussian, test_sigma_dt)
{

   longitudinal_bigaussian(tau_0 / 4, 0, -1, false);

   std::vector<ftype> v;

   util::read_vector_from_file(v, fixed_params + "sigma_dt");
   ftype ref = v[0];
   ftype real = Beam->sigma_dt;
   ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)));
}

TEST_F(testBigaussian, test_dE)
{

   longitudinal_bigaussian(tau_0 / 4, 0, -1, false);

   std::vector<ftype> v;
   util::read_vector_from_file(v, fixed_params + "dE");
   for (unsigned int i = 0; i < v.size(); ++i) {
      ftype ref = v[i];
      ftype real = Beam->dE[i];
      ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)));
   }
}

TEST_F(testBigaussian, test_dt)
{

   longitudinal_bigaussian(tau_0 / 4, 0, -1, false);

   std::vector<ftype> v;
   util::read_vector_from_file(v, fixed_params + "dt");
   for (unsigned int i = 0; i < v.size(); ++i) {
      ftype ref = v[i];
      ftype real = Beam->dt[i];
      ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)));
   }
}


int main(int ac, char *av[])
{
   ::testing::InitGoogleTest(&ac, av);
   return RUN_ALL_TESTS();
}

#include <iostream>

#include <gtest/gtest.h>
#include <blond/math_functions.h>
#include <blond/utilities.h>
#include <blond/input_parameters/GeneralParameters.h>
#include <blond/beams/Beams.h>
#include <blond/input_parameters/RfParameters.h>
#include <blond/beams/Slices.h>

const ftype epsilon = 1e-8;
const std::string params = "../unit-tests/references/RFP/RFP_params/";

GeneralParameters *GP;
Beams *Beam;
RfParameters *RfP;
Slices *Slice;
int n_threads = 1;

class testRFP : public ::testing::Test {

protected:
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
   const long N_b = 1e9;           // Intensity
   const ftype tau_0 = 0.4e-9;          // Initial bunch length, 4 sigma [s]

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


TEST_F(testRFP, test_length_ratio)
{
   std::vector<ftype> v;
   util::read_vector_from_file(v, params + "length_ratio");
   //std::cout << v[0];
   ftype ref = v[0];
   ftype real = RfP->length_ratio;
   ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)));
}

TEST_F(testRFP, test_E_increment)
{
   std::vector<ftype> v;
   util::read_vector_from_file(v, params + "E_increment");
   //std::cout << v.size() << std::endl;
   for (unsigned int i = 0; i < v.size(); ++i) {
      ftype ref = v[i];
      ftype real = RfP->E_increment[i];
      ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)));
   }
}

TEST_F(testRFP, test_phi_s)
{
   std::vector<ftype> v;
   util::read_vector_from_file(v, params + "phi_s");
   //std::cout << v.size() << std::endl;
   for (unsigned int i = 0; i < v.size(); ++i) {
      ftype ref = v[i];
      ftype real = RfP->phi_s[i];
      ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)));
   }
}

TEST_F(testRFP, test_Qs)
{
   std::vector<ftype> v;
   util::read_vector_from_file(v, params + "Qs");
   //std::cout << v.size() << std::endl;
   for (unsigned int i = 0; i < v.size(); ++i) {
      ftype ref = v[i];
      ftype real = RfP->Qs[i];
      ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)));
   }
}

TEST_F(testRFP, test_omega_s0)
{
   std::vector<ftype> v;
   util::read_vector_from_file(v, params + "omega_s0");
   //std::cout << v.size() << std::endl;
   for (unsigned int i = 0; i < v.size(); ++i) {
      ftype ref = v[i];
      ftype real = RfP->omega_s0[i];
      ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)));
   }
}

TEST_F(testRFP, test_omega_RF_d)
{
   std::vector<ftype> v;
   util::read_vector_from_file(v, params + "omega_RF_d[0]");
   //std::cout << v.size() << std::endl;
   for (unsigned int i = 0; i < v.size(); ++i) {
      ftype ref = v[i];
      ftype real = RfP->omega_RF_d[i];
      ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)));
   }
}

TEST_F(testRFP, test_omega_RF)
{
   std::vector<ftype> v;
   util::read_vector_from_file(v, params + "omega_RF[0]");
   //std::cout << v.size() << std::endl;
   for (unsigned int i = 0; i < v.size(); ++i) {
      ftype ref = v[i];
      ftype real = RfP->omega_RF[i];
      ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)));
   }
}

TEST_F(testRFP, test_t_RF)
{
   std::vector<ftype> v;
   util::read_vector_from_file(v, params + "t_RF");
   //std::cout << v.size() << std::endl;
   for (unsigned int i = 0; i < v.size(); ++i) {
      ftype ref = v[i];
      ftype real = RfP->t_RF[i];
      ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)));
   }
}

int main(int ac, char *av[])
{
   ::testing::InitGoogleTest(&ac, av);
   return RUN_ALL_TESTS();
}
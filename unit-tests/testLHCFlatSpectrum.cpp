#include "globals.h"
#include "utilities.h"
#include "math_functions.h"
#include <stdio.h>
#include "../llrf/LHCFlatSpectrum.h"
#include <gtest/gtest.h>


// Simulation parameters --------------------------------------------------------

// Bunch parameters
const long int N_b = (long int) 1e9;            // Intensity
const ftype tau_0 = 0.4e-9;                     // Initial bunch length, 4 sigma [s]

// Machine and RF parameters
const ftype C = 26658.883;                      // Machine circumference [m]
const ftype p_i = 450e9;                        // Synchronous momentum [eV/c]
const long h = 35640;                           // Harmonic number
const ftype V = 6e6;                            // RF voltage [V]
const ftype dphi = 0;                           // Phase modulation/offset
const ftype gamma_t = 55.759505;                // Transition gamma
const ftype alpha = 1.0 / gamma_t / gamma_t;    // First order mom. comp. factor
const int alpha_order = 1;
const int n_sections = 1;
// Tracking details

unsigned N_t = 1000;          // Number of turns to track
unsigned N_p = 10001;         // Macro-particles

int n_threads = 1;
unsigned N_slices = 1 << 8;   // = (2^8)

GeneralParameters *GP;
Beams *Beam;
Slices *Slice;
RfParameters *RfP;

class testLHCFlatSpectrum : public ::testing::Test {

protected:

   virtual void SetUp()
   {
      f_vector_t momentum(N_t + 1);
      mymath::linspace(momentum.data(), p_i, 1.01 * p_i, N_t + 1);

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
   }


   virtual void TearDown()
   {
      // Code here will be called immediately after each test
      // (right before the destructor).
      delete GP;
      delete Beam;
      delete RfP;
   }


};



TEST_F(testLHCFlatSpectrum, constructor1)
{

   auto lhcfs = new LHCFlatSpectrum(100, 1);

   auto params = std::string("../unit-tests/references/")
                 + "LHCFlatSpectrum/constructor/test1/";
   f_vector_t v;

   util::read_vector_from_file(v, params + "fs.txt");
   //util::dump(lhcfs->fFs, "fs\n");
   //ASSERT_EQ(v.size(), res.size());

   auto epsilon = 1e-8;
   for (unsigned int i = 0; i < v.size(); ++i) {
      auto ref = v[i];
      auto real = lhcfs->fFs[i];

      ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
            << "Testing of fFs failed on i "
            << i << std::endl;
   }

   delete lhcfs;
}



int main(int ac, char *av[])
{
   ::testing::InitGoogleTest(&ac, av);
   return RUN_ALL_TESTS();
}
#include "globals.h"
#include "utilities.h"
#include "math_functions.h"
#include <../beams/Distributions.h>
#include "../llrf/PhaseLoop.h"
#include <gtest/gtest.h>
#include "../trackers/Tracker.h"

// Simulation parameters --------------------------------------------------------

// Bunch parameters
const uint N_b = 0;            // Intensity

// Machine and RF parameters
const ftype radius = 25;
const ftype C = 2 * constant::pi * radius;      // Machine circumference [m]
const ftype p_i = 310891054.809;                // Synchronous momentum [eV/c]
const uint h = 1;                               // Harmonic number
const ftype V = 8000;                           // RF voltage [V]
const ftype dphi = -constant::pi;               // Phase modulation/offset
const ftype gamma_t = 4.076750841;              // Transition gamma
const ftype alpha = 1.0 / gamma_t / gamma_t;    // First order mom. comp. factor
const uint alpha_order = 1;
const uint n_sections = 1;
// Tracking details

uint N_t = 1000;            // Number of turns to track
uint N_p = 100000;         // Macro-particles

int n_threads = 1;
uint N_slices = 200;       // = (2^8)

GeneralParameters *GP;
Beams *Beam;
Slices *Slice;
RfParameters *RfP;
// RingAndRFSection *long_tracker;

class testPLLHCF : public ::testing::Test {

protected:

   virtual void SetUp()
   {
      f_vector_2d_t momentumVec(n_sections, f_vector_t(N_t + 1, p_i));

      f_vector_2d_t alphaVec(n_sections, f_vector_t(alpha_order + 1, alpha));

      f_vector_t CVec(n_sections, C);

      f_vector_2d_t hVec(n_sections , f_vector_t(N_t + 1, h));

      f_vector_2d_t voltageVec(n_sections , f_vector_t(N_t + 1, V));

      f_vector_2d_t dphiVec(n_sections , f_vector_t(N_t + 1, dphi));

      GP = new GeneralParameters(N_t, CVec, alphaVec, alpha_order,
                                 momentumVec, proton);

      Beam = new Beams(N_p, N_b);

      RfP = new RfParameters(n_sections, hVec, voltageVec, dphiVec);

      // long_tracker = new RingAndRfSection();


      Slice = new Slices(N_slices, 0,
                         -constant::pi,
                         constant::pi,
                         cuts_unit_type::rad);

   }


   virtual void TearDown()
   {
      // Code here will be called immediately after each test
      // (right before the destructor).
      delete GP;
      delete Beam;
      delete RfP;
      delete Slice;
      // delete long_tracker;
   }

};



TEST_F(testPLLHCF, track1)
{

   longitudinal_bigaussian(10e-9, 1e6, 10, false);

   auto lhcf = new LHC_F(1.0 / 25e-6, 0, 0);


   auto params = std::string("../unit-tests/references/")
                 + "PL/LHCF/track1/";

   Slice->track();
   f_vector_t domega_RF;

   for (uint i = 0; i < N_t; ++i) {
      lhcf->track();
      domega_RF.push_back(lhcf->domega_RF);
      RfP->counter++;
   }

   f_vector_t v;

   util::read_vector_from_file(v, params + "domega_RF_mean.txt");
   // util::dump(domega_RF, "domega_RF");
   auto epsilon = 1e-2;
   auto ref = v[0];
   auto real = mymath::mean(domega_RF.data(), domega_RF.size());
   ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
         << "Testing of domega_RF_mean failed\n";

   v.clear();
   util::read_vector_from_file(v, params + "domega_RF_std.txt");
   epsilon = 1e-2;
   ref = v[0];
   real = mymath::standard_deviation(domega_RF.data(), domega_RF.size());
   ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
         << "Testing of domega_RF_std failed\n";


   delete lhcf;
}



int main(int ac, char *av[])
{
   ::testing::InitGoogleTest(&ac, av);
   return RUN_ALL_TESTS();
}

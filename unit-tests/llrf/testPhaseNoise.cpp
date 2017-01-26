#include <blond/blond.h>

#include <gtest/gtest.h>

// Simulation parameters
// --------------------------------------------------------

// Bunch parameters
const long long int N_b = 1e9; // Intensity
const double tau_0 = 0.4e-9;              // Initial bunch length, 4 sigma [s]

// Machine and RF parameters
const double C = 25558.883;                   // Machine circumference [m]
const double p_i = 450e9;                     // Synchronous momentum [eV/c]
const long long h = 35640;                   // Harmonic number
const double V = 6e6;                         // RF voltage [V]
const double dphi = 0;                        // Phase modulation/offset
const double gamma_t = 55.759505;             // Transition gamma
const double alpha = 1.0 / gamma_t / gamma_t; // First order mom. comp. factor
const int alpha_order = 1;
const int n_sections = 1;
// Tracking details

unsigned N_t = 1000;  // Number of turns to track
unsigned N_p = 10001; // Macro-particles

unsigned N_slices = 1 << 8; // = (2^8)

class testPhaseNoise : public ::testing::Test {

  protected:
    virtual void SetUp() {}

    virtual void TearDown() {
        // Code here will be called immediately after each test
        // (right before the destructor).
    }
};

/*

TEST_F(testPhaseNoise, spectrum_to_phase_noise_real1)
{

   auto freqV = mymath::arange<double>(0.0, 100.0, 2.0);
   f_vector_t spectrumV(10, 10.0);
   spectrumV.resize(50, 0);

   auto RFNoise = new PhaseNoise(freqV, spectrumV, 1, 2);
   RFNoise->spectrum_to_phase_noise();



   auto params = std::string(TEST_FILES"/")
                 + "PhaseNoise/spectrum_to_phase_noise/test1/";
   f_vector_t v;

   util::read_vector_from_file(v, params + "nt.txt");

   //ASSERT_EQ(v.size(), res.size());

   auto epsilon = 1e-8;
   for (unsigned int i = 0; i < v.size(); ++i) {
      auto ref = v[i];
      auto real = RFNoise->fNt;

      ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
            << "Testing of fNt failed on i "
            << i << std::endl;
   }

   v.clear();
   util::read_vector_from_file(v, params + "dt.txt");

   epsilon = 1e-8;
   for (unsigned int i = 0; i < v.size(); ++i) {
      auto ref = v[i];
      auto real = RFNoise->fDt;

      ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
            << "Testing of fDt failed on i "
            << i << std::endl;
   }

   v.clear();
   util::read_vector_from_file(v, params + "t.txt");
   ASSERT_EQ(v.size(), RFNoise->fT.size());

   epsilon = 1e-8;
   for (unsigned int i = 0; i < v.size(); ++i) {
      auto ref = v[i];
      auto real = RFNoise->fT[i];

      ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
            << "Testing of fT failed on i "
            << i << std::endl;
   }

   v.clear();

   util::read_vector_from_file(v, params + "dphi.txt");
   //ASSERT_EQ(v.size(), RFNoise->fDphi.size());

   auto meanV = mymath::mean(v.data(), v.size());
   auto stdV = mymath::standard_deviation(v.data(), v.size(), meanV);

   auto real = RFNoise->fDphi;

   auto meanR = mymath::mean(real.data(), real.size());
   auto stdR = mymath::standard_deviation(real.data(), real.size(), meanR);

   epsilon = 0.5;
   ASSERT_NEAR(meanV, meanR, epsilon * std::min(fabs(meanV), fabs(meanR)));

   epsilon = 0.5;
   ASSERT_NEAR(stdV, stdR, epsilon * std::min(fabs(stdV), fabs(stdR)));


   delete RFNoise;
}



TEST_F(testPhaseNoise, spectrum_to_phase_noise_complex1)
{

   auto freqV = mymath::arange<double>(0.0, 1000.0, 1.0);
   f_vector_t spectrumV(50, 0.5);
   spectrumV.resize(freqV.size(), 0);

   auto RFNoise = new PhaseNoise(freqV, spectrumV, 100, 200);
   RFNoise->spectrum_to_phase_noise(PhaseNoise::transform_t::c);



   auto params = std::string(TEST_FILES"/")
                 + "PhaseNoise/spectrum_to_phase_noise/test2/";
   f_vector_t v;

   util::read_vector_from_file(v, params + "nt.txt");

   //ASSERT_EQ(v.size(), res.size());

   auto epsilon = 1e-8;
   for (unsigned int i = 0; i < v.size(); ++i) {
      auto ref = v[i];
      auto real = RFNoise->fNt;

      ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
            << "Testing of fNt failed on i "
            << i << std::endl;
   }

   v.clear();
   util::read_vector_from_file(v, params + "dt.txt");

   epsilon = 1e-8;
   for (unsigned int i = 0; i < v.size(); ++i) {
      auto ref = v[i];
      auto real = RFNoise->fDt;

      ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
            << "Testing of fDt failed on i "
            << i << std::endl;
   }

   v.clear();
   util::read_vector_from_file(v, params + "t.txt");
   ASSERT_EQ(v.size(), RFNoise->fT.size());

   epsilon = 1e-8;
   for (unsigned int i = 0; i < v.size(); ++i) {
      auto ref = v[i];
      auto real = RFNoise->fT[i];

      ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
            << "Testing of fT failed on i "
            << i << std::endl;
   }

   v.clear();

   util::read_vector_from_file(v, params + "dphi.txt");
   //ASSERT_EQ(v.size(), RFNoise->fDphi.size());

   auto meanV = mymath::mean(v.data(), v.size());
   auto stdV = mymath::standard_deviation(v.data(), v.size(), meanV);

   auto real = RFNoise->fDphi;

   auto meanR = mymath::mean(real.data(), real.size());
   auto stdR = mymath::standard_deviation(real.data(), real.size(), meanR);

   // epsilon = 0.5;
   // ASSERT_NEAR(meanV, meanR, epsilon * std::min(fabs(meanV), fabs(meanR)));

   epsilon = 0.5;
   ASSERT_NEAR(stdV, stdR, epsilon * std::min(fabs(stdV), fabs(stdR)));


   delete RFNoise;
}

*/

int main(int ac, char* av[]) {
    ::testing::InitGoogleTest(&ac, av);
    return RUN_ALL_TESTS();
}

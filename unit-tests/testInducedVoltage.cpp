#include "globals.h"
#include "utilities.h"
#include "math_functions.h"
#include <omp.h>
#include <stdio.h>
#include "../input_parameters/GeneralParameters.h"
#include "../input_parameters/RfParameters.h"
#include "../beams/Beams.h"
#include "../beams/Slices.h"
#include "../beams/Distributions.h"
#include "../trackers/Tracker.h"
#include "../impedances/InducedVoltage.h"
#include <gtest/gtest.h>
#include <complex>


const std::string datafiles =
   "../tests/input_files/TC5_Wake_impedance/";

GeneralParameters *GP;
Beams *Beam;
Slices *Slice;
RfParameters *RfP;
//RingAndRfSection *long_tracker;
Resonators *resonator;
int n_threads = 1;



class testInducedVoltage : public ::testing::Test {
public:
   // Simulation parameters --------------------------------------------------------
// Bunch parameters
   const long int N_b = (long int) 1e10;                          // Intensity
   const ftype tau_0 = 2e-9;                       // Initial bunch length, 4 sigma [s]
// const particle_type particle = proton;
// Machine and RF parameters
   const ftype C = 6911.56;                        // Machine circumference [m]
   const ftype p_i = 25.92e9;                      // Synchronous momentum [eV/c]
//const ftype p_f = 460.005e9;                  // Synchronous momentum, final
   const long h = 4620;                            // Harmonic number
   const ftype V = 0.9e6;                          // RF voltage [V]
   const ftype dphi = 0;                           // Phase modulation/offset
   const ftype gamma_t = 1 / std::sqrt(0.00192);   // Transition gamma
   const ftype alpha = 1.0 / gamma_t / gamma_t;    // First order mom. comp. factor
   const int alpha_order = 1;
   const int n_sections = 1;
// Tracking details

   int N_t = 2;    // Number of turns to track
   int N_p = 5000000;         // Macro-particles

   int N_slices = 1 << 8; // = (2^8)


protected:

   virtual void SetUp()
   {

      omp_set_num_threads(n_threads);

      ftype *momentum = new ftype[N_t + 1];
      std::fill_n(momentum, N_t + 1, p_i);

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

      //RingAndRfSection *long_tracker = new RingAndRfSection();

      longitudinal_bigaussian(tau_0 / 4, 0, 1, false);

      Slice = new Slices(N_slices, 0, 0, 2 * constant::pi, rad);
      //util::dump(Slice->bin_centers, 10, "bin_centers\n");

      std::vector<ftype> v;
      util::read_vector_from_file(v, datafiles +
                                  "TC5_new_HQ_table.dat");
      assert(v.size() % 3 == 0);

      std::vector<ftype> R_shunt, f_res, Q_factor;

      R_shunt.reserve(v.size() / 3);
      f_res.reserve(v.size() / 3);
      Q_factor.reserve(v.size() / 3);

      for (uint i = 0; i < v.size(); i += 3) {
         f_res.push_back(v[i] * 1e9);
         Q_factor.push_back(v[i + 1]);
         R_shunt.push_back(v[i + 2] * 1e6);
      }

      resonator = new Resonators(R_shunt, f_res, Q_factor);

   }


   virtual void TearDown()
   {
      // Code here will be called immediately after each test
      // (right before the destructor).
      delete GP;
      delete Beam;
      delete RfP;
      delete Slice;
      //delete long_tracker;
   }


};

class testInducedVoltageSmall : public ::testing::Test {

public:
   // Simulation parameters --------------------------------------------------------
// Bunch parameters
   const long int N_b = (long int) 1e10;                          // Intensity
   const ftype tau_0 = 2e-9;                       // Initial bunch length, 4 sigma [s]
// const particle_type particle = proton;
// Machine and RF parameters
   const ftype C = 6911.56;                        // Machine circumference [m]
   const ftype p_i = 25.92e9;                      // Synchronous momentum [eV/c]
//const ftype p_f = 460.005e9;                  // Synchronous momentum, final
   const long h = 4620;                            // Harmonic number
   const ftype V = 0.9e6;                          // RF voltage [V]
   const ftype dphi = 0;                           // Phase modulation/offset
   const ftype gamma_t = 1 / std::sqrt(0.00192);   // Transition gamma
   const ftype alpha = 1.0 / gamma_t / gamma_t;    // First order mom. comp. factor
   const int alpha_order = 1;
   const int n_sections = 1;
// Tracking details

   int N_t = 2;    // Number of turns to track
   int N_p = 1000;         // Macro-particles

   int n_threads = 1;
   int N_slices = 1 << 5; // = (2^8)


protected:

   virtual void SetUp()
   {

      omp_set_num_threads(n_threads);

      ftype *momentum = new ftype[N_t + 1];
      std::fill_n(momentum, N_t + 1, p_i);

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

      //RingAndRfSection *long_tracker = new RingAndRfSection();

      longitudinal_bigaussian(tau_0 / 4, 0, 1, false);

      Slice = new Slices(N_slices, 0, 0, 2 * constant::pi, rad);
      //util::dump(Slice->bin_centers, 10, "bin_centers\n");

      std::vector<ftype> v;
      util::read_vector_from_file(v, datafiles +
                                  "TC5_new_HQ_table.dat");
      assert(v.size() % 3 == 0);

      std::vector<ftype> R_shunt, f_res, Q_factor;

      R_shunt.reserve(v.size() / 3);
      f_res.reserve(v.size() / 3);
      Q_factor.reserve(v.size() / 3);

      for (uint i = 0; i < v.size(); i += 3) {
         f_res.push_back(v[i] * 1e9);
         Q_factor.push_back(v[i + 1]);
         R_shunt.push_back(v[i + 2] * 1e6);
      }

      resonator = new Resonators(R_shunt, f_res, Q_factor);
   }


   virtual void TearDown()
   {
      delete GP;
      delete Beam;
      delete RfP;
      delete Slice;
   }


};


TEST_F(testInducedVoltage, InducedVoltageTime_Constructor)
{

   std::vector<Intensity *> wakeSourceList({resonator});
   //wakeSourceList.push_back(resonator);
   InducedVoltageTime *indVoltTime = new InducedVoltageTime(wakeSourceList);

   std::string params = std::string("../unit-tests/references/Impedances/")
                        + "InducedVoltage/InducedVoltageTime/";

   std::vector<ftype> v;
   util::read_vector_from_file(v, params + "time_array.txt");

   ASSERT_EQ(v.size(), indVoltTime->fTimeArray.size());

   ftype epsilon = 1e-8;

   for (unsigned int i = 0; i < v.size(); ++i) {
      ftype ref = v[i];
      ftype real = indVoltTime->fTimeArray[i];
      ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
            << "Testing of indVoltTime->fTimeArray failed on i "
            << i << std::endl;
   }
   v.clear();

   util::read_vector_from_file(v, params + "total_wake.txt");
   ASSERT_EQ(v.size(), indVoltTime->fTotalWake.size());

   epsilon = 1e-8;
   for (unsigned int i = 0; i < v.size(); ++i) {
      ftype ref = v[i];
      ftype real = indVoltTime->fTotalWake[i];
      ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
            << "Testing of indVoltTime->fTotalWake failed on i "
            << i << std::endl;
   }
   v.clear();

   util::read_vector_from_file(v, params + "cut.txt");

   for (unsigned int i = 0; i < v.size(); ++i) {
      unsigned int ref = v[i];
      unsigned int real = indVoltTime->fCut;
      ASSERT_EQ(ref, real)
            << "Testing of indVoltTime->fCut failed on i "
            << i << std::endl;
   }
   v.clear();

   util::read_vector_from_file(v, params + "fshape.txt");

   for (unsigned int i = 0; i < v.size(); ++i) {
      unsigned int ref = v[i];
      unsigned int real = indVoltTime->fShape;
      ASSERT_EQ(ref, real)
            << "Testing of fShape failed on i "
            << i << std::endl;
   }
   v.clear();


}

TEST_F(testInducedVoltage, InducedVoltageTimeReprocess)
{

   std::vector<Intensity *> wakeSourceList({resonator});
   //wakeSourceList.push_back(resonator);
   InducedVoltageTime *indVoltTime = new InducedVoltageTime(wakeSourceList);

   Slice->track();

   for (int i = 0; i < Slice->n_slices; i++)
      Slice->bin_centers[i] *= 1.1;

   indVoltTime->reprocess();

   std::string params = std::string("../unit-tests/references/Impedances/")
                        + "InducedVoltage/InducedVoltageTime/reprocess/";

   std::vector<ftype> v;
   util::read_vector_from_file(v, params + "time_array.txt");

   ASSERT_EQ(v.size(), indVoltTime->fTimeArray.size());

   ftype epsilon = 1e-8;

   for (unsigned int i = 0; i < v.size(); ++i) {
      ftype ref = v[i];
      ftype real = indVoltTime->fTimeArray[i];
      ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
            << "Testing of indVoltTime->fTimeArray failed on i "
            << i << std::endl;
   }
   v.clear();

   util::read_vector_from_file(v, params + "total_wake.txt");
   ASSERT_EQ(v.size(), indVoltTime->fTotalWake.size());

   epsilon = 1e-8;
   for (unsigned int i = 0; i < v.size(); ++i) {
      ftype ref = v[i];
      ftype real = indVoltTime->fTotalWake[i];
      ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
            << "Testing of indVoltTime->fTotalWake failed on i "
            << i << std::endl;
   }
   v.clear();

   util::read_vector_from_file(v, params + "cut.txt");

   for (unsigned int i = 0; i < v.size(); ++i) {
      unsigned int ref = v[i];
      unsigned int real = indVoltTime->fCut;
      ASSERT_EQ(ref, real)
            << "Testing of indVoltTime->fCut failed on i "
            << i << std::endl;
   }
   v.clear();

   util::read_vector_from_file(v, params + "fshape.txt");

   for (unsigned int i = 0; i < v.size(); ++i) {
      unsigned int ref = v[i];
      unsigned int real = indVoltTime->fShape;
      ASSERT_EQ(ref, real)
            << "Testing of fShape failed on i "
            << i << std::endl;
   }
   v.clear();


}



TEST_F(testInducedVoltage, induced_voltage_generation)
{
   Slice->track();

   std::vector<Intensity *> wakeSourceList({resonator});
   InducedVoltageTime *indVoltTime = new InducedVoltageTime(wakeSourceList);
   std::vector<ftype> res = indVoltTime->induced_voltage_generation();

   std::string params = std::string("../unit-tests/references/Impedances/")
                        + "InducedVoltage/InducedVoltageTime/";

   std::vector<ftype> v;
   util::read_vector_from_file(v, params + "induced_voltage.txt");

   ASSERT_EQ(v.size(), res.size());

   ftype max = *max_element(res.begin(), res.end(),
   [](ftype i, ftype j) {return fabs(i) < fabs(j);});
   ftype epsilon = 1e-9 * max;

   for (unsigned int i = 0; i < v.size(); ++i) {
      ftype ref = v[i];
      ftype real = res[i];

      ASSERT_NEAR(ref, real, epsilon)
            << "Testing of indVoltTime->fInducedVoltage failed on i "
            << i << std::endl;

   }

}


TEST_F(testInducedVoltage, induced_voltage_generation_convolution)
{
   Slice->track();

   std::vector<Intensity *> wakeSourceList({resonator});
   InducedVoltageTime *indVoltTime = new InducedVoltageTime(wakeSourceList, time_or_freq::time_domain);
   std::vector<ftype> res = indVoltTime->induced_voltage_generation();

   std::string params = std::string("../unit-tests/references/Impedances/")
                        + "InducedVoltage/InducedVoltageTime/";

   std::vector<ftype> v;
   util::read_vector_from_file(v, params + "induced_voltage_with_convolution.txt");

   ASSERT_EQ(v.size(), res.size());

   ftype max = *max_element(res.begin(), res.end(),
   [](ftype i, ftype j) {return fabs(i) < fabs(j);});
   ftype epsilon = 1e-9 * max;

   for (unsigned int i = 0; i < v.size(); ++i) {
      ftype ref = v[i];
      ftype real = res[i];

      ASSERT_NEAR(ref, real, epsilon)
            << "Testing of indVoltTime->fInducedVoltage failed on i "
            << i << std::endl;

   }

}

TEST_F(testInducedVoltage, track)
{
   Slice->track();

   std::vector<Intensity *> wakeSourceList({resonator});
   InducedVoltageTime *indVoltTime = new InducedVoltageTime(wakeSourceList);
   //std::vector<ftype> res = indVoltTime->induced_voltage_generation();
   indVoltTime->track();
   std::string params = std::string("../unit-tests/references/Impedances/")
                        + "InducedVoltage/InducedVoltageTime/";

   std::vector<ftype> v;
   util::read_vector_from_file(v, params + "beam_dE.txt");
   // WARNING checking only the fist 100 elems
   std::vector<ftype> res(Beam->dE.data(), Beam->dE.data() + 100);
   ASSERT_EQ(v.size(), res.size());

   ftype epsilon = 1e-8;
   // warning checking only the first 100 elems
   for (unsigned int i = 0; i < v.size(); ++i) {
      ftype ref = v[i];
      ftype real = res[i];
      ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
            << "Testing of Beam->dE failed on i "
            << i << std::endl;
   }

   v.clear();

}


TEST_F(testInducedVoltage, totalInducedVoltageSum)
{
   Slice->track();

   std::vector<Intensity *> wakeSourceList({resonator});
   InducedVoltageTime *indVoltTime = new InducedVoltageTime(wakeSourceList);
   //std::vector<ftype> res = indVoltTime->induced_voltage_generation();
   //indVoltTime->track();
   std::vector<InducedVoltage *> indVoltList({indVoltTime});
   TotalInducedVoltage *totVol = new TotalInducedVoltage(indVoltList);

   std::vector<ftype> res = totVol->induced_voltage_sum(200);
   auto params = std::string("../unit-tests/references/Impedances/")
                 + "InducedVoltage/TotalInducedVoltage/";

   std::vector<ftype> v;
   util::read_vector_from_file(v, params + "extIndVolt.txt");
   // WARNING checking only the fist 100 elems
   //for (auto & :)
   ASSERT_EQ(v.size(), res.size());

   ftype max = *max_element(res.begin(), res.end(),
   [](ftype i, ftype j) {return fabs(i) < fabs(j);});
   ftype epsilon = 1e-9 * max;
   // warning checking only the first 100 elems
   for (unsigned i = 0; i < v.size(); ++i) {
      ftype ref = v[i];
      ftype real = res[i];

      ASSERT_NEAR(ref, real, epsilon)
            << "Testing of extIndVolt failed on i "
            << i << std::endl;
   }

   v.clear();
   res.clear();

   res = totVol->fInducedVoltage;
   util::read_vector_from_file(v, params + "induced_voltage.txt");
   // WARNING checking only the fist 100 elems
   //for (auto & :)
   ASSERT_EQ(v.size(), res.size());

   max = *max_element(res.begin(), res.end(),
   [](ftype i, ftype j) {return fabs(i) < fabs(j);});
   epsilon = 1e-9 * max;
   // warning checking only the first 100 elems
   for (unsigned i = 0; i < v.size(); ++i) {
      ftype ref = v[i];
      ftype real = res[i];

      ASSERT_NEAR(ref, real, epsilon)
            << "Testing of fInducedVoltage failed on i "
            << i << std::endl;
   }

   v.clear();

}


TEST_F(testInducedVoltage, totalInducedVoltageTrack)
{

   //Slice->track(0, Beam->n_macroparticles);

   std::vector<Intensity *> wakeSourceList({resonator});
   InducedVoltageTime *indVoltTime = new InducedVoltageTime(wakeSourceList);
   std::vector<InducedVoltage *> indVoltList({indVoltTime});


   TotalInducedVoltage *totVol = new TotalInducedVoltage(indVoltList);

   totVol->track();
   //std::cout << "made it here\n";

   auto params = std::string("../unit-tests/references/Impedances/")
                 + "InducedVoltage/TotalInducedVoltage/";

   std::vector<ftype> v;
   util::read_vector_from_file(v, params + "beam_dE.txt");

   // WARNING checking only the fist 100 elems
   std::vector<ftype> res(Beam->dE.data(), Beam->dE.data() + 100);
   ASSERT_EQ(v.size(), res.size());

   ftype epsilon = 1e-8;
   // warning checking only the first 100 elems
   for (unsigned int i = 0; i < v.size(); ++i) {
      ftype ref = v[i];
      ftype real = res[i];

      ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
            << "Testing of Beam->dE failed on i "
            << i << std::endl;
   }

}

TEST_F(testInducedVoltage, Freq_constructor1)
{
   std::vector<Intensity *> ImpSourceList({resonator});

   auto indVoltFreq = new InducedVoltageFreq(ImpSourceList, 1e5);

   auto params = std::string("../unit-tests/references/Impedances/")
                 + "InducedVoltage/InducedVoltageFreq/constructor1/";

   std::vector<ftype> v;

   util::read_vector_from_file(v, params + "n_fft_sampling.txt");
   //ASSERT_EQ(v.size(), indVoltTime->fTotalWake.size());

   auto epsilon = 1e-8;
   for (unsigned int i = 0; i < v.size(); ++i) {
      auto ref = v[i];
      auto real = indVoltFreq->fNFFTSampling;
      ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
            << "Testing of indVoltFreq->fNFFTSampling failed on i "
            << i << std::endl;
   }
   v.clear();

   util::read_vector_from_file(v, params + "frequency_resolution.txt");

   epsilon = 1e-8;
   for (unsigned int i = 0; i < v.size(); ++i) {
      auto ref = v[i];
      auto real = indVoltFreq->fFreqResolution;
      ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
            << "Testing of indVoltFreq->fFreqResolution failed on i "
            << i << std::endl;
   }
   v.clear();

   util::read_vector_from_file(v, params + "frequency_array.txt");
   // only first 1k elements of frequency_array are tested
   epsilon = 1e-8;
   for (unsigned int i = 0; i < v.size(); ++i) {
      auto ref = v[i];
      auto real = indVoltFreq->fFreqArray[i];
      ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
            << "Testing of indVoltFreq->fFreqArray failed on i "
            << i << std::endl;
   }
   v.clear();


}

TEST_F(testInducedVoltage, Freq_constructor2)
{
   std::vector<Intensity *> ImpSourceList({resonator});

   auto indVoltFreq = new InducedVoltageFreq(ImpSourceList, 1e5,
         round_option, 100);

   auto params = std::string("../unit-tests/references/Impedances/")
                 + "InducedVoltage/InducedVoltageFreq/constructor2/";

   std::vector<ftype> v;

   util::read_vector_from_file(v, params + "n_turns_memory.txt");
   //ASSERT_EQ(v.size(), indVoltTime->fTotalWake.size());

   auto epsilon = 1e-8;
   for (unsigned int i = 0; i < v.size(); ++i) {
      auto ref = v[i];
      auto real = indVoltFreq->fNTurnsMem;
      ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
            << "Testing of indVoltFreq->fNTurnsMem failed on i "
            << i << std::endl;
   }
   v.clear();

   util::read_vector_from_file(v, params + "len_array_memory.txt");

   epsilon = 1e-8;
   for (unsigned int i = 0; i < v.size(); ++i) {
      auto ref = v[i];
      auto real = indVoltFreq->fLenArrayMem;
      ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
            << "Testing of indVoltFreq->fLenArrayMem failed on i "
            << i << std::endl;
   }
   v.clear();

   util::read_vector_from_file(v, params + "len_array_memory_extended.txt");

   epsilon = 1e-8;
   for (unsigned int i = 0; i < v.size(); ++i) {
      auto ref = v[i];
      auto real = indVoltFreq->fLenArrayMemExt;
      ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
            << "Testing of indVoltFreq->fLenArrayMemExt failed on i "
            << i << std::endl;
   }
   v.clear();


   util::read_vector_from_file(v, params + "n_points_fft.txt");

   epsilon = 1e-8;
   for (unsigned int i = 0; i < v.size(); ++i) {
      auto ref = v[i];
      auto real = indVoltFreq->fNPointsFFT;
      ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
            << "Testing of indVoltFreq->fNPointsFFT failed on i "
            << i << std::endl;
   }
   v.clear();

   util::read_vector_from_file(v, params + "frequency_array_memory.txt");
   // Only fist 1k elements of frequency_array_memory are tested
   epsilon = 1e-8;
   for (unsigned int i = 0; i < v.size(); ++i) {
      auto ref = v[i];
      auto real = indVoltFreq->fFreqArrayMem[i];
      ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
            << "Testing of indVoltFreq->fFreqArrayMem failed on i "
            << i << std::endl;
   }
   v.clear();

   util::read_vector_from_file(v, params + "time_array_memory.txt");
   // Only fist 100 elements of frequency_array_memory are tested
   epsilon = 1e-8;
   for (unsigned int i = 0; i < v.size(); ++i) {
      auto ref = v[i];
      auto real = indVoltFreq->fTimeArrayMem[i];
      ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
            << "Testing of indVoltFreq->fTimeArrayMem failed on i "
            << i << std::endl;
   }
   v.clear();

   util::read_vector_from_file(v, params + "total_impedance_memory.txt");
   // Only fist 1000 elements of total_impedance_memory are tested
   epsilon = 1e-8;
   for (unsigned int i = 0; i < v.size(); ++i) {
      auto ref = v[i];
      auto real = std::abs(indVoltFreq->fTotalImpedanceMem[i]);
      ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
            << "Testing of indVoltFreq->fTotalImpedanceMem failed on i "
            << i << std::endl;
   }
   v.clear();


}

TEST_F(testInducedVoltage, Freq_sum_impedances1)
{
   std::vector<Intensity *> ImpSourceList({resonator});

   auto indVoltFreq = new InducedVoltageFreq(ImpSourceList, 1e5);

   auto freq_array = mymath::rfftfreq(Slice->n_slices);

   indVoltFreq->sum_impedances(freq_array);

   auto params = std::string("../unit-tests/references/Impedances/")
                 + "InducedVoltage/InducedVoltageFreq/sum_impedances/";

   std::vector<ftype> v;

   util::read_vector_from_file(v, params + "total_impedance.txt");

   ASSERT_EQ(v.size(), indVoltFreq->fTotalImpedance.size());

   auto epsilon = 1e-8;
   for (unsigned int i = 0; i < v.size(); ++i) {
      auto ref = v[i];
      auto real = std::abs(indVoltFreq->fTotalImpedance[i]);
      ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
            << "Testing of indVoltFreq->fTotalImpedance failed on i "
            << i << std::endl;
   }
}


TEST_F(testInducedVoltage, Freq_sum_impedances2)
{
   std::vector<Intensity *> ImpSourceList({resonator});

   auto indVoltFreq = new InducedVoltageFreq(ImpSourceList, 1e5);

   auto freq_array = mymath::rfftfreq(Slice->n_slices,
                                      Slice->bin_centers[1] - Slice->bin_centers[0]);

   indVoltFreq->sum_impedances(freq_array);

   auto params = std::string("../unit-tests/references/Impedances/")
                 + "InducedVoltage/InducedVoltageFreq/sum_impedances2/";

   std::vector<ftype> v;

   util::read_vector_from_file(v, params + "total_impedance.txt");
   ASSERT_EQ(v.size(), indVoltFreq->fTotalImpedance.size());

   auto epsilon = 1e-8;
   for (unsigned int i = 0; i < v.size(); ++i) {
      auto ref = v[i];
      auto real = std::abs(indVoltFreq->fTotalImpedance[i]);
      ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
            << "Testing of indVoltFreq->fTotalImpedance failed on i "
            << i << std::endl;
   }
}



TEST_F(testInducedVoltage, Freq_reprocess1)
{
   std::vector<Intensity *> ImpSourceList({resonator});

   auto indVoltFreq = new InducedVoltageFreq(ImpSourceList, 1e5);

   for (int i = 0; i < Slice->n_slices; ++i) {
      Slice->bin_centers[i] = 1.1 * Slice->bin_centers[i];
   }

   indVoltFreq->reprocess();

   auto params = std::string("../unit-tests/references/Impedances/")
                 + "InducedVoltage/InducedVoltageFreq/reprocess1/";

   std::vector<ftype> v;

   util::read_vector_from_file(v, params + "n_fft_sampling.txt");
   auto epsilon = 1e-8;
   for (unsigned int i = 0; i < v.size(); ++i) {
      auto ref = v[i];
      auto real = indVoltFreq->fNFFTSampling;
      ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
            << "Testing of indVoltFreq->fNFFTSampling failed on i "
            << i << std::endl;
   }
   v.clear();

   util::read_vector_from_file(v, params + "frequency_resolution.txt");

   epsilon = 1e-8;
   for (unsigned int i = 0; i < v.size(); ++i) {
      auto ref = v[i];
      auto real = indVoltFreq->fFreqResolution;
      ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
            << "Testing of indVoltFreq->fFreqResolution failed on i "
            << i << std::endl;
   }
   v.clear();

   util::read_vector_from_file(v, params + "frequency_array.txt");
   // only first 1k elements of frequency_array are tested
   epsilon = 1e-8;
   for (unsigned int i = 0; i < v.size(); ++i) {
      auto ref = v[i];
      auto real = indVoltFreq->fFreqArray[i];
      ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
            << "Testing of indVoltFreq->fFreqArray failed on i "
            << i << std::endl;
   }
   v.clear();
}


TEST_F(testInducedVoltageSmall, Freq_induced_voltage_generation1)
{
   std::vector<Intensity *> ImpSourceList({resonator});

   auto indVoltFreq = new InducedVoltageFreq(ImpSourceList, 1e5);
   Slice->track();

   indVoltFreq->induced_voltage_generation();
   //util::dump(indVoltFreq->fInducedVoltage.data(), 10, "irfft\n");
   auto params = std::string("../unit-tests/references/Impedances/")
                 + "InducedVoltage/InducedVoltageFreq/induced_voltage_generation1/";

   std::vector<ftype> v;

   util::read_vector_from_file(v, params + "induced_voltage.txt");

   ASSERT_EQ(v.size(), indVoltFreq->fInducedVoltage.size());

   auto epsilon = 1e-4;
   for (unsigned int i = 0; i < v.size(); ++i) {
      auto ref = v[i];
      ftype real = indVoltFreq->fInducedVoltage[i];
      ASSERT_NEAR(ref, real, epsilon * std::max(fabs(ref), fabs(real)))
            << "Testing of indVoltFreq->fInducedVoltage failed on i "
            << i << std::endl;
   }
   v.clear();


}



int main(int ac, char *av[])
{
   ::testing::InitGoogleTest(&ac, av);
   return RUN_ALL_TESTS();
}
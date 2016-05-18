
#include "InducedVoltage.h"

#include "utilities.h"
#include "constants.h"
#include "math_functions.h"
#include "Ham.h"
#include "globals.h"


InducedVoltageTime::InducedVoltageTime(std::vector<Intensity *> &WakeSourceList,
                                       time_or_freq TimeOrFreq)
{
   // Induced voltage derived from the sum of
   // several wake fields (time domain).*

   // *Wake sources inputed as a list (eg: list of BBResonators objects)*
   fWakeSourceList = WakeSourceList;

   // *Time array of the wake in [s]*
   fTimeArray = std::vector<ftype>();

   // *Total wake array of all sources in* [:math:`\Omega / s`]
   fTotalWake = std::vector<ftype>();

   // *Induced voltage from the sum of the wake sources in [V]*
   fInducedVoltage = std::vector<ftype>();

   // Pre-processing the wakes
   fTimeArray.resize(Slice->n_slices);
   for (unsigned int i = 0; i < fTimeArray.size(); ++i) {
      fTimeArray[i] = Slice->bin_centers[i] - Slice->bin_centers[0];
   }
   sum_wakes(fTimeArray);

   fCut = fTimeArray.size() + Slice->n_slices - 1;
   fShape = next_regular(fCut);

   fTimeOrFreq = TimeOrFreq;

}


void InducedVoltageTime::track() {}

void InducedVoltageTime::sum_wakes(std::vector<ftype> &TimeArray)
{
   // *Summing all the wake contributions in one total wake.*
   fTotalWake.resize(TimeArray.size());
   std::fill(fTotalWake.begin(), fTotalWake.end(), 0);
   for (Intensity *i : fWakeSourceList) {

      i->wake_calc(TimeArray);
      std::transform(fTotalWake.begin(), fTotalWake.end(),
                     i->fWake.begin(), fTotalWake.begin(),
                     std::plus<ftype>());

   }

}

void InducedVoltageTime::reprocess()
{
   // *Reprocess the wake contributions with respect to the new_slicing.*
   // WARNING As Slice is a global variable,
   // users will have to change this variable and call reprocess()
   fTimeArray.resize(Slice->n_slices);
   for (unsigned int i = 0; i < fTimeArray.size(); ++i) {
      fTimeArray[i] = Slice->bin_centers[i] - Slice->bin_centers[0];
   }
   sum_wakes(fTimeArray);

   fCut = fTimeArray.size() + Slice->n_slices - 1;
   fShape = next_regular(fCut);


}

// TODO resolve the use of length
void InducedVoltageTime::induced_voltage_generation(unsigned int length)
{

   // Method to calculate the induced voltage from wakes with convolution.*

   // TODO test this function

   std::vector<ftype> inducedVoltage;

   //util::dump(inducedVoltage.data(), inducedVoltage.size(), "inducedVoltage\n");
   if (fTimeOrFreq == freq) {
      std::vector<complex_t> fft1, fft2;
      std::vector<ftype> in(Slice->n_macroparticles,
                            Slice->n_macroparticles + Slice->n_slices);
      
      mymath::real_to_complex(in, fft1);

      in.clear();
      in = fTotalWake;
      //std::copy(fTotalWake.begin(), fTotalWake.end(), in.begin());
      mymath::real_to_complex(in, fft2);

      mymath::fft(fft1, fShape, fft1);
      mymath::fft(fft2, fShape, fft2);

      //std::vector<ftype> temp;

      //mymath::complex_to_real(fft1, temp);
      //util::dump(temp.data(), 10, "fft1\n");
      //temp.clear();
      //mymath::complex_to_real(fft1, temp);

      //util::dump(temp.data(), 10, "fft2\n");
      //std::vector<complex_t> fft3;
      //fft3.resize(fShape);
      //for (int i = 0; i < 10; ++i) {
      //   std::cout << "fft1: " << fft1[i] << "\n";
      //   std::cout << "fft2: " << fft2[i] << "\n";

      //   std::cout << "fft3: " << fft1[i] * fft2[i] << "\n";
         //fft3[i] = fft1[i] * fft2[i];
         //printf("fft3: %lf\n", std::abs(fft1[i] * fft2[i]));
      //}

      std::transform(fft1.begin(), fft1.end(), fft2.begin(),
                     fft1.begin(), std::multiplies<complex_t>());

      //temp.clear();
      //mymath::complex_to_real(fft3, temp);
      //util::dump(temp.data(), 10, "fft1*fft2\n");
      //for (int i = 0; i < 10; ++i)
      //{
      //   printf("fft1*fft2 = %lf\n",std::abs(fft1[i]));
      //}
      mymath::ifft(fft1, fShape, fft1);

      mymath::complex_to_real(fft1, inducedVoltage);

      const ftype factor = - GP->charge * constant::e *
                           Beam->intensity / Beam->n_macroparticles;
      //printf("GP->charge = %e\n", GP->charge);
      //printf("e = %e\n", constant::e);
      //printf("Beam->intensity = %ld\n", Beam->intensity);
      //printf("Beam->n_macroparticles = %d\n", Beam->n_macroparticles);
      std::transform(inducedVoltage.begin(), inducedVoltage.end(),
                     inducedVoltage.begin(),
                     std::bind1st(std::multiplies<ftype>(), factor));

   } else if (fTimeOrFreq == time) {


   } else {
      dprintf("Error: Only freq or time are allowed\n");
      exit(-1);
   }

   //fInducedVoltage = inducedVoltage;
   //fInducedVoltage(inducedVoltage.begin(),
   //                inducedVoltage.begin() + Slice->n_slices);
   fInducedVoltage = inducedVoltage;
   fInducedVoltage.resize(Slice->n_slices, 0);
}





InducedVoltageFreq::InducedVoltageFreq() {}

void InducedVoltageFreq::track() {}

void InducedVoltageFreq::sum_impedances() {}

void InducedVoltageFreq::reprocess() {}

void InducedVoltageFreq::induced_voltage_generation(unsigned int length) {}

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
   fInducedVoltage = std::vector<complex_t>();

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
   std::vector<complex_t> inducedVoltage;

   //util::dump(inducedVoltage.data(), inducedVoltage.size(), "inducedVoltage\n");
   if (fTimeOrFreq == freq) {
      std::vector<complex_t> fft1, fft2;
      std::vector<ftype> in(Slice->n_macroparticles,
                            Slice->n_macroparticles + Slice->n_slices);
      
      mymath::rfft(in, fShape, fft1);
      mymath::rfft(fTotalWake, fShape, fft2);

      std::transform(in.begin(), in.end(), fTotalWake.begin(),
                    in.begin(), std::multiplies<ftype>());
      
      assert(in.size() == fTotalWake.size());

      mymath::irfft(in, fShape, inducedVoltage);
   } else if (fTimeOrFreq == time) {


   } else {
      dprintf("Error: Only freq or time are allowed\n");
      exit(-1);
   }

   fInducedVoltage = inducedVoltage;
   //fInducedVoltage(inducedVoltage.begin(),
   //                inducedVoltage.begin() + Slice->n_slices);

}





InducedVoltageFreq::InducedVoltageFreq() {}

void InducedVoltageFreq::track() {}

void InducedVoltageFreq::sum_impedances() {}

void InducedVoltageFreq::reprocess() {}

void InducedVoltageFreq::induced_voltage_generation(unsigned int length) {}
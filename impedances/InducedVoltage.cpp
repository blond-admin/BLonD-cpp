
#include "InducedVoltage.h"

#include "utilities.h"
#include "constants.h"
#include "math_functions.h"
#include "Ham.h"
#include "globals.h"


inline void InducedVoltage::linear_interp_kick(
   const ftype *__restrict__ beam_dt,
   ftype *__restrict__ beam_dE,
   const ftype *__restrict__ voltage_array,
   const ftype *__restrict__ bin_centers,
   const int n_slices,
   const int n_macroparticles,
   const ftype acc_kick)
{

   // double a;
   // int i;
   // double fbin;
   // int ffbin;
   // double voltageKick;
   double inv_bin_width = (n_slices - 1) / (bin_centers[n_slices - 1]
                          - bin_centers[0]);

   for (int i = 0; i < n_macroparticles; i++) {
      ftype a = beam_dt[i];
      ftype fbin = (a - bin_centers[0]) * inv_bin_width;
      int ffbin = (int)(fbin);
      ftype voltageKick;
      if ((a < bin_centers[0]) || (a > bin_centers[n_slices - 1]))
         voltageKick = 0.;
      else
         voltageKick = voltage_array[ffbin] + (a - bin_centers[ffbin])
                       * (voltage_array[ffbin + 1] - voltage_array[ffbin])
                       * inv_bin_width;
      beam_dE[i] = beam_dE[i] + voltageKick + acc_kick;
   }


}

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


void InducedVoltageTime::track()
{
   // Tracking Method
   std::vector<ftype> v = induced_voltage_generation();

   std::transform(v.begin(), v.end(), v.begin(),
                  std::bind1st(std::multiplies<ftype>(), GP->charge));

   linear_interp_kick(Beam->dt, Beam->dE, v.data(), Slice->bin_centers,
                      Slice->n_slices, Beam->n_macroparticles, 0.0);
}

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

std::vector<ftype> InducedVoltageTime::induced_voltage_generation(unsigned int length)
{

   // Method to calculate the induced voltage from wakes with convolution.*

   std::vector<ftype> inducedVoltage;

   const ftype factor = - GP->charge * constant::e *
                        Beam->intensity / Beam->n_macroparticles;

   if (fTimeOrFreq == freq_domain) {
      std::vector<complex_t> fft1, fft2;
      std::vector<ftype> in(Slice->n_macroparticles,
                            Slice->n_macroparticles + Slice->n_slices);

      mymath::real_to_complex(in, fft1);

      in.clear();
      in = fTotalWake;
      mymath::real_to_complex(in, fft2);

      mymath::fft(fft1, fShape, fft1);
      mymath::fft(fft2, fShape, fft2);

      std::transform(fft1.begin(), fft1.end(), fft2.begin(),
                     fft1.begin(), std::multiplies<complex_t>());

      mymath::ifft(fft1, fShape, fft1);

      mymath::complex_to_real(fft1, inducedVoltage);

      std::transform(inducedVoltage.begin(),
                     inducedVoltage.end(),
                     inducedVoltage.begin(),
                     std::bind1st(std::multiplies<ftype>(), factor));

   } else if (fTimeOrFreq == time_domain) {
      std::vector<ftype> temp(Slice->n_macroparticles,
                              Slice->n_macroparticles + Slice->n_slices);
      inducedVoltage =
         mymath::convolution(fTotalWake, temp);

      std::transform(inducedVoltage.begin(),
                     inducedVoltage.end(),
                     inducedVoltage.begin(),
                     std::bind1st(std::multiplies<ftype>(), factor));

   } else {
      dprintf("Error: Only freq_domain or time_domain are allowed\n");
      exit(-1);
   }



   inducedVoltage.resize(
      std::max((unsigned int)Slice->n_slices, length), 0);


   // TODO do we really need both of these?
   fInducedVoltage = inducedVoltage;

   return inducedVoltage;

}





InducedVoltageFreq::InducedVoltageFreq() {}

void InducedVoltageFreq::track() {}

void InducedVoltageFreq::sum_impedances() {}

void InducedVoltageFreq::reprocess() {}

std::vector<ftype> InducedVoltageFreq::induced_voltage_generation(unsigned int length)
{
   return std::vector<ftype>();
}
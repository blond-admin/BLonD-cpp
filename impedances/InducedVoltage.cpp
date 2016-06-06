
#include "InducedVoltage.h"

#include "utilities.h"
#include "constants.h"
#include "math_functions.h"
#include "Ham.h"
#include "globals.h"


std::vector<ftype> TotalInducedVoltage::induced_voltage_generation(uint length)
{
   return std::vector<ftype>();
}


inline void InducedVoltage::linear_interp_kick(
   const ftype *__restrict__ beam_dt,
   ftype *__restrict__ beam_dE,
   const ftype *__restrict__ voltage_array,
   const ftype *__restrict__ bin_centers,
   const int n_slices,
   const int n_macroparticles,
   const ftype acc_kick)
{

   //int LOOP_UNROLL = 8;
   //LOOP_UNROLL = atoi(util::GETENV("LOOP_UNROLL")) ? atoi(util::GETENV("LOOP_UNROLL")) : LOOP_UNROLL;


   const ftype binFirst = bin_centers[0];
   const ftype binLast = bin_centers[n_slices - 1];

   const ftype inv_bin_width = (n_slices - 1) / (binLast - binFirst);
   #pragma omp parallel for
   for (int i = 0; i < n_macroparticles; i++) {
      const ftype a = beam_dt[i];
      const int ffbin = static_cast<int>((a - binFirst) * inv_bin_width);
      const ftype voltageKick = ((a < binFirst) || (a > binLast)) ?
                                0 : voltage_array[ffbin] + (a - bin_centers[ffbin])
                                * (voltage_array[ffbin + 1] - voltage_array[ffbin])
                                * inv_bin_width;
      beam_dE[i] += voltageKick + acc_kick;
   }

   //int * __restrict b;
   //ftype *__restrict a = (ftype *) malloc(LOOP_UNROLL * sizeof(ftype));//[LOOP_UNROLL];
   //int *__restrict ffbin = (int *) malloc(LOOP_UNROLL * sizeof(int));;//[LOOP_UNROLL];
   //ftype *__restrict voltageKick = (ftype *) malloc(LOOP_UNROLL * sizeof(ftype));;//[LOOP_UNROLL];

   /*
   #pragma omp parallel for //private(a, ffbin, voltageKick)
   for (int i = 0; i < n_macroparticles; i++) {
      std::vector<ftype> a(LOOP_UNROLL);
      std::vector<int> ffbin(LOOP_UNROLL);
      std::vector<ftype> voltageKick(LOOP_UNROLL);

      for (int j = 0; j < LOOP_UNROLL && i + j < n_macroparticles; ++j) {
         a[j] = beam_dt[i + j];
         //ftype fbin = (a - binFirst) * inv_bin_width;
         ffbin[j] = static_cast<int>((a[j] - binFirst) * inv_bin_width);
         //unsigned ffbin = (unsigned)(fbin);
         voltageKick[j] = ((a[j] < binFirst) || (a[j] > binLast)) ?
                          0 : voltage_array[ffbin[j]] + (a[j] - bin_centers[ffbin[j]])
                          * (voltage_array[ffbin[j] + 1] - voltage_array[ffbin[j]])
                          * inv_bin_width;
         beam_dE[i + j] += voltageKick[j] + acc_kick;
      }
   }
   */

   //ftype inv_bin_width = (n_slices-1) / (bin_centers[n_slices-1] - bin_centers[0]);
   /*
   double inv_bin_width = (n_slices-1) / (bin_centers[n_slices-1] - bin_centers[0]);
   #pragma omp parallel for
   for (int i = 0; i < n_macroparticles; i++) {
      double a = beam_dt[i];
      double fbin = (a - bin_centers[0]) * inv_bin_width;
      int ffbin = (int)(fbin);
      double voltageKick;
      if ((a < bin_centers[0])||(a > bin_centers[n_slices-1]))
         voltageKick = 0.;
      else
         voltageKick = voltage_array[ffbin] + (a - bin_centers[ffbin]) * (voltage_array[ffbin+1]-voltage_array[ffbin]) * inv_bin_width;
      beam_dE[i] = beam_dE[i] + voltageKick + acc_kick;
    }
   */

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
   for (uint i = 0; i < fTimeArray.size(); ++i) {
      fTimeArray[i] = Slice->bin_centers[i] - Slice->bin_centers[0];
   }
   sum_wakes(fTimeArray);

   fCut = fTimeArray.size() + Slice->n_slices - 1;
   fShape = next_regular(fCut);

   fTimeOrFreq = TimeOrFreq;

}


inline void InducedVoltageTime::track()
{
   // Tracking Method
   std::vector<ftype> v = this->induced_voltage_generation();

   //std::cout << "induced v size is " << v.size() << "\n";

   std::transform(v.begin(), v.end(), v.begin(),
                  std::bind1st(std::multiplies<ftype>(),
                               GP->charge));

   linear_interp_kick(Beam->dt.data(), Beam->dE.data(), v.data(),
                      Slice->bin_centers, Slice->n_slices,
                      Beam->n_macroparticles, 0.0);

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
   for (uint i = 0; i < fTimeArray.size(); ++i) {
      fTimeArray[i] = Slice->bin_centers[i] - Slice->bin_centers[0];
   }
   sum_wakes(fTimeArray);

   fCut = fTimeArray.size() + Slice->n_slices - 1;
   fShape = next_regular(fCut);


}

std::vector<ftype> InducedVoltageTime::induced_voltage_generation(uint length)
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

   inducedVoltage.resize((uint) Slice->n_slices);
   fInducedVoltage = inducedVoltage;

   if (length > 0)
      inducedVoltage.resize(
         std::min((uint) Slice->n_slices, length), 0);

   //std::cout << "inducedVoltage size is " << inducedVoltage.size() << "\n";
   return inducedVoltage;

}


InducedVoltageFreq::InducedVoltageFreq(
   std::vector<Intensity *> &impedanceSourceList,
   ftype freqResolutionInput,
   freq_res_option_t freq_res_option,
   uint NTurnsMem,
   bool recalculationImpedance,
   bool saveIndividualVoltages)
{
   fNTurnsMem = NTurnsMem;

   fImpedanceSourceList = impedanceSourceList;
   fFreqResolutionInput = freqResolutionInput;

   // *Length of one slice.*
   auto timeResolution = (Slice->bin_centers[1] - Slice->bin_centers[0]);
   fRecalculationImpedance = recalculationImpedance;
   fFreqResOption = freq_res_option;

   if (fNTurnsMem == 0) {

      if (fFreqResolutionInput == 0) {
         fNFFTSampling = Slice->n_slices;
      } else {
         int a;
         ftype b = 1 / (fFreqResolutionInput * timeResolution);
         switch (fFreqResOption) {
            case freq_res_option_t::round_option:
               a = std::round(b);
               break;
            case freq_res_option_t::ceil_option:
               a = std::ceil(b);
               break;
            case freq_res_option_t::floor_option:
               a = std::floor(b);
               break;
            default:
               std::cerr << "The input freq_res_option is not recognized\n";
               exit(-1);
               break;
         }
         fNFFTSampling = next_regular(a);


         if (fNFFTSampling < (uint) Slice->n_slices) {
            std::cerr << "The input frequency resolution step is too big, and the whole\n"
                      << "bunch is not sliced... The number of sampling points for the\n"
                      << "FFT is corrected in order to sample the whole bunch (and\n"
                      << "you might consider changing the input in order to have\n"
                      << "a finer resolution\n";
            fNFFTSampling = next_regular(Slice->n_slices);
         }
      }

      fFreqResolution = 1 / (fNFFTSampling * timeResolution);

      //self.frequency_array = rfftfreq(self.n_fft_sampling, self.slices.bin_centers[1] - self.slices.bin_centers[0])
      fFreqArray = mymath::rfftfreq(fNFFTSampling, timeResolution);
      sum_impedances(fFreqArray);

      fSaveIndividualVoltages = saveIndividualVoltages;
      if (fSaveIndividualVoltages) {
         // Do I really need to store the length??
         uint length = fImpedanceSourceList.size();
         fMatrixSaveIndividualImpedances = complex_vector_t(length * fFreqArray.size(), 0);
         fMatrixSaveIndividualVoltages = f_vector_t(length * Slice->n_slices, 0);
         for (uint i = 0; i < length; ++i) {
            const uint row_width = fImpedanceSourceList[i]->fImpedance.size();
            for (uint j = 0; j < row_width; ++j) {
               fMatrixSaveIndividualImpedances[i * row_width + j] =
                  fImpedanceSourceList[i]->fImpedance[j];
            }
         }
      }

   } else {
      fNTurnsMem = NTurnsMem;
      fLenArrayMem = (fNTurnsMem + 1) * Slice->n_slices;
      fLenArrayMemExt = (fNTurnsMem + 2) * Slice->n_slices;
      fNPointsFFT = next_regular(fLenArrayMemExt);
      fFreqArrayMem = mymath::rfftfreq(fNPointsFFT, timeResolution);
      fTotalImpedanceMem = complex_vector_t(fFreqArrayMem.size(), complex_t(0, 0));

      fTimeArrayMem.reserve(fNTurnsMem * Slice->n_slices);
      const ftype factor = Slice->edges[Slice->n_slices + 1] - Slice->edges[0];
      for (uint i = 0; i < fNTurnsMem + 1; ++i) {
         for (uint j = 0; j < (uint) Slice->n_slices; ++j) {
            fTimeArrayMem.push_back(Slice->bin_centers[j] + factor * i);
         }
      }

      for (const auto &impObj : fImpedanceSourceList) {
         impObj->imped_calc(fFreqArrayMem);
         std::transform(impObj->fImpedance.begin(),
                        impObj->fImpedance.end(),
                        fTotalImpedanceMem.begin(),
                        fTotalImpedanceMem.begin(),
                        std::plus<complex_t>());
      }


   }

}

// TODO test this function
void InducedVoltageFreq::track()
{
   // Tracking Method
   induced_voltage_generation();
   auto v = fInducedVoltage;
   std::transform(v.begin(), v.end(), v.begin(),
                  std::bind1st(std::multiplies<ftype>(),
                               GP->charge));

   linear_interp_kick(Beam->dt.data(), Beam->dE.data(), v.data(),
                      Slice->bin_centers, Slice->n_slices,
                      Beam->n_macroparticles, 0.0);

}

// TODO test this function
void InducedVoltageFreq::sum_impedances(f_vector_t &freq_array)
{
   fTotalImpedance.resize(freq_array.size());
   std::fill(fTotalImpedance.begin(),
             fTotalImpedance.end(),
             complex_t(0, 0));
   for (auto &i : fImpedanceSourceList) {
      i->imped_calc(freq_array);
      std::transform(fTotalImpedance.begin(), fTotalImpedance.end(),
                     i->fImpedance.begin(), fTotalImpedance.begin(),
                     std::plus<complex_t>());

   }
}


// TODO test this function
void InducedVoltageFreq::reprocess()
{
   auto timeResolution = (Slice->bin_centers[1] - Slice->bin_centers[0]);
   if (fFreqResolutionInput == 0) {
      fNFFTSampling = Slice->n_slices;
   } else {
      int a;
      ftype b = 1 / (fFreqResolutionInput * timeResolution);
      switch (fFreqResOption) {
         case freq_res_option_t::round_option:
            a = std::round(b);
            break;
         case freq_res_option_t::ceil_option:
            a = std::ceil(b);
            break;
         case freq_res_option_t::floor_option:
            a = std::floor(b);
            break;
         default:
            std::cerr << "The input freq_res_option is not recognized\n";
            exit(-1);
            break;
      }
      fNFFTSampling = next_regular(a);


      if (fNFFTSampling < (uint) Slice->n_slices) {
         std::cerr << "The input frequency resolution step is too big, and the whole\n"
                   << "bunch is not sliced... The number of sampling points for the\n"
                   << "FFT is corrected in order to sample the whole bunch (and\n"
                   << "you might consider changing the input in order to have\n"
                   << "a finer resolution\n";
         fNFFTSampling = next_regular(Slice->n_slices);
      }
   }

   fFreqResolution = 1 / (fNFFTSampling * timeResolution);

   Slice->beam_spectrum_generation(fNFFTSampling, true);
   fFreqArray = Slice->fBeamSpectrumFreq;

   fTotalImpedance.clear();
   sum_impedances(fFreqArray);

}

std::vector<ftype> InducedVoltageFreq::induced_voltage_generation(uint length)
{
   if (fRecalculationImpedance)
      sum_impedances(fFreqArray);

   Slice->beam_spectrum_generation(fNFFTSampling);

   const auto n = fImpedanceSourceList.size();
   const auto factor = - GP->charge * constant::e *
                       Beam->ratio * Slice->fBeamSpectrumFreq[1]
                       * 2 * (Slice->fBeamSpectrum.size() - 1);
   if (fSaveIndividualVoltages) {

      for (uint i = 0; i < n; ++i) {
         f_vector_t res;
         complex_vector_t in(Slice->fBeamSpectrum.size());

         for (uint j = 0; j < in.size(); ++j) {
            in[j] = fMatrixSaveIndividualImpedances[j * n + i]
                    * Slice->fBeamSpectrum[j];
         }

         mymath::irfft(in, res);

         assert((int) res.size() >= Slice->n_slices);

         res.resize(Slice->n_slices);

         std::transform(res.begin(),
                        res.end(),
                        res.begin(),
                        std::bind1st(
                           std::multiplies<ftype>(),
                           factor));

         for (uint j = 0; j < (uint) Slice->n_slices; ++j) {
            fMatrixSaveIndividualVoltages[j * n + i] = res[j];
         }
      }

      fInducedVoltage.clear();
      fInducedVoltage.resize(Slice->fBeamSpectrum.size());
      for (uint i = 0; i < Slice->fBeamSpectrum.size(); ++i) {
         ftype sum = 0.0;
         for (uint j = 0; j < n; ++j) {
            sum += fMatrixSaveIndividualVoltages[i * n + j];
         }
         fInducedVoltage[i] = sum;
      }

      return std::vector<ftype>();

   } else {
      f_vector_t res;
      complex_vector_t in(Slice->fBeamSpectrum.size());
      for (uint j = 0; j < in.size(); ++j) {
         in[j] = fTotalImpedance[j] * Slice->fBeamSpectrum[j];
      }
      mymath::irfft(in, res);

      assert((int) res.size() >= Slice->n_slices);

      res.resize((uint)Slice->n_slices);

      fInducedVoltage = res;

      if (length > res.size())
         res.resize(length, 0);
      else
         res.resize(length);

      return res;
   }

}


TotalInducedVoltage::TotalInducedVoltage(
   std::vector<InducedVoltage *> &InducedVoltageList,
   uint NTurnsMemory,
   std::vector<ftype> RevTimeArray)
{
   fInducedVoltageList = InducedVoltageList;
   fNTurnsMemory = NTurnsMemory;
   fInducedVoltage = std::vector<ftype>();
   fTimeArray = std::vector<ftype> (Slice->bin_centers,
                                    Slice->bin_centers + Slice->n_slices);

}

void TotalInducedVoltage::track()
{
   //std::cout << "I am here\n";
   this->induced_voltage_sum();
   auto v = this->fInducedVoltage;
   //std::cout << "total v size is " << v.size() << "\n";

   std::transform(v.begin(), v.end(), v.begin(),
                  std::bind1st(std::multiplies<ftype>(),
                               GP->charge));
   //std::cout << "I am here\n";

   linear_interp_kick(Beam->dt.data(), Beam->dE.data(), v.data(),
                      Slice->bin_centers, Slice->n_slices,
                      Beam->n_macroparticles, 0.0);
   //std::cout << "I am here\n";

}

void TotalInducedVoltage::track_memory() {}

void TotalInducedVoltage::track_ghosts_particles() {}

void TotalInducedVoltage::reprocess()
{
   for (auto &v : fInducedVoltageList)
      v->reprocess();
}

std::vector<ftype> TotalInducedVoltage::induced_voltage_sum(uint length)
{
   // Method to sum all the induced voltages in one single array.
   std::vector<ftype> tempIndVolt;
   std::vector<ftype> extIndVolt;

   for (auto &v : fInducedVoltageList) {
      std::vector<ftype> a;
      a = v->induced_voltage_generation(length);
      //std::cout << "a size is " << a.size() << '\n';
      if (length > 0) {
         extIndVolt.resize(a.size(), 0);
         std::transform(extIndVolt.begin(), extIndVolt.end(),
                        a.begin(), extIndVolt.begin(),
                        std::plus<ftype>());
      }
      tempIndVolt.resize(v->fInducedVoltage.size(), 0);
      std::transform(tempIndVolt.begin(), tempIndVolt.end(),
                     v->fInducedVoltage.begin(), tempIndVolt.begin(),
                     std::plus<ftype>());
   }

   fInducedVoltage = tempIndVolt;
   //std::cout << "extIndVolt size = " << extIndVolt.size() << std::endl;
   return extIndVolt;


}

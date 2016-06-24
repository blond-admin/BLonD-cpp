/*
* @Author: Konstantinos Iliakis
* @Date:   2016-05-04 11:51:39
* @Last Modified by:   Konstantinos Iliakis
* @Last Modified time: 2016-05-04 14:38:41
*/

#ifndef IMPEDANCES_INDUCEDVOLTAGE_H_
#define IMPEDANCES_INDUCEDVOLTAGE_H_

//class InducedVoltage;

#include "configuration.h"
//#include <complex>
#include <vector>
#include "Intensity.h"
#include "fft.h"

//typedef std::complex<float> complex_t;



enum time_or_freq {
   time_domain, freq_domain
};

typedef enum freq_res_option_t {
   round_option, ceil_option, floor_option
} freq_res_option_t;

//uint next_regular(uint target);



class InducedVoltage {
public:
   std::vector<ftype> fInducedVoltage;

   InducedVoltage() {};
   inline void linear_interp_kick(const ftype *__restrict__ beam_dt,
                                  ftype *__restrict__ beam_dE,
                                  const ftype *__restrict__ voltage_array,
                                  const ftype *__restrict__ bin_centers,
                                  const int n_slices,
                                  const int n_macroparticles,
                                  const ftype acc_kick = 0.0);
   virtual void track() = 0;
   virtual void reprocess() = 0;
   virtual std::vector<ftype> induced_voltage_generation(uint length = 0) = 0;
   virtual ~InducedVoltage() {};
};



class InducedVoltageTime: public InducedVoltage {
public:


   std::vector<Intensity *> fWakeSourceList;
   std::vector<ftype> fTimeArray;
   std::vector<ftype> fTotalWake;
   uint fCut;
   uint fShape;
   time_or_freq fTimeOrFreq;


   void track();
   void sum_wakes(std::vector<ftype> &v);
   void reprocess();
   std::vector<ftype> induced_voltage_generation(uint length = 0);
   InducedVoltageTime(std::vector<Intensity *> &WakeSourceList,
                      time_or_freq TimeOrFreq = freq_domain);
   ~InducedVoltageTime()
   {
      fft::destroy_plans();
   };
};


class InducedVoltageFreq: public InducedVoltage {
public:

   // Impedance sources inputed as a list (eg: list of BBResonators objects)*
   std::vector<Intensity *> fImpedanceSourceList;

   // *Input frequency resolution in [Hz], the beam profile sampling for the spectrum
   // will be adapted according to the freq_res_option.*
   ftype fFreqResolutionInput;

   // Number of turns to be considered as memory for induced voltage calculation.*
   uint fNTurnsMem;
   bool fRecalculationImpedance;
   bool fSaveIndividualVoltages;
   // *Real frequency resolution in [Hz], according to the obtained n_fft_sampling.*
   ftype fFreqResolution;
   // *Frequency array of the impedance in [Hz]*
   f_vector_t fFreqArray;
   uint fNFFTSampling;
   freq_res_option_t fFreqResOption;
   // *Total impedance array of all sources in* [:math:`\Omega`]
   complex_vector_t fTotalImpedance;
   complex_vector_t fMatrixSaveIndividualImpedances;
   f_vector_t fMatrixSaveIndividualVoltages;

   uint fLenArrayMem;
   uint fLenArrayMemExt;
   uint fNPointsFFT;

   f_vector_t fFreqArrayMem;
   complex_vector_t fTotalImpedanceMem;
   f_vector_t fTimeArrayMem;


   // *Induced voltage from the sum of the wake sources in [V]*
   //f_vector_t fInducedVoltage;

   // and many more! //

   void track();
   void sum_impedances(f_vector_t &);

   //Reprocess the impedance contributions with respect to the new_slicing.
   void reprocess();
   std::vector<ftype> induced_voltage_generation(uint length = 0);
   InducedVoltageFreq(std::vector<Intensity *> &impedanceSourceList,
                      ftype freqResolutionInput = 0.0,
                      freq_res_option_t freq_res_option = freq_res_option_t::round_option,
                      uint NTurnsMem = 0,
                      bool recalculationImpedance = false,
                      bool saveIndividualVoltages = false);
   ~InducedVoltageFreq()
   {
      fft::destroy_plans();
   };
};


class TotalInducedVoltage : public InducedVoltage {
public:
   std::vector<InducedVoltage *> fInducedVoltageList;
   std::vector<ftype> fTimeArray;
   std::vector<ftype> fRevTimeArray;
   uint fCounterTurn = 0;
   uint fNTurnsMemory;
   bool fInductiveImpedanceOn = false;

   void track();
   void track_memory();
   void track_ghosts_particles();
   std::vector<ftype> induced_voltage_sum(uint length = 0);
   void reprocess();
   
   std::vector<ftype> induced_voltage_generation(uint length = 0)
   {
      return std::vector<ftype>();
   };

   TotalInducedVoltage(std::vector<InducedVoltage *> &InducedVoltageList,
                       uint NTurnsMemory = 0,
                       std::vector<ftype> RevTimeArray = std::vector<ftype>());
   
   ~TotalInducedVoltage()
   {
      fft::destroy_plans();
   };

};

#endif /* IMPEDANCES_INDUCEDVOLTAGE_H_ */
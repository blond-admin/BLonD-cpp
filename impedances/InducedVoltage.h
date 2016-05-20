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

//typedef std::complex<float> complex_t;



enum time_or_freq {
   time_domain, freq_domain
};
//unsigned int next_regular(unsigned int target);

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
   virtual std::vector<ftype> induced_voltage_generation(unsigned int length = 0) = 0;
   virtual ~InducedVoltage() {};
};



class InducedVoltageTime: public InducedVoltage {
public:


   std::vector<Intensity *> fWakeSourceList;
   std::vector<ftype> fTimeArray;
   std::vector<ftype> fTotalWake;
   unsigned int fCut;
   unsigned int fShape;
   time_or_freq fTimeOrFreq;


   void track();
   void sum_wakes(std::vector<ftype> &v);
   void reprocess();
   std::vector<ftype> induced_voltage_generation(unsigned int length = 0);
   InducedVoltageTime(std::vector<Intensity *> &WakeSourceList,
                      time_or_freq TimeOrFreq = freq_domain);
   ~InducedVoltageTime() {};
};


class InducedVoltageFreq: public InducedVoltage {
public:
   void track();
   void sum_impedances();
   void reprocess();
   std::vector<ftype> induced_voltage_generation(unsigned int length = 0);
   InducedVoltageFreq();
   ~InducedVoltageFreq() {};
};

#endif /* IMPEDANCES_INDUCEDVOLTAGE_H_ */